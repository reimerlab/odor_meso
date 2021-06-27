import datajoint as dj

import numpy as np
import os

from . import experiment, lab
from .utils import h5
from .exceptions import PipelineException

dj.config['database.prefix'] = os.environ.get('DJ_PREFIX', '')
schema = dj.schema(dj.config['database.prefix'] + 'pipeline_stimulus', locals())

@schema
class BehaviorSync(dj.Computed):
     definition = """ # syncing scanimage frame times to behavior clock
     -> experiment.Scan
     ---
     frame_times                         : longblob      # start of each scanimage slice in behavior clock
     behavior_sync_ts=CURRENT_TIMESTAMP  : timestamp
     """
     @property
     def key_source(self):
          return experiment.Scan() & experiment.Scan.BehaviorFile().proj()

     def _make_tuples(self, key):
          # Get behavior filename
          print(key)
          stim_type = (experiment.ScanStimType() & key).fetch1('stim_type_id')
          behavior_path = (experiment.Session() & key).fetch1('behavior_path')
          local_path = lab.Paths().get_local_path(behavior_path)
          filename = (experiment.Scan.BehaviorFile() & key).fetch1('filename')
          full_filename = os.path.join(local_path, filename)

          # Read file
          data = h5.read_behavior_file(full_filename)

          # Get counter timestamps and convert to seconds
          timestamps_in_secs = h5.ts2sec(data['ts'], is_packeted=True)

          # Detect rising edges in scanimage clock signal (start of each frame)
          binarized_signal = data['scanImage'] > 2.7 # TTL voltage low/high threshold
          rising_edges = np.where(np.diff(binarized_signal.astype(int)) > 0)[0]
          frame_times = timestamps_in_secs[rising_edges]
          
          ## done for testing purposes 
          if stim_type != 1:
               frame_times = frame_times[10:]

          # Correct NaN gaps in timestamps (mistimed or dropped packets during recording)
          if np.any(np.isnan(frame_times)):
               # Raise exception if first or last frame pulse was recorded in mistimed packet
               if np.isnan(frame_times[0]) or np.isnan(frame_times[-1]):
                    msg = ('First or last frame happened during misstamped packets. Pulses '
                         'could have been missed: start/end of scanning is unknown.')
                    raise PipelineException(msg)

               # Fill each gap of nan values with correct number of timepoints
               frame_period = np.nanmedian(np.diff(frame_times)) # approx
               nan_limits = np.where(np.diff(np.isnan(frame_times)))[0]
               nan_limits[1::2] += 1 # limits are indices of the last valid point before the nan gap and first after it
               correct_fts = []
               for i, (start, stop) in enumerate(zip(nan_limits[::2], nan_limits[1::2])):
                    correct_fts.extend(frame_times[0 if i == 0 else nan_limits[2 * i - 1]: start + 1])
                    num_missing_points = int(round((frame_times[stop] - frame_times[start]) /
                                                  frame_period - 1))
                    correct_fts.extend(np.linspace(frame_times[start], frame_times[stop],
                                                  num_missing_points + 2)[1:-1])
               correct_fts.extend(frame_times[nan_limits[-1]:])
               frame_times = correct_fts

          # Check that frame times occur at the same period
          frame_intervals = np.diff(frame_times)
          frame_period = np.median(frame_intervals)
          if np.any(abs(frame_intervals - frame_period) > 0.15 * frame_period):
               raise PipelineException('Frame time period is irregular')

          # Drop last frame time if scan crashed or was stopped before completion
          valid_times = ~np.isnan(timestamps_in_secs[rising_edges[0]: rising_edges[-1]]) # restricted to scan period
          binarized_valid = binarized_signal[rising_edges[0]: rising_edges[-1]][valid_times]
          frame_duration = np.mean(binarized_valid) * frame_period
          falling_edges = np.where(np.diff(binarized_signal.astype(int)) < 0)[0]
          last_frame_duration = timestamps_in_secs[falling_edges[-1]] - frame_times[-1]
          if (np.isnan(last_frame_duration) or last_frame_duration < 0 or
               abs(last_frame_duration - frame_duration) > 0.15 * frame_duration):
               frame_times = frame_times[:-1]

          self.insert1({**key, 'frame_times': frame_times})