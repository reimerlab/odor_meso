# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

import os
import lzma
import pickle
from tqdm import tqdm
import datajoint as dj

dj.__version__

# # Define database prefix


from pipeline import meso, stack, treadmill, odor, mice, experiment, shared#, pupil - version of deeplabcut doesn't work with upgraded version of scipy, which was upgraded for griddata interpolation of masks in stitched space

# # Define animal_id data to be ingested

# animal_id = 1571
# animal_id  = 14864
animal_id = 122

# # Function to insert metadata

# +
directory = f'{os.environ.get("INGESTION_STORAGE")}/animal_id_{animal_id}/'

def insert_data(filename, table, skip_duplicates, ignore_extra_fields, allow_direct_insert):
    with lzma.open(directory + filename, "rb") as f:
        loaded_data = pickle.load(f)
        f.close()

#     print(loaded_data)

    for i in tqdm(range(len(loaded_data))):
        table.insert1(loaded_data[i],
                      skip_duplicates=skip_duplicates,
                      ignore_extra_fields=ignore_extra_fields,
                      allow_direct_insert=allow_direct_insert)


# -

# # Mice pipeline

insert_data(filename='mice.Mice.xz', 
            table=mice.Mice, 
            skip_duplicates=False, 
            ignore_extra_fields=False, 
            allow_direct_insert=False)

# # Shared pipeline

insert_data(filename='shared.PipelineVersion.xz', 
            table=shared.PipelineVersion, 
            skip_duplicates=True, 
            ignore_extra_fields=False, 
            allow_direct_insert=False)

# # Experiment pipeline

insert_data('experiment.Rig.xz', experiment.Rig, True, False, False)
insert_data('experiment.Lens.xz', experiment.Lens, True, False, False)
insert_data('experiment.Aim.xz', experiment.Aim, True, False, False)
insert_data('experiment.Software.xz', experiment.Software, True, False, False)
insert_data('experiment.Session.xz', experiment.Session, False, True, False)
insert_data('experiment.Scan.xz', experiment.Scan, False, False, False)
insert_data('experiment.Scan.EyeVideo.xz', experiment.Scan.EyeVideo, False, False, False)
insert_data('experiment.Scan.BehaviorFile.xz', experiment.Scan.BehaviorFile, False, False, False)
insert_data('experiment.Scan.Laser.xz', experiment.Scan.Laser, False, False, False)
insert_data('experiment.TreadmillSpecs.xz', experiment.TreadmillSpecs, False, False, False)

experiment.ExperimentalIdentifier.insert(experiment.Scan.fetch('KEY'),skip_duplicates=True)

# # Meso pipeline

insert_data(filename='meso.Version.xz', 
            table=meso.Version, 
            skip_duplicates=False, 
            ignore_extra_fields=False, 
            allow_direct_insert=False)

meso.ScanInfo.populate()

meso.Quality.populate()

# + jupyter={"outputs_hidden": true}
# insert_data('meso.CorrectionChannel.xz', meso.CorrectionChannel, False, False, False)

# + jupyter={"outputs_hidden": true}
meso.RasterCorrection.populate()
# -

meso.MotionCorrection.populate()

meso.SummaryImages.populate()

meso.Stitch.populate()

insert_data('meso.SegmentationTask.xz', meso.SegmentationTask, False, False, False)
insert_data('meso.Segmentation.xz', meso.Segmentation, False, False, True)
insert_data('meso.Segmentation.Manual.xz', meso.Segmentation.Manual, False, False, True)
insert_data('meso.Segmentation.Mask.xz', meso.Segmentation.Mask, False, False, True)

meso.Fluorescence.populate()

meso.MaskClassification.populate()

meso.ScanSet.populate()

meso.Activity.populate()

meso.ScanDone.populate()

# # WIP - Automated glomeruli segmentation

# +
entry_segmentation_task=[]

for field in [1,2,3]:
    entry_segmentation_task.append({'animal_id': animal_id,
      'session': 1,
      'scan_idx': 1,
      'field': field,
      'channel': 1,
      'segmentation_method': 2,
      'compartment': 'unknown'})
meso.SegmentationTask.insert(entry_segmentation_task)  
# -

meso.Segmentation.populate()

# +
# (meso.SegmentationTask & 'segmentation_method=2').delete()
# (meso.Segmentation & 'segmentation_method=2').delete()
# -

# # Odor pipeline

insert_data('odor.Odorant.xz', odor.Odorant, False, False, False)
insert_data('odor.OdorSolution.xz', odor.OdorSolution, False, False, False)
insert_data('odor.OdorSession.xz', odor.OdorSession, False, False, False)
insert_data('odor.OdorConfig.xz', odor.OdorConfig, False, False, False)
insert_data('odor.OdorRecording.xz', odor.OdorRecording, False, False, False)
insert_data('odor.MesoMatch.xz', odor.MesoMatch, False, False, False)

odor.OdorTrials.populate()

odor.OdorSync.populate()

odor.Respiration.populate()

odor.OdorAnalysis.populate()

# # Treadmill pipeline

treadmill.Treadmill.populate()

# # TODO: Stack pipeline

# + jupyter={"outputs_hidden": true}
experiment.Stack.insert1([121,8,1,'meso','unset','stack','scanimage','2017b',0,0,0,'','2020-12-09 11:46:48'])#unsure of depths
experiment.Stack.Filename.insert1([121,8,1,1,'121_8_00001',0])#{animal_id, session, stack_idx, filename_idx, filename, surf_depth}) #unsure of depths

# + jupyter={"outputs_hidden": true}
stack.StackInfo.populate()
stack.Quality.populate()
# stack.CorrectionChannel.insert1([121,8,1,1])#correct?
stack.RasterCorrection.populate()
stack.MotionCorrection.populate()
stack.Stitching.populate()
stack.CorrectedStack.populate()
stack.PreprocessedStack.populate()#dj.conn.ping required
# stack.Surface.populate() # did not populate
# stack.SegmentationTask.insert1([121,8,1,1,1,2,'soma'])#correct?

# + jupyter={"outputs_hidden": true}
stack.Segmentation.populate()
# stack.RegistrationTask.insert1([121,8,1,1,,5])#correct?
stack.Registration.populate()
stack.FieldSegmentation.populate()
stack.RegistrationOverTime.populate()
stack.Drift.populate()
stack.StackSet.populate()
stack.Area.populate()
# -

# # TODO: Pupil pipeline

# + jupyter={"outputs_hidden": true}
pupil.Eye.populate()


# + jupyter={"outputs_hidden": true}
class TrackingTask(dj.Manual):
pupil.TrackedVideo.populate()
class ManuallyTrackedContours(dj.Manual, AutoPopulate):
pupil.FittedContour.populate()
class ConfigDeeplabcut(dj.Manual):
pupil.Tracking.populate()
pupil.FittedPupil.populate()
# -

# # Drop schemas

# + jupyter={"outputs_hidden": true}
from pipeline import meso, reso, stack, treadmill, odor, mice, \
                     experiment, shared, pupil, map, injection, \
                     anatomy, lab, notify

# + jupyter={"outputs_hidden": true}
notify.schema.drop()
pupil.schema.drop()
stack.schema.drop()
treadmill.schema.drop()
odor.schema.drop()
meso.schema.drop()
reso.schema.drop()
injection.schema.drop()
anatomy.schema.drop()
map.schema.drop()
experiment.schema.drop()
shared.schema.drop()
mice.schema.drop()
lab.schema.drop()

# + jupyter={"outputs_hidden": true}

