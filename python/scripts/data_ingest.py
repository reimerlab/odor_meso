import os
import lzma
import pickle
from tqdm import tqdm
import datajoint as dj
from pipeline import meso, stack, treadmill, pupil, odor, mice, experiment, shared

directory = '/data/odor_meso/' + os.environ.get('INGESTION_STORAGE')

animal_id = 1571

#%% Function to insert metadata
def insert_data(filename, table, skip_duplicates, ignore_extra_fields, allow_direct_insert, print_data=False):
    with lzma.open(filename, "rb") as f:
        loaded_data = pickle.load(f)
        f.close()
    
    if print_data:
        print(loaded_data)

    for i in tqdm(range(len(loaded_data))):
        table.insert1(loaded_data[i],
                      skip_duplicates=skip_duplicates,
                      ignore_extra_fields=ignore_extra_fields,
                      allow_direct_insert=allow_direct_insert)

#%% Insert or populate for each schema
insert_data(directory + f'/animal_id_{animal_id}' + '/mice.Mice.xz', mice.Mice, False, False, False, True)

insert_data(directory + f'/animal_id_{animal_id}' + '/shared.PipelineVersion.xz', shared.PipelineVersion, True, False, False, True)

insert_data(directory + f'/animal_id_{animal_id}' + '/experiment.Rig.xz', experiment.Rig, True, False, False, True)
insert_data(directory + f'/animal_id_{animal_id}' + '/experiment.Lens.xz', experiment.Lens, True, False, False, True)
insert_data(directory + f'/animal_id_{animal_id}' + '/experiment.Aim.xz', experiment.Aim, True, False, False, True)
insert_data(directory + f'/animal_id_{animal_id}' + '/experiment.Software.xz', experiment.Software, True, False, False, True)
insert_data(directory + f'/animal_id_{animal_id}' + '/experiment.Session.xz', experiment.Session, False, True, False, True)
insert_data(directory + f'/animal_id_{animal_id}' + '/experiment.Scan.xz', experiment.Scan, False, False, False, True)
insert_data(directory + f'/animal_id_{animal_id}' + '/experiment.Scan.EyeVideo.xz', experiment.Scan.EyeVideo, False, False, False, True)
insert_data(directory + f'/animal_id_{animal_id}' + '/experiment.Scan.BehaviorFile.xz', experiment.Scan.BehaviorFile, False, False, False, True)
insert_data(directory + f'/animal_id_{animal_id}' + '/experiment.Scan.Laser.xz', experiment.Scan.Laser, False, False, False, True)
insert_data(directory + f'/animal_id_{animal_id}' + '/experiment.TreadmillSpecs.xz', experiment.TreadmillSpecs, False, False, False, True)

experiment.ExperimentalIdentifier.insert(experiment.Session.fetch('KEY'),skip_duplicates=True)

insert_data(directory + f'/animal_id_{animal_id}' + '/odor.Odorant.xz', odor.Odorant, False, False, False, True)
insert_data(directory + f'/animal_id_{animal_id}' + '/odor.OdorSolution.xz', odor.OdorSolution, False, False, False, True)
insert_data(directory + f'/animal_id_{animal_id}' + '/odor.OdorSession.xz', odor.OdorSession, False, False, False, True)
insert_data(directory + f'/animal_id_{animal_id}' + '/odor.OdorConfig.xz', odor.OdorConfig, False, False, False, True)
insert_data(directory + f'/animal_id_{animal_id}' + '/odor.OdorRecording.xz', odor.OdorRecording, False, False, False, True)
insert_data(directory + f'/animal_id_{animal_id}' + '/odor.MesoMatch.xz', odor.MesoMatch, False, False, False, True)
odor.OdorTrials.populate()
odor.OdorSync.populate()
odor.Respiration.populate()

insert_data(directory + f'/animal_id_{animal_id}' + '/meso.Version.xz', meso.Version, False, False, False, True)
meso.ScanInfo.populate()
dj.config["enable_python_native_blobs"] = True
meso.Quality.populate()
insert_data(directory + f'/animal_id_{animal_id}' + '/meso.CorrectionChannel.xz', meso.CorrectionChannel, False, False, False, True)
meso.RasterCorrection.populate()
meso.MotionCorrection.populate()
meso.SummaryImages.populate()
insert_data(directory + f'/animal_id_{animal_id}' + '/meso.SegmentationTask.xz', meso.SegmentationTask, False, False, False, True)
insert_data(directory + f'/animal_id_{animal_id}' + '/meso.Segmentation.xz', meso.Segmentation, False, False, True, False)
insert_data(directory + f'/animal_id_{animal_id}' + '/meso.Segmentation.Manual.xz', meso.Segmentation.Manual, False, False, True, False)
insert_data(directory + f'/animal_id_{animal_id}' + '/meso.Segmentation.Mask.xz', meso.Segmentation.Mask, False, False, True, False)
meso.Fluorescence.populate()
meso.MaskClassification.populate()
meso.ScanSet.populate()
meso.Activity.populate()
meso.ScanDone.populate()

treadmill.Treadmill.populate()

odor.OdorAnalysis.populate()

#%% Fetch data
data = (treadmill.Treadmill() & 'animal_id=1571').fetch('treadmill_vel')

