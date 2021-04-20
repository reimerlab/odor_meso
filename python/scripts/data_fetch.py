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

# # Fetch data from AT database for ingestion on `jr-compute001`

import datajoint as dj
import pickle
import lzma
import os

dj.config['database.prefix'] = ''

# +
# animal_id  = 1571
# session    = 1
# scan_idx   = 1

# animal_id  = 121
# session    = 7
# scan_idx   = 1

animal_id  = 122          # Animal with multiple z planes to be stitched
session    = 4
scan_idx   = 5

# animal_id  = 14864
# segmentation_method = 2

# +
directory = f'{os.environ.get("INGESTION_STORAGE")}/animal_id_{animal_id}/'
if not os.path.exists(directory):
    os.mkdir(directory)

def save_data(filename, data):
     with lzma.open(directory + filename, "wb") as handle:
          pickle.dump(data, handle)


# -

# Instantiate modules
mice       = dj.create_virtual_module('mice', 'common_mice')
shared     = dj.create_virtual_module('shared', 'pipeline_shared')
experiment = dj.create_virtual_module('experiment', 'pipeline_experiment')
odor       = dj.create_virtual_module('odor', 'pipeline_odor')
meso       = dj.create_virtual_module('meso', 'pipeline_meso')
pupil      = dj.create_virtual_module('eye', 'pipeline_eye')

# # Fetch data

mice_Mice                    = (mice.Mice & f'animal_id={animal_id}').fetch(as_dict=True)

shared_pipelineVersion       = shared.PipelineVersion.fetch(as_dict=True)

experiment_Rig               = experiment.Rig.fetch(as_dict=True)
experiment_Lens              = experiment.Lens.fetch(as_dict=True)
experiment_Aim               = experiment.Aim.fetch(as_dict=True)
experiment_Software          = experiment.Software.fetch(as_dict=True)
experiment_Person            = experiment.Person.fetch(as_dict=True)
experiment_Anesthesia        = experiment.Anesthesia.fetch(as_dict=True)
experiment_BrainArea         = experiment.BrainArea.fetch(as_dict=True)
experiment_Session           = (experiment.Session & f'animal_id={animal_id}').fetch(as_dict=True)
experiment_Scan              = (experiment.Scan & f'animal_id={animal_id}').fetch(as_dict=True)
experiment_Scan              = (experiment.Scan & f'animal_id={animal_id}' & f'scan_idx={scan_idx}').fetch(as_dict=True)
experiment_Scan_EyeVideo     = (experiment.Scan.EyeVideo & f'animal_id={animal_id}' & f'scan_idx={scan_idx}').fetch(as_dict=True)
experiment_Scan_BehaviorFile = (experiment.Scan.BehaviorFile & f'animal_id={animal_id}' & f'scan_idx={scan_idx}').fetch(as_dict=True)
experiment_Scan_Laser        = (experiment.Scan.Laser & f'animal_id={animal_id}' & f'scan_idx={scan_idx}').fetch(as_dict=True)
experiment_TreadmillSpecs    = experiment.TreadmillSpecs.fetch(as_dict=True)

odor_Odorant                 = odor.Odorant.fetch(as_dict=True)
odor_OdorSolution            = odor.OdorSolution.fetch(as_dict=True)
odor_OdorSession             = (odor.OdorSession & f'animal_id={animal_id}').fetch(as_dict=True)
odor_OdorConfig              = (odor.OdorConfig & f'animal_id={animal_id}').fetch(as_dict=True)
odor_OdorRecording           = (odor.OdorRecording & f'animal_id={animal_id}').fetch(as_dict=True)
odor_MesoMatch               = (odor.MesoMatch & f'animal_id={animal_id}').fetch(as_dict=True)

meso_Version                 = meso.Version.fetch(as_dict=True)
meso_ScanInfo                = (meso.ScanInfo & f'animal_id={animal_id}').fetch(as_dict=True)
meso_CorrectionChannel       = (meso.CorrectionChannel & f'animal_id={animal_id}' & f'scan_idx={scan_idx}').fetch(as_dict=True)
meso_SegmentationTask        = (meso.SegmentationTask & f'animal_id={animal_id}' & f'scan_idx={scan_idx}').fetch(as_dict=True)
meso_Segmentation            = (meso.Segmentation & f'animal_id={animal_id}' & f'scan_idx={scan_idx}').fetch(as_dict=True)
meso_Segmentation_Manual     = (meso.Segmentation.Manual & f'animal_id={animal_id}' & f'scan_idx={scan_idx}').fetch(as_dict=True)
meso_Segmentation_Mask       = (meso.Segmentation.Mask & f'animal_id={animal_id}' & f'scan_idx={scan_idx}').fetch(as_dict=True)

# +
# WIP
# pupil_TrackingTask                      = (pupil.TrackingTask & f'animal_id={animal_id}' & f'scan_idx={scan_idx}').fetch(as_dict=True)
# pupil_TrackingTask_ManualParameters     = (pupil.TrackingTask.ManualParameters & f'animal_id={animal_id}' & f'scan_idx={scan_idx}').fetch(as_dict=True)
# pupil_TrackingTask_Ignore               = (pupil.TrackingTask.Ignore & f'animal_id={animal_id}' & f'scan_idx={scan_idx}').fetch(as_dict=True)
# pupil_TrackingTask_Mask                 = (pupil.TrackingTask.Mask & f'animal_id={animal_id}' & f'scan_idx={scan_idx}').fetch(as_dict=True)
# pupil_ManuallyTrackedContours           = (pupil.ManuallyTrackedContours & f'animal_id={animal_id}' & f'scan_idx={scan_idx}').fetch(as_dict=True)
# pupil_ManuallyTrackedContours_Frame     = (pupil.ManuallyTrackedContours.Frame & f'animal_id={animal_id}' & f'scan_idx={scan_idx}').fetch(as_dict=True)
# pupil_ManuallyTrackedContours_Parameter = (pupil.ManuallyTrackedContours.Parameter & f'animal_id={animal_id}' & f'scan_idx={scan_idx}').fetch(as_dict=True)
# pupil_ConfigDeeplabcut                  = pupil.ConfigDeeplabcut.fetch(as_dict=True)
# -

# # Save data

save_data('mice.Mice.xz', mice_Mice)

save_data('shared.PipelineVersion.xz', shared_pipelineVersion)

save_data('experiment.Rig.xz', experiment_Rig)
save_data('experiment.Lens.xz', experiment_Lens)
save_data('experiment.Aim.xz', experiment_Aim)
save_data('experiment.Software.xz', experiment_Software)
save_data('experiment.Person.xz', experiment_Person)
save_data('experiment.Anesthesia.xz', experiment_Anesthesia)
save_data('experiment.BrainArea.xz', experiment_BrainArea)
save_data('experiment.Session.xz', experiment_Session)
save_data('experiment.Scan.xz', experiment_Scan)
save_data('experiment.Scan.EyeVideo.xz', experiment_Scan_EyeVideo)
save_data('experiment.Scan.BehaviorFile.xz', experiment_Scan_BehaviorFile)
save_data('experiment.Scan.Laser.xz', experiment_Scan_Laser)
save_data('experiment.TreadmillSpecs.xz', experiment_TreadmillSpecs)

save_data('odor.Odorant.xz', odor_Odorant)
save_data('odor.OdorSolution.xz', odor_OdorSolution)
save_data('odor.OdorSession.xz', odor_OdorSession)
save_data('odor.OdorConfig.xz', odor_OdorConfig)
save_data('odor.OdorRecording.xz', odor_OdorRecording)
save_data('odor.MesoMatch.xz', odor_MesoMatch)

save_data('meso.Version.xz', meso_Version)
save_data('meso.ScanInfo.xz', meso_ScanInfo)
save_data('meso.CorrectionChannel.xz', meso_CorrectionChannel)
save_data('meso.SegmentationTask.xz', meso_SegmentationTask)
save_data('meso.Segmentation.xz', meso_Segmentation)
save_data('meso.Segmentation.Manual.xz', meso_Segmentation_Manual)
save_data('meso.Segmentation.Mask.xz', meso_Segmentation_Mask)

# +
# WIP
# save_data('pupil.TrackingTask.xz', pupil_TrackingTask)
# save_data('pupil.TrackingTask.ManualParameters.xz', pupil_TrackingTask_ManualParameters)
# save_data('pupil.TrackingTask.Ignore.xz', pupil_TrackingTask_Ignore)
# save_data('pupil.TrackingTask.Mask.xz', pupil_TrackingTask_Mask)
# save_data('pupil.ManuallyTrackedContours.xz', pupil_ManuallyTrackedContours)
# save_data('pupil.ManuallyTrackedContours.Frame.xz', pupil_ManuallyTrackedContours_Frame)
# save_data('pupil.ManuallyTrackedContours.Parameter.xz', pupil_ManuallyTrackedContours_Parameter)
# save_data('pupil.ConfigDeeplabcut.xz', pupil_ConfigDeeplabcut)
# -

#
