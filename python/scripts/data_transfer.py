import os
import datajoint as dj
import odor_meso.python.scripts.data_copy as data_copy

dj.config['database.prefix'] = os.environ.get('DJ_PREFIX', 'user_kabilar_')

# ------------------------------------------------------------------------------
# Instantiate virtual module connections to source (at-database) and destination (jr-database) databases
pc = data_copy.PipelineCopy.get()

# Instantiate modules on jr-database for populating tables
from pipeline import meso, odor, treadmill, stack

# Verify connection to jr-database
assert meso.schema.connection.conn_info['host'] == pc.dst_vmods['meso'].schema.connection.conn_info['host']

# Create restriction for database transfer
query_meso = pc.src_vmods['odor'].MesoMatch()
query_experiment = pc.src_vmods['experiment'].AutoProcessing.proj()

populate_settings = {'display_progress': True, 'reserve_jobs': True, 'suppress_errors': False}

# ------------------------------------------------------------------------------
print('Mice pipeline')

mice_Mice = (pc.src_vmods['mice'].Mice & query_meso).fetch(as_dict=True)
for i in range(len(mice_Mice)): 
     if mice_Mice[i]['dob'] == '0000-00-00': mice_Mice[i]['dob'] = None
     if mice_Mice[i]['dow'] == '0000-00-00': mice_Mice[i]['dow'] = None
pc.dst_vmods['mice'].Mice.insert(mice_Mice, skip_duplicates=False)

# ------------------------------------------------------------------------------
print('Experiment pipeline')

# Lookup tables
experiment_Rig = pc.src_vmods['experiment'].Rig.fetch(as_dict=True)
pc.dst_vmods['experiment'].Rig.insert(experiment_Rig, skip_duplicates=True)

experiment_Lens = pc.src_vmods['experiment'].Lens.fetch(as_dict=True)
pc.dst_vmods['experiment'].Lens.insert(experiment_Lens, skip_duplicates=True)

experiment_Aim = pc.src_vmods['experiment'].Aim.fetch(as_dict=True)
pc.dst_vmods['experiment'].Aim.insert(experiment_Aim, skip_duplicates=True)

experiment_Software = pc.src_vmods['experiment'].Software.fetch(as_dict=True)
pc.dst_vmods['experiment'].Software.insert(experiment_Software, skip_duplicates=True)

experiment_Person = pc.src_vmods['experiment'].Person.fetch(as_dict=True)
pc.dst_vmods['experiment'].Person.insert(experiment_Person, skip_duplicates=True)

experiment_BrainArea = pc.src_vmods['experiment'].BrainArea.fetch(as_dict=True)
pc.dst_vmods['experiment'].BrainArea.insert(experiment_BrainArea, skip_duplicates=True)

experiment_TreadmillSpecs = pc.src_vmods['experiment'].TreadmillSpecs.fetch(as_dict=True)
pc.dst_vmods['experiment'].TreadmillSpecs.insert(experiment_TreadmillSpecs, skip_duplicates=False)

# Manual tables
experiment_Session = (pc.src_vmods['experiment'].Session & query_meso).fetch(as_dict=True)
pc.dst_vmods['experiment'].Session.insert(experiment_Session, skip_duplicates=True, ignore_extra_fields=True)

experiment_Session_Fluorophore = (pc.src_vmods['experiment'].Session.Fluorophore & query_meso).fetch(as_dict=True)
pc.dst_vmods['experiment'].Session.Fluorophore.insert(experiment_Session_Fluorophore, skip_duplicates=True, ignore_extra_fields=True)

experiment_Session_TargetStructure = (pc.src_vmods['experiment'].Session.TargetStructure & query_meso).fetch(as_dict=True)
pc.dst_vmods['experiment'].Session.TargetStructure.insert(experiment_Session_TargetStructure, skip_duplicates=True, ignore_extra_fields=True)

experiment_Session_PMTFilterSet = (pc.src_vmods['experiment'].Session.PMTFilterSet & query_meso).fetch(as_dict=True)
pc.dst_vmods['experiment'].Session.PMTFilterSet.insert(experiment_Session_PMTFilterSet, skip_duplicates=True, ignore_extra_fields=True)

experiment_Scan = (pc.src_vmods['experiment'].Scan & query_meso).fetch(as_dict=True)
pc.dst_vmods['experiment'].Scan.insert(experiment_Scan, skip_duplicates=True)

experiment_Scan_EyeVideo = (pc.src_vmods['experiment'].Scan.EyeVideo & query_meso).fetch(as_dict=True)
pc.dst_vmods['experiment'].Scan.EyeVideo.insert(experiment_Scan_EyeVideo, skip_duplicates=True)

experiment_Scan_BehaviorFile = (pc.src_vmods['experiment'].Scan.BehaviorFile & query_meso).fetch(as_dict=True)
pc.dst_vmods['experiment'].Scan.BehaviorFile.insert(experiment_Scan_BehaviorFile, skip_duplicates=True)

experiment_Scan_Laser = (pc.src_vmods['experiment'].Scan.Laser & query_meso).fetch(as_dict=True)
pc.dst_vmods['experiment'].Scan.Laser.insert(experiment_Scan_Laser, skip_duplicates=True)

experiment_Scan = (pc.src_vmods['experiment'].Scan & (query_experiment & query_meso)).fetch('KEY')
pc.dst_vmods['experiment'].ExperimentalIdentifier.insert(experiment_Scan, skip_duplicates=True)

experiment_AutoProcessing = (pc.src_vmods['experiment'].AutoProcessing & query_meso).fetch('KEY')
pc.dst_vmods['experiment'].AutoProcessing.insert(experiment_AutoProcessing, skip_duplicates=True)

# ------------------------------------------------------------------------------
print('Meso pipeline')

meso_Version = pc.src_vmods['meso'].Version.fetch(as_dict=True)
pc.dst_vmods['meso'].Version.insert(meso_Version, skip_duplicates=True)

meso.ScanInfo.populate((query_experiment & query_meso), **populate_settings)

meso.Quality.populate((query_experiment & query_meso), **populate_settings)

meso_CorrectionChannel = (pc.src_vmods['meso'].CorrectionChannel & (query_experiment & query_meso)).fetch(as_dict=True)
pc.dst_vmods['meso'].CorrectionChannel.insert(meso_CorrectionChannel, skip_duplicates=True)

meso.RasterCorrection.populate((query_experiment & query_meso), **populate_settings)

meso.MotionCorrection.populate((query_experiment & query_meso), **populate_settings)

meso.SummaryImages.populate((query_experiment & query_meso), **populate_settings)

meso.Stitch.populate((query_experiment & query_meso), **populate_settings)

# Glomeruli manual and CNMF segmentations from at-database
meso_SegmentationTask = (pc.src_vmods['meso'].SegmentationTask & (query_experiment & query_meso)).fetch(as_dict=True)
pc.dst_vmods['meso'].SegmentationTask.insert(meso_SegmentationTask, skip_duplicates=True)

meso_Segmentation = (pc.src_vmods['meso'].Segmentation & (query_experiment & query_meso)).fetch(as_dict=True)
pc.dst_vmods['meso'].Segmentation.insert(meso_Segmentation, skip_duplicates=True, allow_direct_insert=True)

meso_Segmentation_Manual = (pc.src_vmods['meso'].Segmentation.Manual & (query_experiment & query_meso)).fetch(as_dict=True)
pc.dst_vmods['meso'].Segmentation.Manual.insert(meso_Segmentation_Manual, skip_duplicates=True, allow_direct_insert=True)

meso_Segmentation_CNMF = (pc.src_vmods['meso'].Segmentation.CNMF & (query_experiment & query_meso)).fetch(as_dict=True)
pc.dst_vmods['meso'].Segmentation.CNMF.insert(meso_Segmentation_CNMF, skip_duplicates=True, allow_direct_insert=True)

meso_Segmentation_CNMFBackground = (pc.src_vmods['meso'].Segmentation.CNMFBackground & (query_experiment & query_meso)).fetch(as_dict=True)
pc.dst_vmods['meso'].Segmentation.CNMFBackground.insert(meso_Segmentation_CNMFBackground, skip_duplicates=True, allow_direct_insert=True)

meso_Segmentation_Mask = (pc.src_vmods['meso'].Segmentation.Mask & (query_experiment & query_meso)).fetch(as_dict=True)
pc.dst_vmods['meso'].Segmentation.Mask.insert(meso_Segmentation_Mask, skip_duplicates=True, allow_direct_insert=True)

# # TODO remove segmentations that are trials
# # Glomeruli CNMF segmentation
# for key in query_meso:
# key = {'animal_id':1571, 'session':1, 'scan_idx':1}
# key = {'animal_id':120, 'session':1, 'scan_idx':2}
# pc.dst_vmods['meso'].SegmentationTask.insert(
#                               {'animal_id': key['animal_id'],
#                               'session': key['session'],
#                               'scan_idx': key['scan_idx'],
#                               'field': field,
#                               'channel': 1,
#                               'segmentation_method': 2,
#                               'compartment': 'glomerulus'}
#           for field in range(1, (pc.dst_vmods['meso'].ScanInfo & key).fetch1('nfields')+1))

# pc.dst_vmods['meso'].Segmentation.populate(**populate_settings)

meso.Fluorescence.populate((query_experiment & query_meso), **populate_settings)

meso.MaskClassification.populate((query_experiment & query_meso), **populate_settings)

meso.ScanSet.populate((query_experiment & query_meso), **populate_settings)

meso.Activity.populate((query_experiment & query_meso), **populate_settings)

meso.ScanDone.populate((query_experiment & query_meso), **populate_settings)

# ------------------------------------------------------------------------------
print('Odor pipeline')

odor_Odorant = pc.src_vmods['odor'].Odorant.fetch(as_dict=True)
pc.dst_vmods['odor'].Odorant.insert(odor_Odorant, skip_duplicates=True)

odor_OdorSolution = pc.src_vmods['odor'].OdorSolution.fetch(as_dict=True)
pc.dst_vmods['odor'].OdorSolution.insert(odor_OdorSolution, skip_duplicates=True)

odor_OdorSession = (pc.src_vmods['odor'].OdorSession & (query_meso & query_experiment)).fetch(as_dict=True)
pc.dst_vmods['odor'].OdorSession.insert(odor_OdorSession, skip_duplicates=True)

odor_OdorConfig = (pc.src_vmods['odor'].OdorConfig & (query_meso & query_experiment)).fetch(as_dict=True)
pc.dst_vmods['odor'].OdorConfig.insert(odor_OdorConfig, skip_duplicates=True)

odor_OdorRecording = (pc.src_vmods['odor'].OdorRecording & (query_meso & query_experiment)).fetch(as_dict=True)
pc.dst_vmods['odor'].OdorRecording.insert(odor_OdorRecording, skip_duplicates=True)

odor_MesoMatch = (pc.src_vmods['odor'].MesoMatch & (query_meso & query_experiment)).fetch(as_dict=True)
pc.dst_vmods['odor'].MesoMatch.insert(odor_MesoMatch, skip_duplicates=True)

odor.OdorTrials.populate(query_meso & query_experiment, **populate_settings)

odor.OdorSync.populate(query_meso & query_experiment, **populate_settings)

odor.Respiration.populate(query_meso & query_experiment, **populate_settings)

# TODO fix keysource.  make sure dataman is set up for correct experimental id.
odor.OdorAnalysis.populate(**populate_settings)

# TODO run when meso.stitch is complete
odor.SummaryImageSet.populate(**populate_settings)

# ------------------------------------------------------------------------------
print('Treadmill pipeline')

treadmill.Treadmill.populate((query_experiment & query_meso), **populate_settings)

# ------------------------------------------------------------------------------
print('Stack pipeline')

# Create stack restriction for scan-to-stack registration
query_animal = {'animal_id':124}
query_stack = {**query_animal, 'session':9, 'scan_idx':3}
query_meso = (pc.dst_vmods['experiment'].Scan & query_animal & 'aim="2pScan"' & pc.dst_vmods['odor'].MesoMatch).proj()

experiment_Scan = (pc.src_vmods['experiment'].Scan & query_stack).fetch1()
experiment_Scan_Laser = (pc.src_vmods['experiment'].Scan.Laser & query_stack).fetch1()

pc.dst_vmods['experiment'].Stack.insert1({'animal_id': query_stack['animal_id'], 
                                          'session': query_stack['session'], 
                                          'stack_idx': query_stack['scan_idx'],
                                          'lens': experiment_Scan['lens'],
                                          'brain_area': experiment_Scan['brain_area'],
                                          'aim': experiment_Scan['aim'],
                                          'software': experiment_Scan['software'],
                                          'version': experiment_Scan['version'],
                                          'surf_depth': experiment_Scan['depth'],
                                          'top_depth': 0, # TODO find out this data, not required downstream
                                          'bottom_depth': 0, #TODO find out this data, not required downstream
                                          'stack_notes': experiment_Scan['scan_notes'],
                                          'stack_ts': experiment_Scan['scan_ts']})

pc.dst_vmods['experiment'].Stack.Filename.insert1({'animal_id': query_stack['animal_id'], 
                                                   'session': query_stack['session'], 
                                                   'stack_idx': query_stack['scan_idx'],
                                                   'filename_idx': 1,
                                                   'filename': experiment_Scan['filename']})

pc.dst_vmods['experiment'].Stack.Laser.insert1({'animal_id': query_stack['animal_id'], 
                                                'session': query_stack['session'], 
                                                'stack_idx': query_stack['scan_idx'],
                                                'wavelength': experiment_Scan_Laser['wavelength'],
                                                'max_power': experiment_Scan_Laser['power'],
                                                'gdd': experiment_Scan_Laser['gdd']})

stack.StackInfo.populate(**populate_settings)

stack.Quality.populate(**populate_settings)

# Required if number of channels > 1
# pc.dst_vmods['stack'].CorrectionChannel.insert1()

stack.RasterCorrection.populate(**populate_settings)

stack.MotionCorrection.populate(**populate_settings)

stack.Stitching.populate(**populate_settings)

stack.CorrectedStack.populate(**populate_settings)

stack.PreprocessedStack.populate(**populate_settings)

# TODO error
# stack.Surface.populate(**populate_settings)

stack_CorrectedStack = (pc.dst_vmods['stack'].CorrectedStack & \
                        {'animal_id': query_stack['animal_id'], \
                        'session': query_stack['session'], \
                        'stack_idx': query_stack['scan_idx']}).fetch()

meso_ScanInfo_Field = (pc.dst_vmods['meso'].ScanInfo.Field & query_meso).fetch()

for corrected_stack in stack_CorrectedStack:
     for field in meso_ScanInfo_Field:
          pc.dst_vmods['stack'].RegistrationTask.insert1({'animal_id': corrected_stack['animal_id'],
                                                          'stack_session': corrected_stack['session'],
                                                          'stack_idx': corrected_stack['stack_idx'],
                                                          'volume_id': corrected_stack['volume_id'], 
                                                          'stack_channel': 1,
                                                          'scan_session': field['session'],
                                                          'scan_idx': field['scan_idx'],
                                                          'scan_channel': 1,
                                                          'field': field['field'],
                                                          'registration_method': 5})# TODO Default is 5, clarify which method to use

# TODO running
stack.Registration.populate(**populate_settings)

# TODO add subsequent populate steps?