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

# # Summary notebook
# This notebook is used to determine if all data for a given session is acquired and processed.

import os
import datajoint as dj
dj.config['database.prefix'] = os.environ.get('DJ_PREFIX', '')
from pipeline import experiment, meso, odor, treadmill

# # Enter experimental identifier

experiment_id = 3

# # Summary of all entries for experimental identifier

experiment_id_query = experiment.ExperimentalIdentifier() & f'experiment_id={experiment_id}'
experiment_id_query

# +
if len(experiment_id_query) == 1 and \
   len(meso.QualityManualCuration & experiment_id_query) == 1 and \
   len(meso.Segmentation & experiment_id_query) > 0  and \
   len(meso.Fluorescence & experiment_id_query) > 0 and \
   len(meso.Stitch & experiment_id_query) == 1 and \
   len(odor.Respiration & experiment_id_query) == 1 and \
   len(odor.OdorAnalysis & experiment_id_query) == 1 and \
   len(treadmill.Treadmill & experiment_id_query) == 1:
    print(f'All entries exist for experimental id: {experiment_id}')
else:
    print(f'All entries do NOT exist for experiment id: {experiment_id}.\
          \nPlease investigate missing entries using tables below.')

# TODO: add pupil diameter
# -

# # Entries for experimental identifier in individual tables

meso.QualityManualCuration & experiment_id_query

meso.Segmentation & experiment_id_query

meso.Fluorescence & experiment_id_query

meso.Stitch & experiment_id_query

odor.Respiration & experiment_id_query

odor.OdorAnalysis & experiment_id_query

treadmill.Treadmill & experiment_id_query

# +
# TODO: add pupil diameter
# -
# # Summary of all entries
# Make sure the entries are populated for all downstream tables.

experiment.Session.proj()


experiment.Scan.proj()

experiment.Session.proj() - experiment.Scan.proj() #TODO delete entries from Session





experiment.Scan.proj() - meso.ScanInfo.proj()

experiment.Scan.proj() - meso.Quality()

meso.CorrectionChannel.proj() - meso.RasterCorrection.proj()

meso.CorrectionChannel.proj() - meso.MotionCorrection.proj()

experiment.Scan.proj() - meso.Stitch.proj()

meso.SummaryImages()

meso.SegmentationTask.proj() - meso.Segmentation.proj()

meso.SegmentationTask.proj() - meso.Fluorescence.proj()

meso.MaskClassification()

meso.SegmentationTask.proj() - meso.ScanSet.proj()

meso.Activity()

meso.ScanDone()

odor.OdorRecording()

odor.OdorSession()

odor.MesoMatch()




