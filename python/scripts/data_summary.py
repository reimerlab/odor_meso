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
# This notebook is used to determine if all data is populated for the `at-database` to `jr-database` transfer.

import os
import datajoint as dj
dj.config['database.prefix'] = os.environ.get('DJ_PREFIX', '')
from pipeline import mice, experiment, meso, odor, treadmill

# ---

mice.Mice.proj()

experiment.Session.proj()

experiment.Session.proj() - odor.MesoMatch()

experiment.Scan.proj()

experiment.Scan.proj() - odor.MesoMatch()

experiment.AutoProcessing()

experiment.AutoProcessing() - odor.MesoMatch()

experiment.AutoProcessing.proj() - meso.ScanInfo.proj()

experiment.AutoProcessing.proj() - meso.Quality()

meso.ScanInfo.Field.proj() - meso.CorrectionChannel.proj()

meso.ScanInfo.Field.proj() - meso.RasterCorrection.proj()

meso.ScanInfo.Field.proj() - meso.MotionCorrection.proj()

experiment.AutoProcessing.proj() - meso.Stitch.proj()

meso.CorrectionChannel() - meso.SummaryImages()

meso.SegmentationTask.proj() - meso.Segmentation.proj()

meso.SegmentationTask.proj() - meso.Fluorescence.proj()

meso.MaskClassification()

meso.SegmentationTask.proj() - meso.ScanSet.proj()

meso.Activity()

meso.ScanDone()

experiment.AutoProcessing.proj() - treadmill.Treadmill.proj()

odor.OdorRecording.proj() - odor.MesoMatch()

odor.MesoMatch.proj() - odor.OdorTrials()

odor.MesoMatch.proj() - odor.OdorSync.proj()

odor.MesoMatch.proj() - odor.Respiration.proj()

experiment.ExperimentalIdentifier() - odor.OdorAnalysis()

experiment.ExperimentalIdentifier() - odor.SummaryImageSet()

# +
# TODO: add stack pipeline
# -


