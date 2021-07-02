import os
import datajoint as dj
from . import map, experiment, shared

dj.config['database.prefix'] = os.environ.get('DJ_PREFIX', '')
schema = dj.schema(dj.config['database.prefix'] + 'pipeline_anatomy')

@schema
class Area (dj.Lookup):
     definition = """
          brain_area : varchar(12)   # short brain area name
     """
     contents = [
            ['V1'],
            ['P'],
            ['POR'],
            ['PM'],
            ['AM'],
            ['A'],
            ['RL'],
            ['AL'],
            ['LI'],
            ['LM'],
            ['OB'],
            ['PC']
     ]

@schema
class AreaMask (dj.Manual):
     definition = """
          # Area mask for each scan
          -> experiment.Scan
          -> Area
          -> shared.Field
          ---
          -> map.RetMap
          mask                     : mediumblob            # mask of area
     """

     # Convert MATLAB function here