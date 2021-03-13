import datajoint as dj
from . import mice

schema = dj.schema(dj.config['database.prefix'] + 'pipeline_map')

@schema
class RetMap (dj.Manual):
     definition = """
          # Retinotopy map
          -> mice.Mice
          ret_idx : smallint        # retinotopy map index for each animal
          ---
     """

     # Convert MATLAB function here