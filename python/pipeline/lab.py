import os
import sys
import numpy as np

import datajoint as dj

dj.config['database.prefix'] = os.environ.get('DJ_PREFIX', '')
schema = dj.schema(dj.config['database.prefix'] + 'pipeline_lab')

@schema
class Paths(dj.Lookup):
    definition = """
        global       : varchar(255)               # global path name
        ---
        linux        : varchar(255)               # linux path name
        windows      : varchar(255)               # windows path name
        mac          : varchar(255)               # mac path name
        location     : varchar(255)               # computer path
    """

    def get_local_path(self, path, local_os=None):

        if path.split(':')[0] == 'L':
            root_dir = os.environ.get('SCRATCH08_STORAGE')
        elif path.split(':')[0] == 'N':
            root_dir = os.environ.get('SCRATCH10_STORAGE')

        path = os.path.join(root_dir, path.split(':')[1])
        
        path = path.replace('\\', '/')

        return path