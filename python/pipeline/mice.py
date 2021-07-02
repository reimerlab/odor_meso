"""
The `common_mice` schema is maintained by another package and is included here for ease of reference.
DO NOT create new tables here.
"""
import os
import datajoint as dj

dj.config['database.prefix'] = os.environ.get('DJ_PREFIX', '')
schema = dj.schema(dj.config['database.prefix'] + 'pipeline_mice')

@schema
class Mice(dj.Manual):
    definition = """  # calcium-sensitive indicators
          animal_id           : int                                            # id number
          ---
          other_id=""         : varchar(20)                                    # alternative id number
          dob=null            : date                                           # animal's date of birth
          dow=null            : date                                           # animal's date of weaning
          sex="unknown"       : enum('M','F','unknown')                        # animal's sex
          color="unknown"     : enum('Black','Brown','White','unknown')        # animal's color
          ear_punch="unknown" : enum('None','R','L','RL','RR','LL','unknown')  # animal's ear punch
          owner="none"        : enum('Federico','Fabian','Jiakun','Jake','Manolis','Xiaolong','Dimitri','Shan','Keith','Cathryn','Deumani','Matt','Megan','Paul','Shuang','Other','Available','none') # mouse's owner
          facility="unknown"  : enum('TMF','Taub','Other','unknown')           # animal's curent facility 
          room="unknown"      : enum('VH1','T027','VK3','VD4','T014','T057','T086D','Other','unknown') # animal's current room 
          rack=null           : tinyint                                        # animal's curent rack 
          row=null            : tinyint                                        # animal's curent row
          mouse_notes=""      : varchar(4096)                                  # other comments and distinguishing features
          mouse_ts=CURRENT_TIMESTAMP : timestamp                               # automatic
    """