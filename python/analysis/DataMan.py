import csv
import os, re, sys
import time
import pickle
# from scipy.signal import butter, filtfilt   # savgol_filter
from sklearn import linear_model
from matplotlib import pyplot as plt
# from matplotlib.backend_bases import MouseEvent
from matplotlib import collections  as mc
import numpy as np
from scipy import stats
import helpers   # A local module
from helpers import BADFILENAME
import logging
from filters import smooth, slider, band_pass, data_reduce, high_pass    # A local module
import pkg_resources
import copy
from pprint import pprint
import stats      # A local module
import quantities as pq
try:
    pkg_resources.get_distribution('datajoint')
except pkg_resources.DistributionNotFound:
    datajoint_ok = False
else:
    import datajoint as dj
    datajoint_ok = True
from matplotlib.widgets import CheckButtons
#from scipy.ndimage import gaussian_filter
#from scipy.ndimage.morphology import binary_dilation
#import matplotlib.animation as animation
#from tqdm import tqdm
#import seaborn as sns
from scipy import interpolate, signal
from scipy.stats import norm, zscore
#from matplotlib import cm
import pandas as pd
import math

print("DataJoint Module Installed = ", datajoint_ok)
class BADDATAJOINT (Exception): pass


# Suppresses some annoying matplotlib messages
logging.getLogger().setLevel(logging.CRITICAL)

class DataHandler():
    '''
    Defines the general type and provide a clean way to query what this is using what_am_i property
    '''
    def __init__(self, *args, **kwargs):
        if "silence" in kwargs.keys():
            self.silence = kwargs.pop("silence")
        else:
            self.silence = False
        self.connection_status = {}
        self.dj_objects = []
        self.start_time_shift = pq.Quantity(0.0, pq.sec)  # Default for time alignment between stimuli and fluorescence

        self.expts = None
        self._import_from = None  #"datajoint"   # Source of the data either dj, csv or pickle
        self.db_upload_available = False
        self.my_paths = None
        self.pop_arg = None
        for arg in args:
            if not self.silence: print("DataHandler got", type(arg), arg)
        self._what_i_am = "DataHandler"
        # Which record in expt_id
        self._expt_id = None

        try:
            self.my_paths = helpers.My_paths(*args, **kwargs)
        except BADFILENAME:
            print("Bad path.init DataMan construction aborted")
            path_to_init = input("Provide location of paths.init or enter to abort construction")
            if path_to_init != "":
                self.my_paths = helpers.My_paths(location, **{"paths_init": path_to_init})
            else:
                raise BADFILENAME

        for i, arg in enumerate(args):
            if isinstance(arg, (str)):
                if arg in ("datajoint", "dj", "db", "data_joint", "database", "d"):
                    self._import_from = "dj"
                    if not self.silence: print("Building with DataJoint", self._import_from)
                elif arg in ("pickle", "dict", "p"):
                    self._import_from = "pickle"
                    if not self.silence: print("Building with Pickle")
                elif arg in ("csv", "file", "excel", "c"):
                    self._import_from = "csv"
                    if not self.silence: print("Building with CSV")
                elif arg in ("DataMan", "dm"):
                    self._what_i_am = "DataMan"

            elif isinstance(arg, (int)):
                # if we have the info for which expt then load and prepare to pop it from args list when used
                self._expt_id = arg
                self.pop_arg = i

        for key in kwargs.keys():
            if key in ("make_using", "import_from"):
                self._import_from = kwargs[key]
            if key == "expt_id":
                self._expt_id = kwargs[key]
            if key == "what_i_am":
                self._what_i_am = kwargs[key]

        self.db_connect(**kwargs)


    @property
    def what_am_i(self):
        return self._what_i_am

    @property
    def expt_id(self):
        return self._expt_id

    @expt_id.setter
    def expt_id(self, new_id):
        self._expt_id = new_id
        if self.expts is not None and new_id in self.expts.index:
            self.restriction = self.select_expt(new_id)

    @property
    def how_made(self):
        return self._import_from

    def find_expts(self):
        '''
        Pulls out a list of available experiments by animal, session, recording id's
        :return: loads info into DataMan object
        '''
        self.expts = pd.DataFrame((self.odor.MesoMatch).fetch())


    def select_expt_id(self, expt_id):
        if expt_id in self.expts.index:
            self.restriction = dict(self.expts.loc[expt_id])
            return dict(self.expts.loc[expt_id])
        else:
            return {}

    def select_expt(self, row):
        if row < self.expts.shape[0] and row > 0:
            self.restriction = dict(self.expts.iloc[row,:])
            return dict(self.expts.iloc[row,:])
        else:
            return {}

    def db_connect(self, **kwargs):
        '''
        Sets up connection to DataJoint and returns status
        :param kwargs: Provide alternative config info {"host": "strhost", "user": struser, "pwd": strpwd}
        :return: status of connection put into self.db_upload_available
        '''
        host = self.my_paths.datajoint_db
        user_name = self.my_paths.username
        password = self.my_paths.password

        if kwargs:
            for key, value in kwargs.items():
                if key == "host":
                    host = value
                elif key == "user":
                    user_name = value
                elif key == "pwd":
                    password = value
        try:
            dj.config['database.host']= host
            dj.config['database.user'] = user_name
            dj.config['database.password'] = password
            dj.conn()
        except:
            print("Problem connecting to DataJoint ", sys.exc_info()[0])
            self.db_upload_available = False
        else:
            self.db_upload_available = True

        self.db_upload_available = self.test_db_sockets()

    def test_db_sockets(self):
        '''
        This method tests the connection to the datajoint database
        And if good produces a DataFrame of the available experiments.
        :return: Connection status for db
        '''
        if self.db_upload_available:
            # Attempt to make a connection with the essential odor pipeline of the DataJoint db
            try:
                self.odor = dj.create_virtual_module('', 'pipeline_odor')
            except:
                print("DB_test: Problem with DataJoint Connection.  Recheck to see if VPN or network has a problem.")
                return False

            if not self.silence: print("Datajoint odor connection passes test")
            self.connection_status["odor"] = True
            self.dj_objects.append("odor")

            self.find_expts()
            return True

        else:
            print("DataJoint connection information not loaded. First run self.db_connect()")
            return False


class DataMan(DataHandler):
    '''
    This is the base class for managing odor experiment Data.  Creates a set of tables from a data source
    ultimately linked back to DataJoint Database.  It uses pq.Quantities where ever possible to allow
    unit sanity checks, however some numpy and scipy methods gag on pq.Quantities so they are temporarily
    stripped off.
    '''
    rfu = pq.UnitQuantity('relative_fluorescence', symbol='rfu')
    msec = pq.millisecond
    df_f_unit = pq.dimensionless
    testing = False

    def __init__(self, *args, **kwargs):   #
        '''
        Most parameters in the init must be set on object initialization.  Generally the creation
        assumes that the program has access to the paths.init file to understand where things are
        and uinit.init which provide general initialization parameters
        If no args or kwargs are given then the object will be made but later needs to be loaded
        with data.
        :param *args: if a number will load a specific experiment from DataJoint.  Or a key to which
        source to use to build object.
        :param **kwargs: override values for parameters in init_file.  Also a method to reload a prevous analysis by passing reload=True
        possible values: fstim="set1.csv", fdata="1571df_f.csv", fodors="odorset1.csv", sort_by="duration", use_odor="all", use_glomeruli="all", dt=0.077886473, corner_f=1.0, buffer=20, start_plot=0):

        '''
        # Put in dummies for the datajoint pipeline objects
        self.meso = None
        self.odor = None
        self.stimulus = None
        self.treadmill = None
        self.movie = None
        self.experiment = None
        self.schema = None

        try:
            super().__init__(*args, **kwargs)
        except BADFILENAME:
            print("paths.init information is wrong please fix. Build aborted")
            return
        self._what_i_am = "DataMan"
        # Removes the expt_id from args if there.
        if self.pop_arg is not None and isinstance(self.pop_arg, (int)):
            argt = list(args)
            argt.pop(self.pop_arg)
            args = tuple(argt)
        for arg in args:
            if not self.silence: print(arg)
        # Default list of things available for autoimport if desired
        self.autoimport = ["fluorescence", "odor_trials", "treadmill", "respiration"]
        self.understood_units = ["sec", "msec", "Hz", "rfu", "df_f"]
        self.known_masks = ["unmasked", "pre_stim", "pre_maxodor_post", "pre_odor_post", "odor_post", "odor", "common", "maxodor_post"]
        self.known_f0_methods = ["bkg_smooth", "low_pass", "constant_avg", "constant_med", "linear", "combo"]
        self.available_restrictions = ['animal_id', 'odor_session', 'recording_idx','session', 'scan_idx', 'pipe_version']
        self.restriction = {}
        self.df_fs = {}   # Dict Container for df_f types
        self.f0s = {}    # Dict Container for f0 calculations
        self.combo = ["bkg_smooth", "constant_med"]
        self.exp_dict = {}
        self.print_multiple_stimuli = False
        self.exp_dict["mask_dict"] = {}   # Store for created masks
        self.exp_dict["data_subset"] = {} # Store for created subsets of the experimental data
        self.response_dict = {}
        self.sham_resp_dict = {}
        self.available_data = []  # List containing types of data available to use in this experiment
        self.verify_input_type = False  # Interactive flag to ask whether the data is of the expected type dF/f or raw f
        self.label_conc = True   #Add concentration to label name
        #self.reload = False
        self.sig_prob = 0.99   # Significance Criterion for saying response integral is a potential positive
        #self.data_masked = False
        self.global_baseline = None
        #self.buffer = None
        self.response_mask = "common"  # if "common" will treat all stimuli regardless of duration as the same.  Can change in uinit.init
        self.c_start = 0.5 * pq.s  # common start time after beginning of odor stimulus
        self.c_end = 2.0 * pq.s  # common end time after beginning of odor stimulus
        #Parameter/Variables initialized in init file
        import_from = None
        expt_id = None
        if "uinit" in kwargs.keys():
            init_file = kwargs["unit"]
        else:
            init_file = "uinit.init"
        self.init_dict = helpers.process_init(init_file, int_keys_ok=True, my_paths=self.my_paths)
        self.exp_dict["init_dict"] = self.init_dict
        for key, value in self.init_dict.items():
            if key == "import_from":
                import_from = value
                continue
            elif key == "expt_id":
                expt_id = value
                continue
            if isinstance(value, (list, tuple)):
                if value[-1] in self.understood_units:
                    value = make_quantity(value)
            #print("attribute", key, value)
            setattr(self, key, value)
        # can pass kwargs to override init file or in place of reading an init file
        if not self.silence: print("kwargs", kwargs)
        for key in kwargs:
            value = kwargs[key]
            if key == "import_from":
                import_from = value
                continue
            elif key == "expt_id":
                expt_id = value
                continue
            if isinstance(value, (list, tuple)):
                if value[-1] in self.understood_units:
                    value = make_quantity(value)
            setattr(self, key, value)


        if self._expt_id is None and import_from is not None and expt_id is not None:
            self._import_from = import_from
            self._expt_id = expt_id
        if self._import_from in [None, ""] and self._expt_id in [None, ""]:
            print("DataPlot ready to build.")
            print("no data imported, DataPlot object is empty.  Need to fill before use.")
            return
        if not self.silence: print("Import Instructions", self._import_from, self._expt_id)
        # Import Block.  Choose data import method, or wait for separate import instructions
        if isinstance(self._expt_id, (int)) and (self._import_from is None or self._import_from == ""):
            if not self.silence: print("DataMan ready to build from default DataJoint.")
            self.build_from_dj(*args, **kwargs)
        elif self._import_from.lower() in ("pickle", "dict", "p"):
            if not self.silence: print("DataMan will be built from pickle file")
            self.build_from_pickle()
        elif self._import_from.lower() in ("datajoint", "dj", "db", "data_joint", "database") and datajoint_ok:
            self.build_from_dj(*args, **kwargs)

        elif self._import_from.lower() in ("csv", "file", "excel", "c"):
            self.build_from_csv(*args, **kwargs)

        else:
            print("Don't understand or non-existant import instructions ", self._import_from, self._expt_id)
            print("no data import, DataMan object is empty.  Need to fill before use.")

    @property
    def df_f_raw(self):
        return self.df_fs[self.f0_method]

    def loaded_data(self):
        print("Experiment {0}".format(self.expt_id))
        print("Retrieved using Restriction:")
        pprint(self.restriction)
        print("Retrieved data for: ", self.available_data)

    def top(self, p=None):
        if p is not None and p != self.sig_prob:
            self.exp_dict["pd_candidates"] = self.pd_resp(p)
            self.exp_dict["response_cand"] = self.id_candidates()
        print("Candidate list with p= ", p)
        summary = self.exp_dict["response_cand"]["summary"]
        positive = self.exp_dict["response_cand"]["positive"]
        negative = self.exp_dict["response_cand"]["negative"]
        for odor in summary.keys():
            print("Odor: {0} has {1} positive and {2} negative candidates".format(odor, summary[odor]["n_pos"], summary[odor]["n_neg"]))
        print('The positive candidates are:')
        for odor in summary.keys():
            if len(positive[odor]) > 0:
                print(odor, positive[odor])
        print("")
        print('The negative candidates are:')
        for odor in summary.keys():
            if len(negative[odor]) > 0:
                print(odor, negative[odor])    


    def build_from_dj(self, *args, **kwargs):
        '''
        Runs through the steps in preparing DataMan object for use with DataJoint data import
        Steps are:
            self.import_datajoint_data(**kwargs)  #DataJoint data imported
            self.datajoint_associate_odors_nozzles()  # DataJoint imported odors associated with nozzles
            self.process_stimuli_slice_and_mask()  #Concurrent stimuli merged, datamasks and slices made, and active_set of glomeruli identified
            self.make_df_f(test=self.test, plot_it=self.f0_auto_plot)   # df/f data made
            self.smooth_and_make_baseline()  # smoothing applied, baseline made, ready to plot
            self.initialize_exp_dict()  #exp_dict made, concatenate responses made and integrals completed

        :param kwargs: controls for DJ import. Can make with helpers.new_exp:
                        {"animal_id": animal, "odor_session": session, "recording_idx": recording}
        :return: DataJoint object ready to go
        '''
        db_status_ok = self.connect_dm_to_dj()
        if not db_status_ok:
            print("Halting build, no DatatJoint DB connection")
            return
        import_status_ok = self.import_datajoint_data(*args, **kwargs)
        if import_status_ok:
            if not self.silence: print("DataJoint data imported")
        else:
            print("Import info problems.  Build of DataMan stopped")
            return
        if "odor_trials" in self.available_data and self.pd_odor_trials.shape[0] != 0:
            self.process_odors()
            if not self.silence: print("Odors processed")
        else:
            print("No odor trials, stopping build")
            return
        if self.use_glomeruli == "all":
            self.glomeruli = [i for i in range(1, self.raw_fluorescence.shape[1]+1)]
        else:
            self.glomeruli = self.use_glomeruli
        #self.slice_maker()
        if not self.silence: print("extering mask_maker")
        self.mask_maker()
        if not self.silence: print("Concurrent stimuli merged, and odors linked to nozzles")
        #self.f0s = self.make_f0(self.f0_method)
        #self.df_fs = self.make_df_f(self.f0s)
        if not self.silence: print("Calculating f0 and df/f")
        if self.f0_method is not None and self.f0_method not in ("compare", "all", "test"):
            self.set_f0()
            if not self.silence: print("df/f data made")
        else:
            print("set up dm with preferred method using self.set_f0(method)")
        self._import_from = "datajoint"
        self.complete_dm_build()

    def build_from_pickle(self):
        '''
        rebuilds object from pickle of exp_dict.  Resets any common parameter in DataMan object to value in exp_dict
        :return: DataJoint object ready to go
        '''
        self.reload_from_pickle()
        # Align object attributes to the exp_dict imported
        if "state_dict" in self.exp_dict.keys():
            self.old_state_dict = self.exp_dict.pop("state_dict")
        else:
            self.old_state_dict = {}
        # Load all DataMan object params and flags from the saved state of the pickled exp_dict
        for key, value in self.old_state_dict.items():
            if value is "None":
                setattr(self, key, None)
            else:
                setattr(self, key, value)
        # Make sure object arrays have correct values
        for key, value in self.exp_dict.items():
            #print("key", key)
            try:
                my_attribute = getattr(self, key)
                if my_attribute is not value:
                    setattr(self, key, value)
            except AttributeError:
                print("Object attribute ", key, " is only in self.exp_dict")
            except:
                print("unexpected error", sys.exc_info()[0])
        self._import_from = "pickle"


    def build_from_csv(self, *args, **kwargs):
        '''
        Builds DataMan object using csv files for response data with times, odors used with nozzel info and stimuli
        with (odor_id, start_t, stop_t)
        :return: DataMan loaded and ready to plot
        '''
        print("Building from CSV", len(args), len(kwargs))
        if self._expt_id is None:
            try:
                self._expt_id = self.csv_expt_id
            except:
                self._expt_id == 1001
        if self.expts is not None and self._expt_id in self.expts.index:
            self.restriction = self.select_expt(self._expt_id)
        else:
            self.restriction = {}
        self.import_csv_response_data()
        print("csv response data imported")
        self.import_csv_odors()
        print("csv odors imported")
        self.raw_pd_odor = self.import_csv_stimuli()
        print("csv stimuli imported")
        self.pd_odor = self.map_odortrials_to_frames_and_combine(self.raw_pd_odor)
        #self.slice_maker()
        self.mask_maker()



        print("Concurrent stimuli merged, datamasks made, and active_set of glomeruli identified")
        if self.verify_input_type:
            check_import_data = input("Is data type df_f")
            if check_import_data.lower() in ('y', 'yes'):
                self.import_data_type = "df_f"
                if self.raw_fluorescence.units != DataMan.df_f_unit.units:
                    self.raw_fluorescence = pq.Quantity(np.array(self.raw_fluorescence), DataMan.df_f_unit)
                self.f0_method = "pre_processed"
                self.df_fs["pre_processed"] = copy.deepcopy(self.raw_fluorescence)
                self.f0 = make_quantity(np.ones_like(self.df_fs["pre_processed"]), DataMan.rfu)

            elif check_import_data.lower() in ('n', 'no'):
                self.import_data_type = "rfu"
                if self.raw_fluorescence.units != DataMan.rfu.units:
                    self.raw_fluorescence = pq.Quantity(np.array(self.raw_fluorescence), DataMan.rfu)
                self.df_fs = self.make_df_f()
                print("df/f data made")
            else:
                print("Was expecting yes or no but got {}. Will assume this is df_f and build DataMan accordingly ").format(check_import_data)
                self.raw_fluorescence = pq.Quantity(np.array(self.raw_fluorescence), DataMan.df_f_unit)
                self.f0_method = "pre_processed"
                self.df_fs["pre_processed"] = copy.deepcopy(self.raw_fluorescence)
                self.f0 = make_quantity(np.ones_like(self.df_fs["pre_processed"]), DataMan.rfu)
                self.import_data_type = "df_f"
        else:
            print("Assuming csv input is df/f")
            self.raw_fluorescence = pq.Quantity(np.array(self.raw_fluorescence), DataMan.df_f_unit)
            self.f0_method = "pre_processed"
            self.df_fs["pre_processed"] = copy.deepcopy(self.raw_fluorescence)
            self.f0 = make_quantity(np.ones_like(self.df_fs["pre_processed"]), DataMan.rfu)
            self.import_data_type = "df_f"
        self._import_from = "csv"
        if self.break_early: return
        self.complete_dm_build()

    def set_f0(self, good_f0=None, force=False):
        '''
        This method sets what the f0 trace will be.  If force is set, it will clear the existing f0 trace of this type and
        recalculate f0.  Generally follows the instructions in the self.f0_method parameter.  Likewise, if the f0 of the
        instructed method does not already exist this will calculate it first.
        :param good_f0:
        :param force:
        :return:
        '''
        if good_f0 is None:
            good_f0 = self.f0_method
        if force:
            if good_f0 in self.f0s.keys():
                del self.f0s[good_f0]
                del self.df_fs[good_f0]
            elif good_f0 == "combo" and self.combo[0] in self.f0s.keys:
                del self.f0s[self.combo[0]]
                if "combo" in self.df_fs.keys():
                    del self.df_fs["combo"]
        if good_f0 == "combo" and good_f0 not in self.df_fs.keys():
            if self.combo[0] in self.f0s.keys():
                self.f0 = self.f0s[self.combo[0]]
            else:
                self.f0s = self.make_f0(good_f0)
            self.df_fs = self.make_df_f(good_f0)
            self.f0_method = good_f0

        elif good_f0 == "combo":
            self.f0 = self.f0s[self.combo[0]]
            self.f0_method = good_f0

        elif good_f0 not in self.f0s.keys():
            self.f0_method = good_f0
            self.f0s = self.make_f0(good_f0)
            self.f0 = self.f0s[good_f0]
            self.df_fs = self.make_df_f(good_f0)

        elif good_f0 in self.df_fs.keys():
            if good_f0 != "combo":
                self.f0 = self.f0s[good_f0]
            else:
                self.f0 = self.f0s[self.combo[0]]
            self.f0_method = good_f0
        else:
            print("good_f0 not in available f0s", self.f0s.keys())
            return False


    def complete_dm_build(self):
        s_time = time.time()
        self.smooth_and_make_baseline()

        if not self.silence: print("smoothing applied, baseline made, ready to plot")
        self.initialize_exp_dict()
        if not self.silence: print("exp_dict made, concatenated responses made and integrals completed")

    def connect_dm_to_dj(self, *args, **kwargs):
        '''
        Connects DataMan to DataJoint using animal info imported or provided by kwargs
        :param kwargs: controls for DJ import. Can make with helpers.new_exp:
                        {"animal_id": animal, "odor_session": session, "recording_idx": recording}
        :return:
        '''
        if not self.db_upload_available:
            self.db_connect()
            if not self.db_upload_available:
                print("Problem with DataJoint connection, will have to reconnect to build DataMan object")
                return False
        self.connect_db_sockets()
        if not self.silence: print("Datajoint pipelines connected")
        return True

    def smooth_and_make_baseline(self):
        '''

        :param self:
        :return:
        '''
        if not self.silence: print("Smoothing data and making baseline")
        self.sdata = smooth(self.df_fs[self.f0_method], self.corner_f, self.avg_dt, axis=0)
        if self.global_baseline_method == "low_pass":
            self.global_baseline = smooth(self.df_fs[self.f0_method], self.baseline_corner_f, self.avg_dt, axis=0)

        else:
            masking = self.mask_slice_manager("odor_post")
            kwargs = {"times":self.times, "sdata": self.sdata, "odor_mask": copy.deepcopy(masking["mask"]), "interpolate": True}
            self.global_baseline = bkg_smooth(kwargs, bkg_corner_f=self.bkg_smooth_corner_f, kind=self.bkg_smooth_interp, silence=self.silence)

    def initialize_exp_dict(self):
        '''
        Used to build the initial experimental dict after all essential data is loaded
        :return:
        '''
        ied_start = time.time()
        if not self.silence: print("initializing exp_dict")
        exp_dict = self.initial_fill_experiment_dict()
        t1 = time.time()
        self.exp_dict = {**self.exp_dict, **exp_dict}
        #Concatenate the responses to each odor
        #masker = self.mask_slice_manager(how_to=odor)
        if self.testing: return
        t2 = time.time()
        self.exp_dict["response_integral"] = self.simple_integral()
        self.autoscale_responses = [self.exp_dict["sdata-baseline"].min(), self.exp_dict["sdata-baseline"].max()]
        self.autoscale_integrals = [self.exp_dict["response_integral"].min(), self.exp_dict["response_integral"].max()]

        t3 = time.time()
        self.exp_dict["pd_candidates"] = self.pd_cand(self.sig_prob)
        t4 = time.time()
        self.exp_dict["response_cand"] = self.id_candidates()
        t5 = time.time()
        self.make_response_data(self.response_mask)
        t6 = time.time()
        # Add columns to pd_odor that id problematic stimuli
        self.pd_odor = self.characterize_stims()
        #Print out to monitor how much time each step is taking
        if self.testing: print("ied times", t1-ied_start, t2-t1, t3-t2, t4-t3, t5-t4, t6-t5, time.time()-t6)


    def connect_db_sockets(self):
        '''
        This method builds the connection to the datajoint database
        odor is connected by test_db_sockets method
        :return: important data connection objects for the database
        '''
        if self.db_upload_available:
            self.connection_status["meso"] = self.connect_meso()
            if "treadmill" in self.autoimport:
                self.connection_status["treadmill"] = self.connect_treadmill()
                self.connection_status["stimulus"] = self.connect_stimulus()
            if "experiment" in self.autoimport:
                self.connection_status["experiment"] = self.connect_experiment()
            if "movie" in self.autoimport:
                self.connection_status["movie"] = self.connect_movie()
            if "schema" in self.autoimport:
                self.connection_status["schema"] = self.connect_schema()

        else:
            print("DataJoint connection information not loaded. First run self.db_connect()")

    def connect_meso(self):
        try:
            self.meso = dj.create_virtual_module('', 'pipeline_meso')
            if not self.silence: print("Meso Connected")
            self.dj_objects.append("meso")
            return True
        except:
            print("Problem with meso Connection.")
            return False

    def connect_odor(self):
        try:
            self.odor = dj.create_virtual_module('', 'pipeline_odor')
            if not self.silence: print("Odor Connected")
            self.dj_objects.append("odor")
            return True
        except:
            print("Problem with odor Connection.")
            return False

    def connect_treadmill(self):
        try:
            self.treadmill = dj.create_virtual_module('', 'pipeline_treadmill')
            if not self.silence: print("Treadmill Connected")
            self.dj_objects.append("treadmill")
            return True
        except:
            print("Problem with treadmill Connection.")
            return False

    def connect_stimulus(self):
        try:
            self.stimulus = dj.create_virtual_module('', 'pipeline_stimulus')
            if not self.silence: print("Stimulus Connected")
            self.dj_objects.append("stimulus")
            return True
        except:
            print("Problem with stimulus Connection.")
            return False

    def connect_movie(self):
        try:
            self.movie = dj.create_virtual_module('', 'pipeline_tune')
            if not self.silence: print("Movie Connected")
            self.dj_objects.append("movie")
            return True
        except:
            print("Problem with movie Connection.")
            return False

    def connect_experiment(self):
        try:
            self.experiment = dj.create_virtual_module('', 'pipeline_experiment')
            if not self.silence: print("Experiment Connected")
            self.dj_objects.append("experiment")
            return True
        except:
            print("Problem with experiment Connection.")
            return False

    def connect_schema(self):
        try:
            self.schema = dj.create_virtual_module('', 'pipeline_schema')
            if not self.silence: print("Schema Connected")
            self.dj_objects.append("schema")
            return True
        except:
            print("Problem with schema Connection.")
            return False

    def get_restriction(self, *args, **kwargs):
        '''
        Routine to establish a good restriction on DataJoint to get one expt
        :param args:  possible expt_id's to a specific restriction
        :param kwargs: possible restriction keys and values
        :return: a good restriction or False
        '''
        if self._expt_id is not None and self._expt_id in self.expts.index:
            restriction = self.select_expt_id(self._expt_id)
        else:
            restriction = {}
        # The Default restrictions can be overridden by args and kwargs
        for arg in args:
            if isinstance(arg, (int)) or (isinstance(arg, (str)) and arg.isnumeric()):
                self._expt_id = int(arg)
                restriction = self.select_expt_id(arg)
                break
        tmp_restriction = {}
        for key in kwargs.keys():
            if key in self.available_restrictions:
                tmp_restriction[key] = kwargs[key]
        if len(tmp_restriction) > 0:
            query = ' and '.join([f'{k} == {repr(v)}' for k, v in tmp_restriction.items()])
            id = self.expts.query(query)
            if id.shape[0] == 1:
                self._expt_id = int(id.index[0])
                restriction = tmp_restriction
            elif  id.shape[0] == 0:
                restriction = {}
            else:
                print("Multiple experiments selected with query restrictions")
                print("Please provide restriction to a single expt")
                print("query", kwargs)
                print(id)
                return False
        if len(restriction) == 0:
            print("Data = Restriction Needed")
            print("Available Expts")
            print(self.expts)
            id = input("Input index of experiment to analyze")
            if id.isnumeric():
                restriction = self.select_expt_id(int(id))
                self._expt_id = int(id)
            else:
                return False
        if not restriction:
            return False
        else:
            return restriction

    def import_datajoint_data(self, *args, **kwargs):
        '''
        Main routine of using a restriction to extract data from DataJoint.  First the restriction is confirmed, then data retrieved
        if there is already a restriction in place, but you want to change it, or no restrction is present then args, kwargs or
        interactive input can be used
        :param args: Typically and index for a entry in self.expts
        :param kwargs: "Typically entries for making the selection from the database directly.
        :return:
        '''
        self.restriction = self.get_restriction(*args, **kwargs)
        if not self.restriction:
            print("Problem with restriction returning without loading data")
            return False
        self.exp_dict["restriction"] = self.restriction
        if not self.silence: print("import using restrictions ", self._expt_id, self.restriction)

        if "fluorescence" in self.autoimport:
            try:
                self.frame_times, self.raw_fluorescence = self.get_fluorescence()
                if self.frame_times.shape[0] != self.raw_fluorescence.shape[0]:
                    print("got uneven frame_times and responses", self.frame_times.shape, self.raw_fluorescence.shape)
                    print("Trimming back from end to make even")
                    if self.frame_times.shape[0] > self.raw_fluorescence.shape[0]:
                        self.frame_times = self.frame_times[:self.raw_fluorescence.shape[0]]
                    else:
                        self.raw_fluorescence = self.raw_fluorescence[self.frame_times.shape[0],:]
            except:
                print("Failure importing all fluorescence and timing data")
                print("Halting import early")
                return
        if "odor_trials" in self.autoimport:
            try:
                self.pd_odor_trials, self.odor_info, self.nozzle_to_odor_map = self.get_odors()
            except:
                print("Failure importing odor_trial data")
        if "treadmill" in self.autoimport:

            tread_time, tread_velocity = self.get_behavior(copy.deepcopy(self.frame_times))
            if not tread_time:
                print("Failure importing treadmill data")
        if "respiration" in self.autoimport:
            try:
                self.raw_resp_time, self.raw_resp_trace = self.get_respiration()
            except:
                print("Failure importing respiration data")

        # Add units onto data
        if not self.silence: print("importing data type ", self.import_data_type)
        if "fluorescence" in self.available_data:
            if self.import_data_type == "df_f" and "fluorescence" in self.available_data:
                self.raw_fluorescence = pq.Quantity(self.raw_fluorescence, DataMan.df_f_unit)
            elif self.import_data_type == "rfu":
                self.raw_fluorescence = pq.Quantity(self.raw_fluorescence, DataMan.rfu)

        if "frame_times" in self.available_data:
            if self.start_at_zero and "frame_times" in self.available_data:
                self.start_time_shift = -self.exp_signal_start_time
            else:
                self.start_time_shift = pq.Quantity(0.0, self.time_unit)
            self.times = pq.Quantity(self.frame_times, self.time_unit) + self.start_time_shift
            if not self.silence: print("timepoints before after shift: ",self.frame_times.shape, self.times.shape, " shifted time: ", self.frame_times[0], self.frame_times[-1], self.times[0], self.times[-1])
            self.avg_dt = (self.times[-1] - self.times[0]) / float(self.times.shape[0]-1)
        if "tread_time" in self.available_data:
            tread_time = pq.Quantity(tread_time, self.time_unit) + self.start_time_shift
            pos_ttime = np.where(tread_time >= 0.0)[0]
            self.raw_tread_time = tread_time[pos_ttime]
            self.raw_tread_velocity = tread_velocity[pos_ttime]
            self.prepare_treadmill(**kwargs)
        if "resp_time" in self.available_data:
            self.raw_resp_time = pq.Quantity(self.raw_resp_time, self.time_unit) + self.start_time_shift
            self.prepare_respiration(**kwargs)
        #print("imported time", self.times[0], self.times[-1], self.times.shape[0])
        if "odor_trials" in self.available_data:
            self.pd_odor_trials["trial_start_time"] = self.pd_odor_trials["trial_start_time"] + float(self.start_time_shift)
            self.pd_odor_trials["trial_end_time"] = self.pd_odor_trials["trial_end_time"] + float(self.start_time_shift)
            late_stims = self.pd_odor_trials[self.pd_odor_trials["trial_end_time"] > float(self.times[-1])].index
            if len(late_stims) > 0:
                print("Dropping {0} stimuli after imaging ended".format(len(late_stims)))
                self.dropped_stims = self.pd_odor_trials.loc[late_stims, :]
            self.pd_odor_trials = self.pd_odor_trials.drop(late_stims)
        return True

    def get_fluorescence(self, stamp_start=20):
        if self.meso is None:
            self.connect_meso()
        if self.odor is None:
            self.connect_odor()

        if "frame_times" in self.available_data:
            if not self.silence: print("replacing uploaded fluorescence data")
            self.available_data.remove("fluorescence")
            self.available_data.remove("frame_times")
        responses = np.column_stack((self.meso.Fluorescence.Trace & self.restriction).fetch('trace'))
        try:
            frame_times = (self.odor.OdorSync & self.restriction).fetch1('frame_times')
            self.exp_datetime_stamp = (self.odor.OdorSync & self.restriction).fetch("sync_ts")[0]
            self.exp_signal_start_time = pq.Quantity((self.odor.OdorSync & self.restriction).fetch("signal_start_time")[0], self.time_unit)
            self.exp_recording_duration = pq.Quantity((self.odor.OdorSync & self.restriction).fetch("signal_duration")[0], self.time_unit)
        except:
            print("Timing data not available.  Faking it for the sake of examining the data")
            frame_times = np.arange(0, responses.shape[0]) * pq.Quantity(0.078, pq.s)
            self.exp_datetime_stamp = 0
            self.exp_signal_start_time = frame_times[0]
            self.exp_recording_duration = frame_times[-1]
            raise BADDATAJOINT
        if not self.silence: print("unadjusted frame times, shapes", frame_times[0],frame_times[-1], frame_times.shape, responses.shape)
        if stamp_start:
            responses[0:stamp_start,:] = responses[stamp_start:2*stamp_start,:]
        if not self.silence: print("Fluorescence data Downloaded")
        self.available_data.append("fluorescence")
        self.available_data.append("frame_times")

        return frame_times, responses

    def get_odors(self):
        if self.odor is None:
            self.connect_odor()
        if "odor_trials" in self.available_data:
            if not self.silence: print("replacing uploaded odor trial data")
            self.available_data.remove("odor_trials")
        pd_odor_trials = pd.DataFrame(((self.odor.OdorTrials & self.restriction) * self.odor.OdorConfig).fetch(order_by='channel'))
        n_stims = pd_odor_trials["trial_idx"].shape[0]
        if self.label_conc:
            for i in range(n_stims):
                indx = pd_odor_trials.index[i]
                conc_str = "_{0}".format(pd_odor_trials.loc[indx, "concentration"])
                pd_odor_trials.loc[indx, "odorant"] = pd_odor_trials.loc[indx, "odorant"] + conc_str
                #print(pd_odor_tmp.loc[indx, "odorant"], conc_str)

        nozzle_to_odor_map = pd_odor_trials[["channel","odorant"]].set_index('channel').to_dict()["odorant"]
        odor_info = pd_odor_trials.drop_duplicates(["channel"])[["channel","odorant", "concentration",
                                                        "solution_date"]].set_index('channel').to_dict("index")
        if not self.silence: print("Odor info Downloaded")
        self.available_data.append("odor_trials")
        return pd_odor_trials, odor_info, nozzle_to_odor_map

    def get_behavior(self,image_scan_times=None):
        if image_scan_times is None:
            image_scan_times = self.times
        if self.stimulus is None:
            self.connect_stimulus()
        if self.treadmill is None:
            self.connect_treadmill()
        if "tread_time" in self.available_data:
            if not self.silence: print("replacing uploaded treadmill data")
            self.available_data.remove("tread_time")
            self.available_data.remove("tread_velocity")
        try:
            beh_scan_times = (self.stimulus.BehaviorSync & self.restriction).fetch1('frame_times')
            beh_tread_times, tread_vel = (self.treadmill.Treadmill & self.restriction).fetch1('treadmill_time', 'treadmill_vel')
        except:
            print("Problem getting treadmill data from Datajoint, returning False")
            return False, False

        if not self.silence: print("behavior fetched")
        # Step to remove nan values from data.  So far only seen in beh_tread_times
        gd = np.where(~np.isnan(beh_tread_times))[0]
        btt = beh_tread_times[gd]
        otv = tread_vel[gd]
        beh_odor_interp = interpolate.interp1d(beh_scan_times, image_scan_times, bounds_error=False)
        ott = beh_odor_interp(btt)
        # More Nans to remove
        gd2 = np.where(~np.isnan(ott))[0]
        ott = ott[gd2]
        otv = otv[gd2]
        # Also remove values occurring before imaging begins
        bht_mask = np.where(ott>image_scan_times[0])[0]

        if not self.silence: print("Treadmill data downloaded")
        self.available_data.append("tread_time")
        self.available_data.append("tread_velocity")
        return ott[bht_mask], otv[bht_mask]

    def get_respiration(self):
        if self.odor is None:
            self.connect_odor()
        if "resp_time" in self.available_data:
            if not self.silence: print("replacing uploaded respiration data")
            self.available_data.remove("resp_time")
            self.available_data.remove("resp_trace")
        try:
            resp_trace, odor_resp_times = (self.odor.Respiration & self.restriction).fetch1('trace', 'times')
        except:
            print("Unable to get_respiration.  Returning False")
            return False, False
        if not self.silence: print("Respiration data downloaded")
        self.available_data.append("resp_time")
        self.available_data.append("resp_trace")

        return odor_resp_times, resp_trace

    def prepare_respiration(self, force=False, remove_offset=False, **kwargs):
        '''
        Method to process respiration data to smooth, reduce size and ID inhalation and exhalation periods
        :param remove_offset: if the data is uploaded from csv there is often a time offset that was already removed
                from the times data.  If positive this will remove the same offset from respiration data
        :param kwargs:
        :return:  inserts results into self.exp_time["respiration"]
        '''
        offset = pq.Quantity(0.0, self.time_unit)
        if "respiration" not in self.exp_dict.keys():
            self.exp_dict["respiration"] = {}
        keys = ["y_red", "delta", "sensitivity", "f_lower", "f_upper"]
        if "resp_time" not in self.available_data or force:
            # Returns False if fetch of respiration data fails
            if not self.silence: print("Fetching respiration data")
            raw_resp_time, raw_resp_trace = self.get_respiration()
            if not self.silence: print("Respiration data fetched")
            if raw_resp_time is not False:
                if remove_offset:
                    offset = self.start_time_offset
                raw_resp_time = pq.Quantity(raw_resp_time, self.time_unit) + self.start_time_shift + offset
                t = raw_resp_time
                y = raw_resp_trace
            else:  # If can't get respiration data, then nothing to prepare
                return
        elif "respiration" in self.exp_dict.keys():
            t = self.exp_dict["respiration"]["raw_time"]
            y = self.exp_dict["respiration"]["raw_trace"]

        else:
            print("Unable to access respiration")
            return
        dt = (t[-1] - t[0]) / float(len(t) - 1)

        for key, value in kwargs.items():
            if key in keys:
                self.resp_processing[key] = value

        y_red = self.resp_processing["y_red"]
        delta = self.resp_processing["delta"]
        sensitivity = self.resp_processing["sensitivity"]
        f_lower = self.resp_processing["f_lower"]
        f_upper = self.resp_processing["f_upper"]
        if isinstance(y, (pq.Quantity)) and not isinstance(sensitivity, (pq.Quantity)):
            unit = y.units
            sens = make_quantity(sensitivity, unit)
        else:
            sens = sensitivity
        if f_lower == "raw" and f_upper == "raw":
            bp_y = y
        elif f_lower == "raw":
            bp_y = smooth(y, corner_f=(f_upper), dt=dt)
        elif f_upper == "raw":
            bp_y = high_pass(y, corner_f=(f_upper), dt=dt)
        else:
            bp_y = band_pass(y, corner_f=(f_lower, f_upper), dt=dt)
        if y_red:
            resp_t, resp_y = data_reduce(t, bp_y, y_red)
        else:
            resp_t = t
            resp_y = bp_y
        breath = np.zeros_like(resp_y)

        for i in range(delta, len(resp_y) - delta):
            slope = resp_y[i + delta] - resp_y[i - delta]
            if slope > sens:
                breath[i] = 1
            elif slope < -sens:
                breath[i] = -1
            else:
                breath[i] = 0

        resp_time = resp_t
        resp_trace = resp_y
        breaths = breath
        binterp = interpolate.interp1d(resp_time, breaths, kind="nearest" )
        # Give frametimes where animal was inhaling the odor.  Same number of points as self.times
        resp_phase = binterp(self.times)
        self.exp_dict["respiration"] = {"raw_time": raw_resp_time, "raw_trace": raw_resp_trace, "time": resp_time, "trace": resp_trace, "breaths": breath,
                                        "phases": resp_phase, "processing": self.resp_processing}
        self.odor_inhale_offset()


    def odor_inhale_offset(self, force=True):
        '''
        Determine how many timepoints elapsed after odor onset before an inhalation began
        :param force: used to prevent an infinite loop if method is run from within prepare_respiration
        :return: self.odor_inhale_delays dictionary
        '''
        if not force:
            if "resp_time" not in self.available_data:
                self.prepare_respiration()
        odor_offsets = []

        for i in self.pd_odor.index:
            start = self.pd_odor.at[i, "start"]
            count = 0
            resp_phase = self.exp_dict["respiration"]["phases"]
            phase = resp_phase[start]
            # Counting up frames from odor start until an inhalation takes place
            while phase < 1:
                count += 1
                phase = resp_phase[start + count]
            odor_offsets.append(count)
        self.pd_odor = self.pd_odor.assign(delay_to_inhale=odor_offsets)


    def prepare_treadmill(self, remove_offset=False, **kwargs):
        '''
        Method to process respiration data to smooth, reduce size and ID inhalation and exhalation periods
        :param kwargs:
        :param remove_offset: if the data is uploaded from csv there is often a time offset that was already removed
                from the times data.  If positive this will remove the same offset from treadmill data
        :return:  inserts results into self.tread_time and self.tread_velocity
        '''
        offset = pq.Quantity(0.0, self.time_unit)
        if "treadmill" not in self.exp_dict.keys():
            self.exp_dict["treadmill"] = {}
        keys = ["y_red", "delta", "sensitivity", "f_lower", "f_upper"]
        if "tread_time" not in self.available_data:
            # Returns False if fetch of respiration data fails
            tread_time, tread_velocity = self.get_behavior()
        if "tread_time" in self.available_data:
            if remove_offset:
                offset = self.start_time_offset
            tread_time = pq.Quantity(tread_time, self.time_unit) + self.start_time_shift + offset
            pos_ttime = np.where(tread_time >= 0.0)[0]
            raw_tread_time = tread_time[pos_ttime]
            raw_tread_velocity = tread_velocity[pos_ttime]
        else: # If can't get treadmill data, then nothing to prepare
            return

        t = raw_tread_time
        y = raw_tread_velocity

        dt = (t[-1] - t[0]) / float(len(t) - 1)

        # Treadmill processing data can be entered by uinit.init file or handed over during build or run of method
        if "tread_processing" in kwargs.keys():
            self.tread_processing = kwargs["tread_processing"]
        else:
            for key, value in kwargs.items():
                if key in keys:
                    self.tread_processing[key] = value

        y_red = self.tread_processing["y_red"]

        f_lower = self.tread_processing["f_lower"]
        f_upper = self.tread_processing["f_upper"]

        if f_lower == "raw" and f_upper == "raw":
            bp_y = y
        elif f_lower == "raw":
            bp_y = smooth(y, corner_f=(f_upper), dt=dt)
        elif f_upper == "raw":
            bp_y = high_pass(y, corner_f=(f_upper), dt=dt)
        else:
            bp_y = band_pass(y, corner_f=(f_lower, f_upper), dt=dt)
        if y_red:
            tread_t, tread_y = data_reduce(t, bp_y, y_red)
        else:
            tread_t = t
            tread_y = bp_y

        self.tread_time = tread_t
        self.tread_velocity = tread_y
        self.exp_dict["treadmill"] = {"raw_time": raw_tread_time, "raw_velocity": raw_tread_velocity, "time": tread_t, "velocity": tread_y,
                                        "processing": self.tread_processing}

    def process_odors(self):
        self.pd_odor_trials = self.pd_odor_trials.sort_values(by=["trial_start_time", "channel"])
        self.odor_to_nozzle_map = matching_noz(self.nozzle_to_odor_map)
        self.pd_odor = self.map_odortrials_to_frames_and_combine(self.pd_odor_trials)
        odor_labels = pd.unique(self.pd_odor["odorant"])
        odor_labels_dict = {}
        self.n_unique_odor_combos = len(odor_labels)
        n_unique = self.n_unique_odor_combos + 1
        for odorant in odor_labels:
            # Tags loaded odors first by their entry in the odor_to_nozzle_map
            # For same odor coming from multiple nozzles, only the lowest nozzle id is used
            # Actual nozzle id for a trial are in self.pd_odor["nozzles']
            odorid = [v[0] for k, v in self.odor_to_nozzle_map.items() if k == odorant]
            if len(odorid) > 0:
                odorid = min(odorid)
            else:
                # if it is not a loaded odor, then it must be a odor combo of 2 or more nozzles
                odorid = n_unique
                n_unique += 1
            odor_labels_dict[odorid] = odorant
        self.odor_labels = {}
        odor_id = 1
        # Sort odors by n_unique value, then label from 1 to len(odor_labels_dict)
        for n, odorant in sorted(odor_labels_dict.items()):
            self.odor_labels[odor_id] = odorant
            odor_id += 1

    def import_csv_response_data(self):
        '''
        Builds the DataMan object using csv files.  Does not establish a connection to the DataJoint db
        :return: loads into object: self.: data_type, times, responses, avg_dt, glomeruli, odors, raw_stim
        '''
        #import from csv files
        # Main data import rout ine brings in recording time and odor responses of the set of glomeruli of interest
        times, responses = import_csv_data(self.fdata, self.use_glomeruli, path=self.my_paths.data_path)
        if self.start_at_zero:
            self.start_time_shift = pq.Quantity(- times[0], self.time_unit)
        else:
            self.start_time_shift = pq.Quantity(0.0, self.time_unit)
        self.times = pq.Quantity(times , self.time_unit) + self.start_time_shift

        self.avg_dt = pq.Quantity((times[-1] - times[0]) / float(times.shape[0]-1), self.time_unit)
        # Remove artifact at start of traces
        responses[0:20,:] = responses[20:40,:]
        #self.times = pq.Quantity(times, self.time_unit) - self.start_time
        #self.avg_dt = pq.Quantity(avg_dt, self.time_unit)

        if self.import_data_type == "df_f":
            self.raw_fluorescence = pq.Quantity(responses, DataMan.df_f_unit)
        elif self.import_data_type == "rfu":
            self.raw_fluorescence = pq.Quantity(responses, DataMan.rfu)

        if not self.silence: print("analyzing data type ", self.import_data_type)

        # Converts "all" to a list of glomeruli assuming starting from 1 to n
        if self.use_glomeruli == "all":
            self.glomeruli = [i for i in range(1, self.raw_fluorescence.shape[1]+1)]
        else:
            self.glomeruli = self.use_glomeruli

    def import_csv_odors(self):
        if isinstance(self.fodors, (str)):
            odor_import = self.my_paths.data_path + "\\" + self.fodors
            l_str_odors = list(csv.reader(open(odor_import)))
        else:
            l_str_odors = list(fodors)
            # convert IDs to int and remove any leading or trailing spaces from odor names
        self.odor_labels = {int(id): odor.strip() for id, odor in l_str_odors}
        self.n_unique_odor_combos = len(self.odor_labels)

    def import_csv_stimuli(self):
        '''
        Imports the raw stimulus info with odor ID, start, stop for each puff
        Needs to have odor info imported first
        :param self.fstim:  file string for stimuli to use of a list of odor stimuli
        :param self.use_odors:  output from active_odors that IDs odors for study
        :return:  all stimuli calls that have the odors in use_odors
        '''
        stim_file = self.my_paths.data_path + "\\" + self.fstim
        pd_odor = pd.read_csv(stim_file)
        time_shift = self.start_time_shift + self.start_time_offset
        pd_odor.loc[:,"trial_start_time"] = pd_odor.loc[:,"trial_start_time"] + time_shift
        pd_odor.loc[:,"trial_end_time"] = pd_odor.loc[:,"trial_end_time"] + time_shift
        for i in range(pd_odor.shape[0]):
            pd_odor.loc[i,"odorant"] = self.odor_labels[pd_odor.loc[i,"channel"]]
        return pd_odor

    def make_f0(self, f0_method=None, force=False):
        known_masks = self.known_masks
        return_f0s = self.f0s
        # if data is already df_f then nothing else to do, and label the f0_method
        if self.raw_fluorescence.units == DataMan.df_f_unit.units:
            print("Not raw fluorescence data, skipping f0 creation")
            f_f0_data = copy.deepcopy(self.raw_fluorescence)
            f0 = pq.Quantity(np.ones_like(f_f0_data), DataMan.rfu)
            return_f0s["sham"] = f0
            self.f0_method = "pre_processed"
            return return_f0s

        # If we have raw fluorescence rather than dF/F data input then it needs to be processed
        if f0_method is not None:
            if isinstance(f0_method, (str)):
                if f0_method == "combo":
                    f0s = self.combo
                elif f0_method in self.known_f0_methods:
                    f0s = [f0_method]
                elif f0_method in ("compare", "all", "test"):
                    f0s = self.known_f0_methods
            elif isinstance(f0_method, (list, tuple)):
                f0s = f0_method
                if "combo" in f0s:
                    f0s.remove("combo")
                    for f0 in self.combo:
                        if f0 not in f0s:
                            f0s.append(f0)
        else:
            f0s = self.known_f0_methods
            f0s.remove("combo")
        existing_f0s = [key for key in self.f0s.keys()]
        # Don't remake unless force is in play
        if force:
            print("forced creation of f0 ", f0s)
        else:
            f0s = [f0 for f0 in f0s if f0 not in existing_f0s]

        if not self.silence: print("making f0 from input rfu data", f0s)
        # Use regions outside areas of interest for background, so where mask == False
        if self.f0_mask in ["use_all", "unmasked"]:
            f0_masking = self.mask_slice_manager(how_to="unmasked")["mask"]
            interpolate = False
        elif self.f0_mask in known_masks:
            f0_masking= self.mask_slice_manager(how_to=self.f0_mask)["mask"]
            print("using mask", self.f0_mask, f0_masking.shape, f0_masking.sum())

            interpolate = True
        else:
            print("Desired mask is not in known_masks, using default")
            f0_masking = self.mask_slice_manager(how_to="odor_post")["mask"]
            interpolate = True
        for newf0 in f0s:
            if newf0 == "bkg_smooth":
                f_sdata = smooth(self.raw_fluorescence, self.corner_f, self.avg_dt, axis=0, silence=self.silence)
                kwargs = {"times": self.times, "sdata": f_sdata, "odor_mask": f0_masking, "interpolate": interpolate}
                f0 = bkg_smooth(kwargs, bkg_corner_f=self.bkg_smooth_corner_f, kind=self.bkg_smooth_interp, silence=self.silence)
                if not self.silence: print("f0 by bkg smooth", f0.shape)
                # NOTE: Consider whether better to use f_sdata or self.raw_fluorescence here
            elif newf0 == "low_pass":

                f0_gapped = smooth(self.raw_fluorescence[f0_masking==False,:], self.baseline_corner_f, self.avg_dt, axis=0, silence=self.silence)
                if not interpolate:
                    f0 = f0_gapped
                else:
                    f0 = interpolate_masked_regions(self.times, f0_gapped, f0_masking, kind="linear", axis=0)
            elif newf0 == "constant_avg":
                shape = self.raw_fluorescence.shape
                f0 = np.multiply(np.mean(self.raw_fluorescence[f0_masking==False,:], axis=0), np.ones([shape[0],shape[1]], dtype=np.float32))
            elif newf0 == "constant_med":
                shape = self.raw_fluorescence.shape
                f0 = np.multiply(np.median(self.raw_fluorescence[f0_masking==False,:], axis=0), np.ones([shape[0],shape[1]], dtype=np.float32))
            elif newf0 == "linear":

                f0_reg = linear_regression(self.times, self.raw_fluorescence, f0_masking)
                f0 = smooth(f0_reg, self.baseline_corner_f, self.avg_dt, axis=0, silence=self.silence)
            return_f0s[newf0] = make_quantity(f0, DataMan.rfu)
        return return_f0s

    def make_df_f(self, f0_method=None, force=False):
        '''
        Routine that does the actual df/f calculation by creating an f0 data and doing (f-f0)/f0
        Returns a dict of all available df_f with the new one added to the dict
        :param f0_method: Option for f0 to use in calculating df/f
        :return: the df/f data
        '''
        return_df_f = self.df_fs
        # Data already in correct form.
        if self.raw_fluorescence.units == DataMan.df_f_unit.units:
            return_df_f["sham"] = self.raw_fluorescence
            return return_df_f

        df_fs = [key for key in return_df_f.keys()]
        if force:
            if f0_method is not None:
                print("forced created of df_f ", f0_method)
                f0s = f0_method
            else:
                f0s = [key for key in self.f0s.keys()]
        else:
            f0s = [key for key in self.f0s.keys() if key not in df_fs]
        if self.f0_method == "combo" and "combo" not in f0s:
            if force:
                f0s.append("combo")
            elif "combo" not in df_fs:
                f0s.append("combo")
        for f0 in f0s:
            if f0 == "combo":
                f0_dynamic = self.f0s[self.combo[0]]
                f0_static = self.f0s[self.combo[1]]
                return_df_f["combo"] = make_quantity(np.divide(np.subtract(np.array(self.raw_fluorescence), np.array(f0_dynamic)), np.array(f0_static)), DataMan.df_f_unit)
            else:
                return_df_f[f0] = make_quantity(np.divide(np.subtract(np.array(self.raw_fluorescence), np.array(self.f0s[f0])), np.array(self.f0s[f0])), DataMan.df_f_unit)
        return return_df_f

    def make_response_data(self, region="pre_odor_post", force=False):
        '''
        Creates concatenated and overlap odor response data for plotting. If response data already exists, it will only remake if region has
        changed or force is True
        :return: dict entries for responses and odor application times for each odor for all glomeruli
         self.exp_dict["response_data"][anodor]["concat"]
         Keys: {"resp": response_data, "time": concat_resp_time, "stim_time": response_stim_time, "true_stim_time": concat_true_stim_time}
         self.exp_dict["response_data"][anodor]["overlap"]
         Keys: {"resp": olap_signal, "time": overlap_resp_time, "stim_time": olap_stim_time, "true_stim_time": olap_true_stim_time}
        '''
        if not self.silence: print("making concat and overlap response data")
        if "response_data" not in self.exp_dict.keys():
            self.exp_dict["response_data"] = {"region": region}
        elif self.exp_dict["response_data"]["region"] != region:
            self.exp_dict["response_data"] = {"region": region}
        elif force:
            self.exp_dict["response_data"] = {"region": region}
        else:
            return

        if self.plt_baselinesubtracted:
            baseline_method = self.exp_dict["pulse_baseline_method"]
        else:
            baseline_method = "no_baseline"
        for n, anodor in sorted(self.odor_labels.items()):
            if anodor not in self.exp_dict["response_data"].keys():
                self.exp_dict["response_data"][anodor] = {}
            t1 = time.time()
            concat_data, overlap_data = self.vector_isolate_responses(anodor, region=region, force=force)
            t2= time.time()
            self.exp_dict["response_data"][anodor]["concat"] = concat_data
            t3 = time.time()
            self.exp_dict["response_data"][anodor]["overlap"] = overlap_data
            t4 = time.time()
            if self.testing: print(anodor, t2-t1,t3-t2,t4-t3)

    def reload_from_pickle(self, file_name=None, store_path=None):
        '''
        TODO: Add ability to reload pickle from DataJoint
        :param file_name: if None will use self.pickle_file name unless that is None in which case it will open a PyQt file dialog
        :param store_path: a potential alternative storage path
        :return: DataMan rebuilt from stored pickle state
        '''
        if store_path is not None and store_path in ("datajoint_db", "dj_db"):
            print("Upload from DataJoint of pickle files not yet implemented")
            return
        # helpers.unpickle_it will activate a file dialog if no name or path provided
        if file_name is None and self.pickle_file is not None:
            file_name = self.pickle_file
        self.exp_dict = helpers.unpickle_it(file_name, store_path, self.my_paths)

    def save_as_pickle(self, store_path=None, filename="olfactory_exp_dict", prefix= None, compress=False):
        '''
        Packages up items in self.__dict__ into state_dict and includes state_dict in the pickling of self.exp_dict
        which contains all the critical values for working with the data
        :param store_path:
        :param filename:
        :param compress:
        :return:
        '''
        if prefix is None:
            prefix = '_'.join([str(v) for k, v in self.restriction.items()])
        obj_dict = self.__dict__
        if "state_dict" in self.exp_dict.keys():
            del(self.exp_dict["state_dict"])
        state_dict = {}
        kwargs = {"filename": filename, "my_paths": self.my_paths, "path": store_path, "compress":compress, "prefix": prefix}
        kwargs["path_ready"] = helpers.prepare_pickle_name(self.exp_dict,**kwargs)
        if not self.silence: print("saving pickle at ", kwargs["path_ready"])
        for key, value in obj_dict.items():
            if key in ("exp_dict", "my_paths", "ispolar", "old_state_dict"):
                pass
            elif key in self.dj_objects:
                pass
            elif value is None:
                state_dict[key] = "None"
            elif key in self.exp_dict.keys():
                state_dict[key] = "reload"
            else:
                state_dict[key] = value
        # Store a copy of the current state of object parameters in state_dict
        # All DataJoint objects have been removed
        state_dict["dj_objects"] = []
        state_dict["db_upload_available"] = False
        self.exp_dict["state_dict"] = state_dict
        # print("pickling", type(self.exp_dict), type(kwargs))
        helpers.pickle_exp_data(self.exp_dict, **kwargs)
        if "state_dict" in self.exp_dict.keys():
            self.exp_dict["old_state_dict"] = self.exp_dict.pop("state_dict")

    def id_candidates(self, by=""):
        '''
        Id's candidates whose integral responses exceed the criterion for significance
        :param by: whether analysis by odor or by glomerulus  by = "by_odor" or "by_glom"
        :return:
        '''
        cand_set = {}
        which_dict = "response_cand"
        if by != "":
            which_dict = which_dict + "_" + by
        cand_dicts = self.exp_dict["pd_candidates"][which_dict]
        odors = self.odor_labels
        cand_set["stats"] = cand_dicts["stats"]
        cand_set["positive"] = {}
        cand_set["negative"] = {}
        cand_set["summary"] = {}
        pos = cand_dicts["positive"]
        neg = cand_dicts["negative"]
        for odor in odors.values():
            pcands = pos[pos[odor].notnull()]
            ncands = neg[neg[odor].notnull()]
            pcand_ids = list(pcands.index)
            pcand_values = list(pcands[odor])
            cand_set["positive"][odor] = [(indx, val) for indx, val in zip(pcand_ids, pcand_values)]
            ncand_ids = list(ncands.index)
            ncand_values = list(ncands[odor])
            cand_set["negative"][odor] = [(indx, val) for indx, val in zip(ncand_ids, ncand_values)]
            cand_set["summary"][odor]= {"n_pos": len(pcand_ids), "n_neg": len(ncand_ids)}
            

        return cand_set

    def pd_cand(self, sig_prob=None):
        '''
        Method to build up the pd DataFrames that will hold the integrals and significances of the responses
        :param sig_prob: P value to use when determining if a response is significant
        :return:
        '''
        if sig_prob is None:
            sig_prob = self.sig_prob
        responses_dict = {}
        odor_key = [value for key, value in sorted(self.odor_labels.items())]
        df = pd.DataFrame(data= self.exp_dict["response_integral"].T, index=self.glomeruli, columns=odor_key)
        responses_dict["integrals"] = df
        responses_dict["response_cand"] = {}
        std = df.stack().std()
        mean = df.stack().mean()
        criterion = std * norm.ppf(sig_prob)
        if not self.silence: print("criterion to sig m/sd/crit", mean, std, criterion)
        pos = mean+criterion
        neg = mean - criterion
        pos_sig = df[df > pos]
        neg_sig = df[df < neg]
        if not self.silence: print(pos_sig.head(6))
        stats = {"sig_prob": self.sig_prob, "criterion": criterion, "mean": mean, "stdev": std, "pos_criterion": pos, "neg_criterion": neg}
        responses_dict["response_cand"] = {"positive": pos_sig, "negative": neg_sig, "stats": stats}
        
        
        std = df.std(axis=0)
        mean = df.mean(axis=0)
        criterion = std * norm.ppf(sig_prob)
        pos_sig = df[df>mean+criterion]
        neg_sig = df[df < mean - criterion]
        stats = {"sig_prob": self.sig_prob, "criterion": criterion, "mean": mean, "stdev": std, "pos_criterion": mean+criterion, "neg_criterion": mean-criterion}
        responses_dict["response_cand_by_odor"] = {"positive": pos_sig, "negative": neg_sig,  "stats": stats}
    
        std = df.std(axis=1)
        mean = df.mean(axis=1)
        criterion = std * norm.ppf(sig_prob)
        pos_sig = df[df>mean+criterion]
        neg_sig = df[df < mean - criterion]
        stats = {"sig_prob": self.sig_prob, "criterion": criterion, "mean": mean, "stdev": std, "pos_criterion": mean+criterion, "neg_criterion": mean-criterion}
        responses_dict["response_cand_by_glom"] = {"positive": pos_sig, "negative": neg_sig,  "stats": stats}
        
        return responses_dict

    def adj_criterion(self, sig_prob):
        '''
        Changes probability for id'ing significant responses and re-evaluates the candidates that pass criterion using id_candidates() method
        :param criterion:
        :return:
        '''
        if sig_prob != self.sig_prob:
            self.sig_prob = sig_prob
            self.exp_dict["pd_candidates"] = self.pd_cand(self.sig_prob)
            self.exp_dict["response_cand"] = self.id_candidates()

    def re_sort(self, sort_by):
        if sort_by not in ("duration", "time", "odorant"):
            print("Do not understand sort_by: ", sort_by)
            return
        if self.sort_by == sort_by:
            print("Already sorted by ", sort_by)
            return
        self.sort_by =  sort_by
        if self.sort_by == "duration":
            sort_by = ["duration", "trial_start_time"]
        elif self.sort_by == "time":
            sort_by = ["trial_start_time", "duration"]
        elif self.sort_by == "odorant":
            sort_by = ["odorant", "trial_start_time"]
        # TODO: Fix this call
        self.pd_odor =self.pd_odor.sort_values(by=self.sort_by)

        #self.exp_dict["slices"] = self.generate_slices(times=self.exp_dict["times"], stimuli=self.stimuli, pre_buffer=self.pre_buffer, buffer=self.buffer)
        self.make_response_data(self.response_mask, force=True)

    def simple_integral(self):
        '''
        Approximates an "integral" by summing responses above baseline and dividing by total time of stimulus.
        NOTE: Code is written to allow easy modification to save intermediate net_signal and net_time data
        :return: self.exp_dict["response_integral"]
        '''
        max_key = max(self.odor_labels.keys())
        response_sum = np.zeros([max_key,self.exp_dict["sdata-baseline"].shape[1]])

        for n, odor in sorted(self.odor_labels.items()):
            masking = self.mask_slice_manager(how_to=self.response_mask, by_odor=odor)
            resp_mask = masking["mask"]
            resp_slices = masking["slices"]

            if self.scale_integrals == "n_trials":
                scale = float(len(resp_slices))
            elif self.scale_integrals == "stim_time":
                scale = masking["net_stim_time"]
                if self.testing: print("scale", scale)
            else:
                scale = float(self.scale_integrals)
            integral = np.sum(self.exp_dict["sdata-baseline"] * np.expand_dims((np.abs(resp_mask)), axis=1), axis=0)
            response_sum[n-1,:] = integral / scale

        return response_sum

    def refilter(self, corner_f):
        '''
        subroutine for changing filtering.  Can be called at any time since original raw data not affected
        :param corner_f:  corner frequency in Hz for filter or "raw" for unfiltered data
        :return:
        '''
        if corner_f == self.corner_f:
            return
        else:
            self.exp_dict["corner_f"] = corner_f
            self.corner_f = corner_f

        self.exp_dict["sdata"] = smooth(self.df_fs[self.f0_method], self.exp_dict["corner_f"], self.exp_dict["avg_dt"], axis=0)
        self.exp_dict["sdata-baseline"] = np.subtract(self.exp_dict["sdata"], self.exp_dict["global_baseline"])
        self.exp_dict["response_integral"] = self.simple_integral()
        self.exp_dict["pd_candidates"] = self.pd_cand(self.sig_prob)
        self.exp_dict["response_cand"] = self.id_candidates()

        self.make_response_data(self.response_mask, force=True)

    def vector_isolate_responses(self, odor_id, region="pre_maxodor_post", force=False):
        '''
        Isolates baseline subtracted slices of data in response to an odor
        :param odor_id: odor of interest
        :param sdata: smoothed data
        :param ls_data: globally baseline subtracted data
        :param d_times: delta times for the experiment
        :param slices: slices dictionary
        :return: sliced baseline subtracted signals and delta times for each slice
        '''
        # Select region to generate baseline based on baseline method
        points = self.exp_dict["sdata-baseline"].shape[0]
        segments = self.glomeruli
        pd_responses = pd.DataFrame(data=self.exp_dict["sdata-baseline"], index=np.arange(0,points,1), columns = segments)

        masking = self.mask_slice_manager(how_to=region, by_odor=odor_id, force=force)
        response_data = np.array(pd_responses[masking["mask"]==1])
        concat = {"resp": response_data, "time": masking["time"], "stim_time": masking["stim_time"], "true_stim_time":  masking["true_stim_time"]}
        pd_responses["true_time"] = self.times
        pd_responses["delta_ts"] = self.exp_dict["delta_ts"]
        olap_signal = []
        olap_stim_time = []
        olap_true_stim_time = []
        overlap_resp_time = []
        biggest_slice = 0
        for i, aslice in enumerate(masking["slices"]):
            indx = [i for i in range(aslice[0], aslice[1], 1)]
            sliced_at = pd_responses.loc[indx]
            biggest_slice = max(biggest_slice, len(indx))
            olap_signal.append(np.array(sliced_at.loc[:,segments]))
            stim_indx = masking["stim_time"].index.isin(indx)
            olap_stim_time.append(masking["stim_time"].loc[stim_indx])
            olap_true_stim_time.append(masking["true_stim_time"].loc[stim_indx])
            overlap_resp_time.append(masking["time"].loc[indx])
        overlap = {"resp": olap_signal, "time": overlap_resp_time, "stim_time": olap_stim_time, "true_stim_time": olap_true_stim_time}

        return concat, overlap


    def initial_fill_experiment_dict(self):
        '''
        Initial build of a dictionary of experiment data with responses segmented by odor used and masks to grab different parts
        of the data.  Will be built upon as the analysis grows to include results and saved plots
        :param times: time points for the data
        :param data: the dF/F data from the selected glomeruli
        :param stimuli: the stimuli we are working with
        :param dt: time step
        :param buffer: how many points to take before and after the stimulus
        :param glomeruli: index for the glomeruli in the data
        :param kwargs: Rebuilding of an exisitng exp_dict
        :return: dictionary of responses for each stimuli
        '''
        if not self.silence: print("Filling new exp_dict")
        exp_dict = {}
        exp_dict["buffer"] = self.buffer
        exp_dict["pre_buffer"] = self.pre_buffer
        #exp_dict["stimuli"] = copy.deepcopy(self.stimuli)
        #print("number of stimuli", len(self.stimuli))
        exp_dict["avg_dt"] = self.avg_dt
        exp_dict["baseline_corner_f"] = self.baseline_corner_f
        exp_dict["times"] = self.times
        exp_dict["delta_ts"] = self.delta_ts
        #exp_dict["glomeruli"] = self.glomeruli
        #exp_dict["odors"] = self.odors
        exp_dict["corner_f"] = self.corner_f
        exp_dict["pulse_baseline_method"] = self.pulse_baseline_method
        exp_dict["global_baseline_method"] = self.global_baseline_method
        exp_dict["global_baseline"] = self.global_baseline
        exp_dict["sdata"] = self.sdata
        exp_dict["precision"] = self.precision
        # Can extract indexes by filtering >=0
        #exp_dict["slices"] = copy.deepcopy(self.slices)
        odor_masking =  self.mask_slice_manager(how_to="odor_post")["mask"]
        exp_dict["no_odor_times"] = exp_dict["times"][odor_masking == False]

        exp_dict["sdata-baseline"] = np.subtract(exp_dict["sdata"], exp_dict["global_baseline"])

        exp_dict["no_odor-baseline"] = exp_dict["sdata-baseline"][odor_masking==False, :]

        return exp_dict

    def map_odortrials_to_frames_and_combine(self, raw_pd_odor):
        '''
        Generates abilty to slice data around odor trials several ways and to isolate slices for specific trials
        :param times: The frametime np.array in odor time
        :param stimuli: The odor stimuli after aggregation by combine_stimuli() as [odor_id, start, end, duration]
        :param buffer: Time to buffer after stimuli to acquire entire odor response
        :param pre_buffer: time block of data to acquire before odor stimuli.  if None, use buffer value
        :return: dictionary of slicing options both for all odors and _by_odor:
                    pre_stim_slice: slices taken before odor delivery (start-pre_buffer, start)
                    pre_odor_post_slice: slices over whole odor block (start-pre_buffer, end+buffer)
                    pre_maxodor_post_slice: slices over whole odor block (start-pre_buffer, maxend+buffer).
                    so all slices are of the same duration even if the odor deliveries are variable duration
                    odor_slice: slice taken during odor delivery: (start, end)
                    odor_post_slice: slice taken during odor response: (start, end+buffer)

        '''
        if not self.silence: print("Mapping Odor Trials")
        times = self.times
        #print("Using times", times[0], times[-1], len(times))
        self.pre_buffer_slice = int((self.pre_buffer / self.avg_dt) + 0.5)
        self.post_buffer_slice = int((self.buffer / self.avg_dt) + 0.5)
        self.common_slice_start = int((self.c_start / self.avg_dt) + 0.5)
        self.common_slice_end = int((self.c_end / self.avg_dt) + 0.5)
        pd_odor = combine_stim(raw_pd_odor, combine_common=True, print_combined=self.print_multiple_stimuli, sort_by=self.sort_by, silence=self.silence)
        if not self.silence: print("Repeated Stimuli Combined")
        pd_odor = pd_odor.assign(start='', stop='')
        n_slices = pd_odor.shape[0]
        time_start = time.time()
        pop_starts = []
        pop_ends = []
        dur = []
        st = time.time()
        for aslice in range(n_slices):
            start = pq.Quantity(pd_odor["trial_start_time"].iloc[aslice], self.time_unit)
            stop = pq.Quantity(pd_odor["trial_end_time"].iloc[aslice], self.time_unit)
            duration = stop - start
            #print("t1", time.time() - st)
            pop_start, pop_end = index_frames(times, start, stop, min_frame_time=self.min_frame_time, precision=self.precision)
            pop_starts.append(pop_start)
            # 1 is added so np.slicing is inclusive of the end point
            pop_ends.append(pop_end + 1)
            dur.append(pop_end - pop_start + 1)
            #print("t2", time.time() - st)

        pd_odor["start"] = pop_starts
        pd_odor["stop"] = pop_ends
        pd_odor["duration"] = dur
        self.max_stimulus_slice = max(dur)
        pd_odor = self.slice_maker(pd_odor)
        if not self.silence: print("slice timing", time.time() - st)
        return pd_odor

    def characterize_stims(self, tight=1.0):
        '''
        Identifies stimuli that are either overlapping or too close together to be treated independently.  Stimuli that can be treated
        independently are
        :param tight:
        :return:
        '''
        self.pd_odor["problems"] = ""
        pd_odor = copy.deepcopy(self.pd_odor.sort_values(by="trial_start_time"))
        n_stims = self.pd_odor.shape[0]
        n_unique_trials = len(self.odor_labels)
        n_st = pd_odor.columns.get_loc("start")
        n_stop = pd_odor.columns.get_loc("stop")
        n_prob = pd_odor.columns.get_loc("problems")
        for i in range(1, n_stims):
            if pd_odor.iloc[i, n_st] < pd_odor.iloc[i-1, n_stop]:
                pd_odor.iloc[i, n_prob] = "overlap"
            elif pd_odor.iloc[i, n_st] - self.post_buffer_slice < pd_odor.iloc[i-1, n_stop]:
                pd_odor.iloc[i, n_prob] = "tight"
            else:
                pd_odor.iloc[i, n_prob] = "ok"
        overlapping = pd_odor[pd_odor["problems"] == "overlap"].shape[0]
        tight = pd_odor[pd_odor["problems"] == "tight"].shape[0]
        ok = pd_odor[pd_odor["problems"] == "ok"].shape[0]

        self.exp_dict["quality_report"] = {"overlapping": overlapping, "tight": tight, "independent": ok, "unique": n_unique_trials}
        self.quality_report()

        if self.sort_by is not None:
            if self.sort_by == "time":
                pd_odor = self.pd_odor.sort_values(by="trial_start_time")
            elif self.sort_by == "duration":
                pd_odor = self.pd_odor.sort_values(by="duration")
        return pd_odor

    def quality_report(self):
        '''
        Summarizes the analysis quality of the odor trials and how many distinct trial variants were performed
        :return:
        '''
        print("Quality Report: assessment of independent unique stimuli in this experiment")
        if "restriction" in self.exp_dict.keys():
            print("Experiment {0} performed on animal {1} on date {2}".format(self.expt_id, self.exp_dict["restriction"]["animal_id"], self.exp_datetime_stamp))
        print("There are {0} overlapping stimuli in this experiment".format(self.exp_dict["quality_report"]["overlapping"]))
        print("There are {0} tight stimuli in this experiment".format(self.exp_dict["quality_report"]["tight"]))
        print("There are {0} independent stimuli in this experiment".format(self.exp_dict["quality_report"]["independent"]))
        print("There are {0} unique odor stimuli in this experiment".format(self.exp_dict["quality_report"]["unique"]))

    def slice_maker(self, pd_odor=None, force=False):
        print("Making slices")
        if pd_odor is None:
            pd_odor = self.pd_odor
        if "common" in pd_odor.columns and not force:
            return pd_odor
        elif "common" in pd_odor.columns and force:
            masks = copy.deepcopy(self.known_masks)
            masks.remove("unmasked")
            pd_odor.drop(masks)
        start_stop = pd_odor.loc[:,["start", "stop"]]
        com_s = start_stop.loc[:,"start"] + self.common_slice_start
        com_e = start_stop.loc[:,"start"] + self.common_slice_end
        com = pd.concat([com_s, com_e], axis=1)
        pd_odor["common"] = com.values.tolist()
        pre_s = start_stop.loc[:,"start"] - self.pre_buffer_slice
        pre = pd.concat([pre_s, start_stop["start"]], axis=1)
        pd_odor["pre_stim"] = pre.values.tolist()
        post = start_stop.loc[:,"stop"] + self.post_buffer_slice
        op = pd.concat([start_stop["start"], post], axis=1)
        pd_odor["odor_post"] = op.values.tolist()
        pre_s2 = start_stop.loc[:,"start"] - self.pre_buffer_slice
        mo = start_stop.loc[:,"stop"] + self.max_stimulus_slice + self.post_buffer_slice
        pmp = pd.concat([pre_s2, mo], axis=1)
        pd_odor["pre_maxodor_post"] = pmp.values.tolist()
        pre_s3 = start_stop.loc[:,"start"] - self.pre_buffer_slice
        pop = pd.concat([pre_s3, post], axis=1)
        pd_odor["pre_odor_post"] = pmp.values.tolist()
        pd_odor["odor"] = start_stop.values.tolist()
        mo2 = start_stop.loc[:,"stop"] + self.max_stimulus_slice + self.post_buffer_slice
        mp = pd.concat([start_stop["start"], mo2], axis=1)
        pd_odor["maxodor_post"] = mp.values.tolist()
        return pd_odor

    def mask_slice_manager(self, how_to="odor_post", by_odor=None, force=False):
        '''
        Generates useful data for specific odorant trials useful for statistical analyses or plotting
        Calculation done on as needed basis
        :param how_to: which mask to use.
        :param by_odor: which odorant, if None equivalent to "all_odors"
        :param force: Force remake of data before delivery.  Useful if buffer periods or filtering change
        :return: times, masks and slices.  Times are for mask and odor delivery.  Times as concat times and true times.
        '''
        if self.testing: print("retrieving mask", how_to, by_odor)
        odors = {**{len(self.odor_labels)+1: "all_odors"}, **self.odor_labels}
        if by_odor is None:
            by_odor = "all_odors"
        if isinstance(by_odor, (int)):
            try:
                by_odor = odors[by_odor]
            except:
                print("unknown odor_id", by_odor)
                return
        filter_mask_name = how_to + "_" + by_odor

        if filter_mask_name in self.exp_dict["data_subset"].keys() and not force:
            return self.exp_dict["data_subset"][filter_mask_name]

        if how_to not in self.known_masks:
            print("Unknown mask type", how_to)
            return
        if how_to == "unmasked":  # Don't bother storing just return False array and empty slices
            print("no_masking")
            mask_data = {"mask": np.zeros([self.times.shape[0]], dtype=np.bool), "slices": []}
            return mask_data
        # Grab the by_odor mask information and select out the how_to data and the timing data
        if how_to not in ["odor", "common"]:
            odor_pd_mask = copy.deepcopy(self.exp_dict["mask_dict"][by_odor].loc[:, ["time", "delta_ts", "odor", "common", how_to]])
        else:
            odor_pd_mask = copy.deepcopy(self.exp_dict["mask_dict"][by_odor].loc[:, ["time", "delta_ts", how_to]])
        # Filter down so only the True values in the how_to column remain
        tmp_pd_mask = odor_pd_mask[odor_pd_mask[how_to]==True]
        dt_list = list(tmp_pd_mask["delta_ts"])
        sum_times = []
        sum_time = 0
        for i in range(len(dt_list)):
            sum_time += dt_list[i]
            sum_times.append(sum_time)

        tmp_pd_mask = tmp_pd_mask.assign(concat_time = sum_times)
        if how_to =="common":
            odor_col = "common"
        else:
            odor_col = "odor"
        times = tmp_pd_mask["concat_time"]
        true_times = tmp_pd_mask["time"]
        np_mask = np.array(odor_pd_mask[how_to])
        # filter down so only the True values in the odor_col remain to get odor delivery info
        tmp_odor_mask = tmp_pd_mask[tmp_pd_mask[odor_col]==True]
        stim_times = tmp_odor_mask["concat_time"]
        total_stim = tmp_odor_mask["delta_ts"].sum()
        true_stim_times = tmp_odor_mask["time"]
        if by_odor == "all_odors":
            slice_pd = self.pd_odor
        else:
            slice_pd = self.pd_odor[self.pd_odor["odorant"] == by_odor]
        slice_out = list(slice_pd[how_to])
        mask_data = {"time": times, "true_time": true_times, "stim_time": stim_times, "true_stim_time": true_stim_times, "mask": np_mask, "slices": slice_out, "net_stim_time":total_stim}

        self.exp_dict["data_subset"][filter_mask_name] = mask_data
        return mask_data

    def mask_maker(self):
        '''
        Creates a collection of pandas DataFrames for masking fluorescence.  first column is times, second is delta times and
        the remaining columns are timepoints where odorant trials are occurring masked as indicated.
        :return: dictionary of the masks for each odorant
        '''
        additional_masks_to_make = {}
        self.pd_masks = pd.DataFrame(data=self.times, index=np.arange(0,self.times.shape[0],1), columns=["time"])
        delta_ts = self.times[1:] - self.times[:-1]
        self.delta_ts = np.concatenate([[0.0], delta_ts], axis=0)
        self.pd_masks["delta_ts"] = self.delta_ts

        self.exp_dict["mask_dict"] = {}
        additional_masks_to_make[len(self.odor_labels)+1] = "all_odors"
        mask_labels = {**additional_masks_to_make, **self.odor_labels}
        masks = self.known_masks
        masks.remove("unmasked")
        for n, n_odor in mask_labels.items():
            mask_0 = np.zeros([self.times.shape[0], len(masks)], dtype=np.bool)
            pd_mask = pd.DataFrame(data=mask_0, index=np.arange(0,self.times.shape[0],1), columns=masks)
            pd_masks = pd.concat([self.pd_masks, pd_mask], axis=1)
            if n_odor == "all_odors":
                tmp_pd_odor = copy.deepcopy(self.pd_odor)
            else:
                tmp_pd_odor = self.pd_odor[self.pd_odor["odorant"]==n_odor]
            start_stop = tmp_pd_odor.loc[:,"odor"]

            for i in range(start_stop.shape[0]):
                slicer = start_stop.iloc[i]
                pd_masks.loc[slice(slicer[0], slicer[1]),"odor"] = 1
                pd_masks.loc[slice(slicer[0] + self.common_slice_start, slicer[0] + self.common_slice_end),"common"] = 1
                pd_masks.loc[slice(slicer[0] - self.pre_buffer_slice, slicer[0]),"pre_stim"] = 1
                pd_masks.loc[slice(slicer[0], slicer[1] + self.post_buffer_slice),"odor_post"] = 1
                pd_masks.loc[slice(slicer[0] - self.pre_buffer_slice, slicer[1] + self.max_stimulus_slice + self.post_buffer_slice),"pre_maxodor_post"] = 1
                pd_masks.loc[slice(slicer[0] - self.pre_buffer_slice, slicer[1] + self.post_buffer_slice),"pre_odor_post"] = 1
                pd_masks.loc[slice(slicer[0], slicer[1] + self.max_stimulus_slice + self.post_buffer_slice),"maxodor_post"] = 1
            self.exp_dict["mask_dict"][n_odor] = pd_masks


    def reset_buffers(self, pre_buffer=None, buffer=None):
        '''
        Routine to change amount of buffering before and after odor responses
        and remake slices, masks and concat response data
        :param pre_buffer: time to include before odor delivery
        :param buffer: time to include after odor delivery ends
        :return: parent object is updated in self.exp_dict for slices and masks and concat plotting data is remade
        '''
        if pre_buffer is not None:
            self.exp_dict["pre_buffer"] = pq.Quantity(pre_buffer, self.time_unit)
        if buffer is not None:
            self.exp_dict["buffer"] = pq.Quantity(buffer, self.time_unit)

        self.exp_dict["mask_dict"] = {}
        self.mask_maker()
        self.make_response_data(self.response_mask, force=True)
    # TODO: Redo
    def make_sham_responses_dict(self, region="common", n_shams=20):
        '''
        Creates a dictionary of fake responses to a given odor that match the duration of the actual odor deliveries
        for that odor but are taken from the data outside odor delivery times
        :param self:
        :return:
        '''
        sham_resps = {}
        for key, odor in self.odor_labels.items():
            print(type(odor), type(n_shams), type(region))
            self.sham_masks_dict[odor], self.sham_slices_dict[odor] = stats.match_stimuli(self, [key, odor], n_shams, mask_region="odor_post", match_region=region)
        for odor in self.odor_labels.values():
            sham_resps[odor] = {}
            sham_odor_slices = self.sham_slices_dict[odor]
            for key, value in sham_odor_slices.items():
                sham_resps[odor][key] = resp_dict(self, value["odor"], value["response"], sham=True)
        self.sham_resp_dict = sham_resps
    # TODO: Redo
    def make_response_dict(self, region="odor_post"):
        region = region + "_slice_by_odor"
        for odor in self.odor_labels.values():
            if region == "common":
                odor_stim = self.exp_dict["slices_dict"]["common_slice_by_odor"][odor]
                signal_region = self.exp_dict["slices_dict"][region][odor]
            else:
                odor_stim = self.exp_dict["slices_dict"]["odor_slice_by_odor"][odor]
                signal_region = self.exp_dict["slices_dict"][region][odor]

            self.response_dict[odor] = resp_dict(self, odor_stim, signal_region)

    def get_response_data(self, odor, duration="all"):
        if isinstance(odor, (int)):
            odor = self.odor_labels[odor]
        pd_resp = self.pd_dict[odor]
        if not duration == "all":
            if isinstance(duration, (str, int, float)):
                duration = [duration]
            q_list = ["duration == " + str(i) for i in duration]
            my_query = " or ".join(q_list)
            pd_resp = pd_resp.query(my_query)
        return pd_resp
    # TODO: Redo
    def make_pandas_dict(self, region="odor_post", sensitivity=0.25, sig_fig=1):

        pd_dict = {}
        for key, value in self.response_dict.items():
            n = len(value)
            shape = value[0]["integral"].shape[0]
            col_lab = ["duration"]
            col_lab.extend([i for i in range(1, shape+1)])
            index_labels = [j for j in range(n)]
            pd_data = pd.DataFrame(index=index_labels, columns=col_lab)     #, dtype=pd.float32)
            durs = []
            for i, results in value.items():
                durs.append(float(results["duration"]))
            durs = find_modes(durs, sensitivity=sensitivity, sig_figs=sig_fig)
            for i, results in value.items():
                if "common" in region:   # If common don't separate responses by stimulus duration
                    dlabel = [1.0]
                else:
                    dlabel = [d for d in durs if abs(d-float(results["duration"]))<sensitivity]
                pd_data.loc[i,'duration'] = dlabel[0]
                #print(pd_data.iloc[i,1:].shape, results["integral"].shape)
                pd_data.iloc[i,1:] = results["integral"]
            pd_dict[key] = pd_data
        self.pd_dict = pd_dict

    # TODO: Redo
    def make_sham_pandas(self, odor, region="odor_post", sensitivity=0.25, sig_fig=1):
        '''
        Converts the sham_response dict into a pandas DataFrames
        :param odor: odor_id to select which sham_responses to analyze
        :param durs: whether and how to separate by duration of the stimulus if "common" treat all odor deliveries as the same duration
        :return: a dictionary of sham response panda DataFrames matching the stimuli for this odor
        '''

        pd_dict = {}
        for key, value in self.sham_resp_dict[odor].items():
            n_shams = len(value)
            shape = value[0]["integral"].shape[0]
            col_lab = ["duration"]
            col_lab.extend([i for i in range(1, shape+1)])
            index_labels = [j for j in range(n_shams)]
            pd_data = pd.DataFrame(index=index_labels, columns=col_lab, dtype=float)     #, dtype=pd.float32)
            durs = []
            for i, results in value.items():
                durs.append(float(results["duration"]))
            durs = find_modes(durs, sensitivity=sensitivity, sig_figs=sig_fig)
            for i, results in value.items():
                if "common" in region:   # If common don't separate responses by stimulus duration
                    dlabel = [1.0]
                else:
                    dlabel = [d for d in durs if abs(d-float(results["duration"]))<sensitivity]
                pd_data.loc[i,'duration'] = dlabel[0]
                #print(pd_data.iloc[i,1:].shape, results["integral"].shape)
                pd_data.iloc[i,1:] = results["integral"]
            pd_dict[key] = pd_data
        return pd_dict




def combine_stim(pd_odor, combine_common=True, print_combined=False, sort_by=None, min_gap=0.2, silence=False):
    '''
    Combines consecutive odor puffs into one stimulus for stimuli in stims.
    If next odor puff of same odor occurs within 2 x min_gap they are combined into one
    odor delivery.
    :param stims: file with stimuli for this odor(odorID, starttime, endtime)
    :param dt: average time step for data
    :param sort_by: change to "d" or "duration" to sort by stim length
    :param t_adj:  Use to adjust the start time of these data to match other start times.
    :return: sorted stimuli with (odorID, start, end, duration) tuple in list
    '''
    start_time = time.time()
    empty_list = ["Empty", "Mineral Oil", "Blank", "Empty_0.00", "Mineral Oil_0.00", "Blank_0.00"]
    # Sorts stims by start time to combine if consecutive stims are the same odorant
    pd_odor_tmp = copy.deepcopy(pd_odor.sort_values(by=["channel","trial_start_time"]))
    n_trials = pd_odor["trial_idx"].max() + 1
    pd_odor_tmp["to_del"] = 0
    pd_odor_tmp["nozzles"] = ''
    pd_odor_tmp["duration"] = 1
    pd_odor_tmp["dt"] = 0.0
    # This loop examines the different odor trials to combine multi_odor stimuli into a single stimulus
    # These combos are reinserted into the database as a single line entry and the other lines for
    # The same stimulus are deleted
    for i in range(0,n_trials):
        active_trial = pd_odor_tmp[pd_odor_tmp["trial_idx"] == i]
        if len(active_trial) == 0: continue
        odors_in_trial = active_trial.shape[0]
        indxs = active_trial.index
        nozzles = sorted([n for n in list(active_trial["channel"])])
        if odors_in_trial > 1:
            if combine_common:
                unique_odors = list(set(active_trial["odorant"]))
                if len(unique_odors) > 1:
                    unique_odors = [odr for odr in unique_odors if odr not in empty_list]
            if len(unique_odors) > 1:
                odor_combo = "|".join(unique_odors)
            else:
                odor_combo = list(unique_odors)[0]
            pd_odor_tmp.loc[indxs[0],"odorant"] = odor_combo
            pd_odor_tmp.at[indxs[0],"nozzles"] = nozzles
            if print_combined:
                print("Multiple_stimuli", i, pd_odor_tmp[pd_odor_tmp.index == indxs[0]]["odorant"], indxs)

            pd_odor_tmp = pd_odor_tmp.drop(indxs[1:])
        else:
            pd_odor_tmp.at[indxs[0],"nozzles"] = nozzles

    time_loop_1 = time.time() - start_time
    if not silence: print("time loop 1 ", time_loop_1)
    pd_odor_tmp = pd_odor_tmp.sort_values(by="trial_start_time")
    unique_stims = list(set(pd_odor_tmp["odorant"]))
    pd_odor_tmp["dt"] = 0.0
    pd_odor_tmp["to_del"] = 0
    #n_dt = pd_odor_tmp.columns.get_loc("dt")
    #n_st = pd_odor_tmp.columns.get_loc("trial_start_time")
    n_et = pd_odor_tmp.columns.get_loc("trial_end_time")
    n_td = pd_odor_tmp.columns.get_loc("to_del")
    criterion = 2. * min_gap
    cum_time_start = time.time()

    for odor in unique_stims:
        active_trials = copy.deepcopy(pd_odor_tmp[pd_odor_tmp["odorant"] == odor])
        n_o_stims = active_trials.shape[0]
        ndx_id = np.arange(0,n_o_stims)
        active_trials.loc[:,"ndx"] = ndx_id
        start = active_trials["trial_start_time"].to_numpy()
        end = active_trials["trial_end_time"].to_numpy()
        dt = np.zeros(n_o_stims)
        dt[1:] = start[1:] - end[:-1]
        active_trials.loc[:,"dt"] = dt

        new_stims = active_trials[active_trials["dt"] > criterion]["ndx"].to_numpy(dtype=np.int)
        size_new_stims = new_stims.shape[0]
        start_stims = np.zeros(size_new_stims + 1, dtype=np.int)
        end_stims = np.ones(size_new_stims + 1, dtype=np.int) * (n_o_stims -1)
        start_stims[1:] = new_stims
        end_stims[:-1] = new_stims - 1
        end_info = list(active_trials.iloc[end_stims, n_et])
        active_trials.iloc[start_stims, n_et] = end_info
        active_trials.iloc[start_stims, n_td] = 1

        pd_odor_tmp.loc[active_trials.index] = active_trials.loc[active_trials.index]
    pd_return = pd_odor_tmp.drop(pd_odor_tmp[pd_odor_tmp["to_del"] == 0].index)
    if not silence: print("for loop cum time", time.time() - cum_time_start)
    pd_return = pd_return.drop(columns=["to_del", "animal_id", "odor_session", "recording_idx", "channel", "concentration", "solution_date", "dt"])
    if sort_by is not None:
        if sort_by == "time":
            pd_return = pd_return.sort_values(by="trial_start_time")
        elif sort_by == "duration":
            pd_return = pd_return.sort_values(by="duration")
    return pd_return

def import_csv_data(fdata="glom5_19.csv", glomeruli="all", path=None):
    '''
    Takes a data file in csv format and extract traces for the selected glomeruli
    :param fdata: data file name
    :param glomeruli: a list of IDs if restricting glomeruli
    :param start_time: an offset time to remove to align odor delivery and responses
    :return: times array, responses array
    '''
    if path is None:
        load_csv = fdata
    else:
        load_csv = path + "\\" + fdata
    rdata = np.loadtxt(load_csv, dtype=float, delimiter=",")

    x_data = rdata[:,0]
    if isinstance(glomeruli, (str)):
        y_data = rdata[:,1:]
    else:
        y_data = rdata[:,glomeruli]
    return x_data, y_data

def bkg_smooth(kwargs, bkg_corner_f=0.1, kind="linear", interpolate=True, silence=False):
    '''
    This routine creates a smooth background from parts of the y data that do not have stimulus responses
    typically it will receive a times, sdata and mask.
    :param key_name: what to call the new background function
    :param bkg_corner_f:  corner frequency for smoothing the background data
    :param kind:  method of interpolation to use with scipy.interpolate.interp1d
    :param kwargs: typically a dict that has x, y, and a odor_post_mask mask as input Also override for interpolate if needed

    :return: kwarg dict with calculated baseline as key_name
    '''
    print("bkg_smooth", kwargs["times"].shape, kwargs["sdata"].shape, kwargs["odor_mask"].shape, kwargs["odor_mask"].sum())

    x = kwargs["times"]
    y = kwargs["sdata"]
    bkg = kwargs["odor_mask"]
    ybkg = y[bkg==False,:]
   # Smoothed to connect pre and post stimulus points
    sybkg = smooth(ybkg,bkg_corner_f, silence=silence)
    #Interpolate across stimuli using
    if "interpolate" in kwargs.keys():
        interpolate = kwargs["interpolate"]

    if interpolate:
        if ybkg.shape[0] > 0:
            mybaseline = interpolate_masked_regions(x, sybkg, bkg, kind=kind, axis=0)
        else:
            mybaseline = sybkg
    else:
        mybaseline = sybkg


    return mybaseline

def interpolate_masked_regions(x, y_gapped, mask, kind, axis):
    '''
    Takes an x array and a partial y array computed over the region mask >=0 and interpolates across masked regions
    :param x: full x values
    :param y_gapped: y values outside the mask region
    :param mask: mask used to generate y
    :param kind: parameter of interpolate.interp1d.  Typically linear
    :param axis: which axis of y has the gaps
    :return: y values for all of x
    '''
    #Interpolate across stimuli using method specified by kind.  Create interpolation function
    x_gapped = x[mask==False]
    x_gaps = x[np.abs(mask)==True]
    stm_interp = interpolate.interp1d(x_gapped, y_gapped, kind=kind, axis=axis)
    # Apply interpolation function to missing time points
    y_interp2 = stm_interp(x_gaps)
    #Merge the gaps and gapped functions
    if axis == 0:
        ax0 = x.shape[0]
        ax1 = y_gapped.shape[1]
        print("shape 0", ax0,ax1)
    else:
        ax0 = y_gapped.shape[0]
        ax1 = x.shape[0]
        print("shape 1", ax0,ax1)

    interpolation_result = np.zeros([ax0,ax1], dtype=np.float32)
    if axis == 0:
        interpolation_result[mask==False,:]= y_gapped
        interpolation_result[np.abs(mask)==True,:] = y_interp2
    else:
        interpolation_result[:, mask==False]= sybkg
        interpolation_result[:, np.abs(mask)==True] = y_interp2
    return interpolation_result

def linear_regression(t, y, mask):
    '''
    Performs multiple linear regression on y matrix
    :param t: one d time points
    :param y: samples array with one dimension matching t
    :param mask: regions to exclude from regression by setting to non-zero
    :return:
    '''
    # Create linear regression object
    regr = linear_model.LinearRegression()
    # Train the model using the training sets
    big_times =  np.ones_like(y) * np.expand_dims(t, axis=1)
    regr.fit(big_times[mask==False,:], y[mask==False,:])
    # Make predictions using the testing set
    y_reg = regr.predict(big_times)
    return y_reg

def parse_runline(arg_list):
    '''
    Takes key=value repeats in arg_list and converts to dictionary or returns None
    :param arg_list: list of values after run DataMan key1=value1 key2=value2.
    :return:
    '''
    load_dict = {}
    for n, arg in enumerate(arg_list):
        if n > 0:   # Ignore 'DataMan'
            load_dict = helpers.line_spitter(arg, load_dict, int_keys_ok=True)
    if len(load_dict) > 0:
        if "is_value" in load_dict.keys() and len(load_dict) == 1:
            return load_dict["is_value"]
        else:
            return load_dict
    else:
        return None

def index_frames(times, start_t, stop_t, min_frame_time="calculate", precision=1.0e-6):
    '''
    Creates indexes for slicing/extracting data based on intevals where events start and stop
    :param times: frame times from imaging
    :param start_t:  odor delivery start time +/- buffer time
    :param stop_t:  odor delivery stop time +/- buffer time
    :param frame_time: maximum time to consider for the current image is start to start + frame_time
    :return:  Indexes to the imaging data that correspond to the start and stop times.  Add +1 to i_end to include i_end in slice
    '''
    precision = pq.Quantity(precision, pq.sec)
    if min_frame_time == "calculate":
        min_frame_time = safe_frame_time(times, precision)
    try:
        i_start = int(np.where(times + min_frame_time >= start_t)[0][0])
        #print("start", i_start)
        i_end = int(np.where(times + min_frame_time >= stop_t)[0][0])
        #print("stop", i_end)
    except:
        print("inframe", start_t, stop_t, len(times), min_frame_time)
        print("Hit end of available times, using end of time series to end odor delivery")
        i_end = len(times)-1
    return i_start, i_end

def safe_frame_time(times, precision = 1.0e-6):
    '''
    Returns a frame time to use in slicing data that is slightly below the shortest inter-frame time in the data times
    :param times: the frame start times
    :param precision: How much below the minimum interframe time to set the frame_time.
    :return: a safe frame_time to use for combining data with previous timepoints.
    '''
    frame_time = np.min(times[1:] - times[:-1]  - precision)
    return frame_time

def make_quantity(n_id, return_unit=None, verbose=False):
    '''
    Takes an input and outputs a Quantity using the unit return_unit if provided.  If already a quantity,
    will strip the current unit and apply the return_unit.
    :param n_id: Input to use to construct quantity.  Can be a Quantity, float, list, array, tuple or string
    :param return_unit: requested return unit if given.  Must make sense.
    :return:  quantity in requested unit if applicable.
    '''
    return_quant = None
    if isinstance(n_id, pq.Quantity):   #Used to strip unit from a Quantity and force to new unit
        if verbose: print("Got a Quantity")
        if return_unit is None:
            return_quant = n_id
        else:
            return_quant = pq.Quantity(np.array(n_id), return_unit)
    elif isinstance(n_id, (float, int)):  #Requires a return unit to convert
        if verbose: print("No unit provided, assume", return_unit)
        try:
            if return_unit:
                return_quant = pq.Quantity(float(n_id), return_unit)
            else:
                return_quant = pq.Quantity(float(n_id), '')   #Return dimensionless quantity
        except:
            return_quant = None

    elif isinstance(n_id, (list, tuple, np.ndarray)):
        if return_unit is None and len(n_id) == 2:    #Handles case of list that has a number and unit
            return_quant = pq.Quantity(float(n_id[0]), n_id[1])
        elif return_unit is not None:              #Handles case of combining array, list , tuple with provided unit
            return_quant = pq.Quantity(n_id, return_unit)
    elif isinstance(n_id, str):   #Handles the case of a string with a number and unit, like "2.0 pA"
        rn = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", n_id)
        numb = float(rn[0])   #rn[0] is number part of string, rest is a unit
        unit = n_id.replace(rn[0], "").strip()   #pops off number and eliminates white space isolating the unit string
        try:
            if unit:   #Unit in string and potentially a rescale is also needed
                return_quant = pq.Quantity(numb, unit)
                if return_unit is not None:
                    return_quant = return_quant.rescale(return_unit)
            elif return_unit is not None:    #Unit was not in string, so use return_unit to make quantity.
                return_quant = pq.Quantity(numb,return_unit)
            else:
                return_quant = pq.Quantity(numb,'')  #Return dimensionless quantity for ID'd number
        except:
            return_quant = pq.Quantity(float(n_id), '')  #Return dimensionless quantity for input string
        if  verbose: print("New Quantity", return_quant)
    else:
        raise ValueError("ERROR:", n_id, "NOT CONVERTABLE")
    return return_quant

def restart_docker(docker_dir=None, docker_bat=None):
    '''
    Can refresh the docker database if it is out of date.
    :param docker_dir: root folder for the Docker image
    :param docker_bat: program to run in Windows 10.
    :return:
    '''
    cwd = os.getcwd()
    if docker_dir is None:
        my_paths = helpers.My_paths()
        docker_dir = my_paths("docker")
    os.chdir(docker_dir)
    if docker_bat is None:
        start_up_bat = docker_dir + "\\start_docker.bat"
    else:
        start_up_bat = docker_dir + "\\" + docker_bat
    os.startfile(start_up_bat)
    # Return to the original directory
    os.chdir(cwd)

def resp_dict(dp, odor_stim, signal_region, sham=False):
    '''
    This creates a response dictionary summarizing the responses of all glomeruli to a given odor
    :param dp:
    :param odor_id:
    :param region:
    :return:
    '''
    lsdata = dp.exp_dict["sdata-baseline"]
    t_data = dp.times

    resp_dict = {}
    for i, (o_slice, po_slice) in enumerate(zip(odor_stim, signal_region)):
        ls_slice = lsdata[slice(po_slice[0],po_slice[1]),:]
        stim_slice = t_data[o_slice[1]-1] - t_data[o_slice[0]]
        # start time, duration, response curve
        resp_dict[i] = {"start":t_data[o_slice[0]], "duration": stim_slice, "response": ls_slice, "integral": np.sum(ls_slice,axis=0)}
        if sham:
            large = np.sum(np.where(resp_dict[i]["integral"] > 2.0, 1, 0))
            if large > 0.5:
                print("large_sham", i, po_slice[0],po_slice[1], np.where(resp_dict[i]["integral"] > 2.0)[0])
    return resp_dict


def matching_noz(odor_info):
    '''
    Feed in a dm.nozzle_to_odor_map dict and will return nozzles that share a common odor
    of the single value of the nozzle with that odor
    :param odor_info:
    :return:
    '''
    odors = odor_info
    matched_nozzles = []
    identical_odors = {}
    unique_nozzles = {}
    for key, value in odors.items():
        for key2, value2 in odors.items():
            unique_key = True
            if key != key2 and value == value2:
                keys = [key, key2]
                matched_nozzles.append(sorted(keys))
                unique_key = False
            if unique_key:
                unique_nozzles[key] = [key]
    for i in range(1, max(odors.keys())+1):
        for match in matched_nozzles:
            #print("match", i, match)
            if i in match:
                if i not in identical_odors.keys():
                    found = False
                    for k, val in identical_odors.items():
                        if i in val:
                            found = True
                            #print("extending", k, i, match)
                            val.extend(match)
                    if not found:
                        identical_odors[i] = match
                        #print("new", i, match)

                elif i in identical_odors.keys():
                    #print("key",i, match)
                    identical_odors[i].extend(match)
    return_odors = {}
    for key, value in unique_nozzles.items():
        return_odors[odor_info[key]] = value
    for key, value in identical_odors.items():
        return_odors[odor_info[key]] = list(set(value))
    return return_odors


def summarize_results(dmf):
    '''
    Examines the data from DataMan objects and extracts out summary info, like actual integrals, ones that were significant,
    number of significant glomeruli, etc.
    :param dmf: Normally a DataManFactory but can take integers, integer lists, or DataMan objects or lists of DataMan Objects
    :return: summary dictionary
    '''
    expt_summary = {}
    id = None
    try:
        id = dmf.what_am_i
    except:
        if isinstance(dmf, (int)):
            dmf = DataManFactory(dmf)
        elif isinstance(dmf, (list, tuple)):
            dmf = DataManFactory(*dmf)
        else:
            print("Requires a DataManFactory, not ", type(dmf))
            return
    if id != "DataManFactory":
        if id == "DataMan":
            dmf = DataManFactory(dmf)
        else:
            print("Please provide a DataManFactory, not ", id)
            return
    for dm1 in dmf:
        i = dm1.expt_id
        if "fluorescence" in dm1.available_data and 'odor_trials' in dm1.available_data and dm1.pd_odor_trials.shape[0] != 0:
            expt_summary[i] = {}
            expt_summary[i]["restriction"] = dm1.select_expt(i)
            summ = dm1.exp_dict["response_cand"]["summary"]
            expt_summary[i]["pd_candidates"] = dm1.exp_dict["pd_candidates"]
            expt_summary[i]["response_cand"] = dm1.exp_dict["response_cand"]
            n_pos = 0
            n_neg = 0
            for summary in summ.values():
                n_pos += summary["n_pos"]
                n_neg += summary["n_neg"]
            expt_summary[i]["n_pos"] = n_pos
            expt_summary[i]["n_neg"] = n_neg
        else:
            expt_summary[i] = "No Fluorescent Odor Responses"
    return expt_summary    
            
class Factory():
    _what_i_hold = None
    '''
    Defines the general type and provide a clean way to query what this is using what_am_i property
    '''
    def __init__(self, *args, **kwargs):
        self._warehouse = {}
        self._what_i_am = "Factory"

        # Not yet implemented. Needs asynchrony to be useful and haven't figured that out yet
        self.build_on_demand = False
        self.active_key = None
        self.unbuilt_args = []
        self.unbuilt_kwargs = {}
        self.add_in(*args, **kwargs)

    @property
    def what_am_i(self):
        return self._what_i_am

    @property
    def stored_ids(self):
        return list(self._warehouse.keys())

    @property
    def holding_what(self):
        return self._what_i_hold

    # provides output to len() function.  Might represent number of objects contained within this object.  Allows enumerate to be used in iterations on the instance
    def __len__(self):
        return len(self._warehouse)

    # Returns instance which is used for the iteration.  Not sure why it is needed, seems redundant but it is.

    def __iter__(self):
        for key, dm in self._warehouse.items():
            self.active_key = key
            yield dm

    def __next__(self):
        for key, dm in self._warehouse.items():
            self.active_key = key
            yield dm

    def __call__(self, *args, **kwargs):
        self.add_in(*args, **kwargs)

    def get_by_id(self, id):
        return_items = []
        if id in self._warehouse.keys():
            return self._warehouse[id]
        else:
            print(id, "not found in warehouse")

    def order_up(self, order_num):
        tmp_pd = self.expts.sort_index()
        return tmp_pd.iloc[order_num,:]

    def add_in(self, *args, **kwargs):
        '''
        Adds items to the warehouse with a unique id
        :param args:
        :param kwargs:
        :return:
        '''
        for arg in args:
            returned_arg = self._add_check(arg)
            if returned_arg:
                self.unbuilt_args.append(returned_arg)
        for key, value in kwargs.items():
            returned_kwarg = self._add_check(value)
            if returned_kwarg:
                self.unbuilt_kwargs[key] = returned_kwarg

    def replace(self, *args, **kwargs):
        '''
        Replaces items in the warehouse with the same key
        :param args:
        :param kwargs:
        :return:
        '''
        for arg in args:
            returned_arg = self._add_check(arg, replace=True)
            if returned_arg:
                self.unbuilt_args.append(returned_arg)
        for key, value in kwargs.items():
            returned_kwarg = self._add_check(value, replace=True)
            if returned_kwarg:
                self.unbuilt_kwargs[key] = returned_kwarg

    def remove(self, *args):
        '''
        Provides a method to remove items from warehouse.  Provides them as a return,
        but return can be ignored if you only want to delete
        :param args: expt_ids for items to remove
        :return:
        '''
        pop_list = []
        for arg in args:
            if arg in self._warehouse.keys():
                pop_list.append(self._warehouse.pop(arg))
        if len(pop_list) < 2:
            return pop_list[0]
        else:
            return pop_list

    def clean_warehouse(self, retrieve=False):
        '''
        Method to remove all items stored in warehouse.
        :param retrieve: If True will return the filled warehouse before cleaning in the object.
        :return: Optional return of what had been in the warehouse
        '''
        if retrieve:
            return_dict = copy.deepcopy(self._warehouse)
            self._warehouse = {}
            return return_dict
        else:
            self._warehouse = {}

    def _add_check(self, item, replace=False):
        '''
        Checks to make sure the item is appropriate for the warehouse and gives it a unique id
        if parameter replace if False.
        :param item: candidate to add to warehouse
        :return: problem entering so False if added correctly or the item if inappropriate
        '''
        try:
            item.what_am_i
        except:
            # Build on demand not implemented yet
            if self.build_on_demand:
                self._warehouse[item] = "on_demand"
                self.active_key = item
                return False
            return item
        if isinstance(item.what_am_i, (str)):
            if self._what_i_hold is None:
                self._what_i_hold = item.what_am_i
            if self._what_i_hold == item.what_am_i:
                key = item.expt_id
                if not replace:
                    while key in self._warehouse.keys():
                        key = float(key) + 0.1
                self._warehouse[key] = item
                self.active_key = key
            else:
                print("Warehouse only hold one type of item, {0}, not {1}".format(self._what_i_hold, item.what_am_i))
                return item
        return False

class DataManFactory(Factory):
    '''
    This class expands upon the Factory class to build and supply DataMan objects to other analyses
    Like other Factory objects, it supports access of the dm objects by iteration, ID, or order,

    '''
    _what_i_hold = "DataMan"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._what_i_am = "DataManFactory"

        location = None   # local Flag for computer being used to correctly set up paths
        self.import_from = None   # Source of the data either dj, csv or pickle
        self.db_upload_available = False
        if len(self._warehouse) > 0:
            dm_return = next(self)
            tmp_dm= next(dm_return)
            self.my_paths = tmp_dm.my_paths
            try:
                self.expts = tmp_dm.expts
            except:
                self.db_connect()
                if self.db_upload_available:
                    self.find_expts()
        else:
            for arg in args:
                if arg in ("home", "work_1", "work_2"):
                    location = arg
            try:
                self.my_paths = helpers.My_paths(location, **kwargs)
            except BADFILENAME:
                print("Bad path.init DataMan construction aborted")
                path_to_init = input("Provide location of paths.init or enter to abort construction")
                if path_to_init != "":
                    self.my_paths = helpers.My_paths(location, **{"paths_init": path_to_init})
                else:
                    raise BADFILENAME
            self.db_connect()
            if self.db_upload_available:
                self.find_expts()
        args_to_build = copy.deepcopy((self.unbuilt_args))
        self.unbuilt_args = []
        kwargs_to_build = copy.deepcopy((self.unbuilt_kwargs))
        self.unbuilt_kwargs = {}
        self.build_it(*args_to_build, **kwargs_to_build)

    def __call__(self, *args, **kwargs):
        args = list(args)
        args.append("return_build")
        return_build = self.build_it(*args, **kwargs)
        if len(return_build) == 1:
            return_build = return_build[0]
        return return_build

    def find_expts(self):
        '''
        Pulls out a pandas dataframe of available experiments by animal, session, recording id's, etc
        from the DataJoint Database
        :return: loads info into pandas DataFrame object
        '''
        self.expts = pd.DataFrame((self.odor.MesoMatch).fetch())


    def db_connect(self, **kwargs):
        '''
        Sets up connection to DataJoint and returns status
        :param kwargs: Provide alternative config info {"host": "strhost", "user": struser, "pwd": strpwd}
        :return: status of connection put into self.db_upload_available
        '''
        host = self.my_paths.datajoint_db
        user_name = self.my_paths.username
        password = self.my_paths.password

        if kwargs:
            for key, value in kwargs.items():
                if key == "host":
                    host = value
                elif key == "user":
                    user_name = value
                elif key == "pwd":
                    password = value
        try:
            dj.config['database.host']= host
            dj.config['database.user'] = user_name
            dj.config['database.password'] = password
            dj.conn()
        except:
            print("Problem connecting to DataJoint ", sys.exc_info()[0])
            self.db_upload_available = False
        else:
            self.db_upload_available = True

        self.db_upload_available = self.test_db_sockets()

    def test_db_sockets(self):
        '''
        This method tests the connection to the datajoint database
        and gets the list of available experiments if the connection is good
        :return: Connection status for db
        '''
        if self.db_upload_available:
            # Attempt to make a connection with odor part of the DataJoint db
            try:
                self.odor = dj.create_virtual_module('', 'pipeline_odor')
            except:
                print("DB_test: Problem with DataJoint Connection.  Recheck to see if VPN or network has a problem.")
                return False

            self.find_expts()
            print("Datajoint connection passes test")
            return True

        else:
            print("DataJoint connection information not loaded. First run self.db_connect()")
            return False


    def build_it(self, *args, **kwargs):
        '''
        Builds DataMan objects in the Factory and stores in warehouse.  Initial build code will
        initialize the warehouse.  Later additions will require args or kwargs
        :param args:  Expt_ids to add to warehouse
        :param kwargs: instructions for adding more dm objects by DataJoint
        :return:
        '''
        return_dms = []
        return_build = False
        for arg in args:
            if isinstance(arg, (int)):
                if arg in self._warehouse.keys():
                    self.active_key = arg
                    return_dms.append(self._warehouse[arg])
                else:
                    self.add_in(DataMan(arg, "dj"))
                    return_dms.append(self._warehouse[self.active_key])
            elif isinstance(arg, (str)) and arg in ("all", "everything"):
                for id in self.expts.index():
                    self.add_in(DataMan(id, "dj"))
                    return_dms.append(self._warehouse[self.active_key])
            elif isinstance(arg, (str)) and arg == "return_build":
                return_build = True
            else:
                self.add_in(arg)
                return_dms.append(self._warehouse[self.active_key])

        for key, values in kwargs.items():
            if key == "return_build":
                return_build = values
            if isinstance(values, (str, int, float)):
                self.add_in(DataMan(values, "dj"))
                return_dms.append(self._warehouse[self.active_key])
            elif isinstance(values, (list, tuple)):
                for value in values:
                    if isinstance(value, (int)):
                        self.add_in(DataMan(value, "dj"))
                        return_dms.append(self._warehouse[self.active_key])
                    else:
                        self.add_in(value)
                        return_dms.append(self._warehouse[self.active_key])
        if return_build:
            return return_dms



class TestFactory(Factory):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.build_it()

    def build_it(self):
        for arg in self.build_args:
            if isinstance(arg, (int)):
                self._warehouse[arg] = TestItem()
        for key, values in self.build_kwargs.items():
            if key in ("test"):
                for value in values:
                    self._warehouse[value] = TestItem()

class TestItem():
    def __init__(self):
        self._what_i_am = "TestItem"
        self._exp_id = float(np.random.randint(0,100))

    @property
    def what_am_i(self):
        return self._what_i_am

    @property
    def exp_id(self):
        return self._exp_id

class PlaceHolder():
    _what_i_am = "Placeholder"
    def __init__(self, for_what=None):
        self._for_what = for_what
    @property
    def what_am_i(self):
        return self._what_i_am

    @property
    def holding(self):
        return self._for_what

def process_respiration(t, y, y_red=11, delta=20, sensitivity=0.01, dt = 1./5000.):

    if isinstance(y, (pq.Quantity)) and not isinstance(sensitivity, (pq.Quantity)):
        unit = y.units
        sens = make_quantity(sensitivity, unit)
    else:
        sens = sensitivity
    bp_y = band_pass(y, corner_f=(0.5,10.), dt=dt)
    if y_red:
        resp_t, resp_y = data_reduce(t,bp_y, y_red)
    else:
        resp_t = t
        resp_y = bp_y
    breath = np.zeros_like(resp_y)

    for i in range(delta,len(resp_y)-delta):
        slope = resp_y[i+delta] - resp_y[i-delta]
        if slope > sens:
            breath[i] = 1
        elif slope < -sens:
            breath[i] = -1
        else:
            breath[i] = 0

    return resp_t, resp_y, breath

def process_respold(t, y, f_corner=[0.5, 20.], by=11, delta=100, sensitivity=0.01, dt = 1./5000., interp=True):
    interp_y = None
    interp_breath = None
    if isinstance(y, (pq.Quantity)) and not isinstance(sensitivity, (pq.Quantity)):
        unit = y.units
        sens = make_quantity(sensitivity, unit)
    else:
        sens = sensitivity
    breath = np.zeros_like(y)
    resp_y = band_pass(y, corner_f=f_corner, dt=dt)
    resp_t, resp_y = data_reduce(t, resp_y, by)
    for i in range(delta,len(resp_y)-delta):
        slope = resp_y[i+delta] - resp_y[i-delta]
        if slope > sens:
            breath[i] = 1
        elif slope < -sens:
            breath[i] = -1
        else:
            breath[i] = 0
    if interp:
        interp_y_func = interpolate.interp1d(t, resp_y, kind="linear",bounds_error=False)
        interp_resp_func = interpolate.interp1d(t, breath, kind="linear", bounds_error=False)
        interp_y = interp_y_func(t_interp)
        interp_breath = interp_resp_func(t_interp)
    return_dict = {"t": resp_t, "y": resp_y, "breaths": breath, "interp_resp": interp_y, "interp_breath": interp_breath}
    return return_dict

def testisolate_responses(self, odor_id, leak_method, region="pre_maxodor_post"):
    '''
    Isolates leak subtracted slices of data in response to an odor
    :param odor_id: odor of interest
    :param sdata: smoothed data
    :param ls_data: globally leak subtracted data
    :param d_times: delta times for the experiment
    :param slices: slices dictionary
    :return: sliced leak subtracted signals and delta times for each slice
    '''
    # Select region to generate leak based on leak method
    olap_signal = []
    olap_stim_time = []
    olap_true_stim_time = []
    biggest_slice = 0
    if region == "common":
        stim_mask, stim_region = self.mask_slice_manager(how_to="common", by_odor=odor_id)
    else:
        stim_mask, stim_region = self.mask_slice_manager(how_to="odor", by_odor=odor_id)

    leaker = np.ones_like(self.exp_dict["sdata"], dtype=np.float32)
    # buffer_pts method uses average of values in pre-buffer region to approximate baseline
    __, leak_slices = self.mask_slice_manager(how_to="pre_stim", by_odor=odor_id)
    signal_mask, signal_region = self.mask_slice_manager(how_to=region, by_odor=odor_id)
    first_slice = True
    concat_true_time = self.times[np.abs(signal_mask) > 0]
    concat_true_stim_time = self.times[np.abs(stim_mask) > 0]
    concat_resp_time = self.exp_dict["avg_dt"] * np.linspace(0, concat_true_time.shape[0], concat_true_time.shape[0], endpoint=False)

    cum_time = 0
    for a_slice, po_slice, stim_slice in zip(leak_slices, signal_region, stim_region):
        if not self.silence: print("pa", a_slice, po_slice)
        avg_leak = np.average(self.exp_dict["sdata"][slice(a_slice[0], a_slice[1]),:], axis=0)
        leak = leaker * avg_leak
        if leak_method == 'buffer_pts':
            ls_slice = self.exp_dict["sdata"][slice(po_slice[0],po_slice[1]),:] - leak[slice(po_slice[0],po_slice[1]),:]
        elif leak_method == "no_leak":
            ls_slice = self.exp_dict["sdata"][slice(po_slice[0],po_slice[1]),:]
        else:
            ls_slice = self.exp_dict["sdata-leak"][slice(po_slice[0],po_slice[1]),:]
        olap_true_stim_time.append(self.times[slice(po_slice[0],po_slice[1])])
        stim_start = stim_slice[0] - po_slice[0] + cum_time
        stim_end = stim_slice[1] - po_slice[0] + cum_time
        cum_time += (po_slice[1] - po_slice[0])
        t_slice = concat_resp_time[slice(stim_start, stim_end)]
        tr_slice = concat_true_time[slice(stim_start, stim_end)]

        if first_slice:
            response_data = ls_slice
            first_slice = False
            response_stim_time = t_slice
        else:
            response_data = np.append(response_data, ls_slice, axis=0)
            response_stim_time = np.append(response_stim_time, t_slice, axis=0)
        biggest_slice = max(biggest_slice, ls_slice.shape[0])
        olap_signal.append(ls_slice)
        olap_stim_time.append(t_slice)

    overlap_resp_time = self.exp_dict["avg_dt"] * np.linspace(0, biggest_slice, biggest_slice, endpoint=False)

    concat = {"resp": response_data, "time": concat_resp_time, "stim_time": response_stim_time, "true_stim_time": concat_true_stim_time}
    overlap = {"resp": olap_signal, "time": overlap_resp_time, "stim_time": olap_stim_time, "true_stim_time": olap_true_stim_time}
    return concat, overlap

def find_modes(num_list, sensitivity=0.25, sig_figs=1):
    modes = []
    for num in num_list:
        if len(modes) == 0:
            modes.append(num)
        else:
            diff = [abs(n - num) for n in modes]
            min_diff = min(diff)
            if min_diff > sensitivity:
                modes.append(num)
    modes = sorted(modes)
    modes = [n_decimals(num, sig_figs) for num in modes]
    return sorted(modes)

def n_sig_figs(x, n=1):
    '''
    Performs sig_fig rounding to n places
    :param x: number to use
    :param n: number of sig figs
    :return: x truncated to n sig figs
    '''
    rounded = lambda x, n: x if x == 0 else round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))
    return rounded(x, n)

def n_decimals(x, n=1):
    '''
    Truncates float numbers to have n decimal.
    :param x: float to truncate
    :param n:
    :return:
    '''
    # Numbers near 0 are treated as 0 which has no sign

    sign = np.sign(np.array([x]))[0]
    int_x = int(abs(x) * (10.**n) + 0.5)
    return sign * float(int_x) / 10.**n

def hr_min_sec(time_in_sec):
    t = float(time_in_sec)
    min, sec = divmod(t, 60)
    hr, min = divmod(min, 60)
    #subsec = sec - int(sec)
    return_str = "{0:.0f}:{1:.0f}:{2:.3f}".format(hr, min, sec)
    return return_str

def pickle_dump(*args):
    try:
        id = args[0].who_am_i
        if id == "DataMan":
            dmf = DataManFactory(*args)
        elif id == "DataManFactory":
            dmf = args[0]
    except:
        dmf = DataManFactory(*args)
    for dm in dmf:
        dm.save_as_pickle()



if __name__ == "__main__":
    np.set_printoptions(precision=3)
    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 1000)
