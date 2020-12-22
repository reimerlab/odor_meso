import quantities as pq
from matplotlib import pyplot as plt
# from matplotlib.backend_bases import MouseEvent
from matplotlib import collections  as mc
import numpy as np
import logging
from filters import smooth, band_pass   # A local module
# Suppresses some annoying matplotlib messages
logging.getLogger().setLevel(logging.CRITICAL)
from DataMan import DataMan, DataHandler, make_quantity, DataManFactory
import math
import copy
from matplotlib.widgets import CheckButtons
from PyQt5.QtWidgets import QMainWindow, QFileDialog, QApplication, QTableView
from PyQt5 import QtCore
app = QApplication.instance()
if app is None:
    # if it does not exist then a QApplication is created
    app = QApplication(sys.argv)


class plot_dm():
    '''
    plots:
    1) Integrals                    plot_integrals()
    2) Concat responses for a glom  pg(glom)
    3) Plot All                     plot_all(glom)
    4)
    5) plot overlapping responses  plot_overlap
    6) self.f0_fig_id = 6
    7)  self.df_f_fig_id = 7
    '''
    _what_i_am = "Plot_dm"
    n_figures = 10
    def __init__(self, *args, **kwargs):
        '''
        Typical args and kwargs are to shift plot indexes so multiple plots of the same type can exist together or to provide
        parameters to build a DataMan or DataManFactory object
        :param args:
        :param kwargs:
        '''
        self.data_loaded = False
        # Don't close preexisting figs unless instructed to by uinit.init override from DataMan
        self.close_preexisting_figs = False
        # Will create a DataMan object if not provided
        #id = False
        self.dm = None
        self.dm_factory = None
        id = 12   # Default experiment to open
        for arg in args:
            if not isinstance(arg, (int, str, list, tuple, float, dict)):
                try:
                    id = arg.what_am_i
                except:
                    continue
            print("id is ", id)
            if id == "DataMan":
                if self.dm is None:
                    self.dm = arg
                if self.dm_factory is None:
                    self.dm_factory = DataManFactory(self.dm, arg)
                else:
                    self.dm_factory.add_in(arg)
            elif id in ("Factory", "DataManFactory", "DataManConstructor", "dmf", "dmc"):
                self.dm_factory = arg
                self.dm_gen = next(self.dm_factory)
                self.dm = next(self.dm_gen)

        if self.dm is None:
            for key, value in kwargs.items():
                if key in ("DataHandler", "DataMan", "DataManager", "dm"):
                    self.dm = value
                elif key in ("Factory", "DataManFactory", "DataManConstructor", "dmf", "dmc"):
                    self.dm_factory = value
                    self.dm_gen = next(self.dm_factory)
                    self.dm = next(self.dm_gen)
        if self.dm is None:  # if all else fails...
            self.dm = DataMan(*args, **kwargs)

        self.lines = {}
        def register_polar():
            fig = plt.figure(99)
            ax = plt.subplot(111, projection='polar')
            ispolar = ax.__class__
            plt.clf()
            plt.close(99)
            return ispolar
        # Container needed if multiple plot_dm objects are plotting at the same time to ensure
        # Unique figure #'s
        self.fig_shift = 0
        # TODO: See if there is a more elegant way to do this.
        self.ispolar = register_polar()
        try:
            self.extract_from_dm(self.dm)
            self.data_loaded = True
        except:
            print("Unable to extract from DataMan, probably need to load with data")
        #Not from init.  These are housekeeping parameters
        self.active_plots = [0] * self.n_figures
        self.polar1 = "inactive"   # Keeps track of whether there is a current polar plot in place
        # Index of the glomerulus that is currently active
        self.plot_id = 0
        self.polar_max = 0.0
        self.active_odors = {}
        self.active_odor_set = 0
        self.x_range_all = [0.0, 400.0]
        self.add_vbars = True
        self.vbar_trans = 0.4
        # Boolean to monitor whether plots are cleared for pen choice
        self.clear_plot = [True] * self.n_figures
        self.manual_baseline = False  # Not yet implemented manual way to set baseline
        if self.close_preexisting_figs:
            plt.close("all")
        # fig containers for main plots
        self.figs = [0] * self.n_figures
        self.active_figs = [False] * self.n_figures
        if self.data_loaded:
            self.color = [self.default_color] * self.n_figures   # Current pen color
            self.active_odors = self.select_active_odors(**kwargs)
        else:
            self.color = [0] * self.n_figures   # Current pen color

    @property
    def what_am_i(self):
        return self._what_i_am

    @property
    def odor_set(self):
        return self.active_odor_set


    def next_dm(self):
        try:
            id = self.dm_factory.what_am_i
        except:
            print("No DataManFactory loaded")
            return
        self.dm = next(self.dm_gen)
        print("Loading experiment ", self.dm.exp_id, self.dm.exp_dict["restriction"])

    def new_set(self, n=None):
        '''
        Which odor set to plot
        :return:
        '''
        if n is None:
            self.active_odor_set += 1
            if len(self.active_odors.keys()) <= self.active_odor_set:
                self.active_odor_set = 0
        elif n >= len(self.active_odors.keys()):
            print("Maximum n_set is", len(self.active_odors.keys()))
        else:
            self.active_odor_set = n

    def reload(self, dm=None):
        if dm is not None:
            self.dm = dm
            try:
                self.extract_from_dm(self.dm)
                self.data_loaded = True
            except:
                print("Unable to extract from DataMan, probably need to load with data")
                return
        self.color = [self.default_color] * self.n_figures   # Current pen color
        self.active_odors = self.select_active_odors()

    def select_active_odors(self, **kwargs):
        '''
        If odor_labels is greater than can be fit on one figure, then split into multiple groups
        :param kwargs: a defined splitting of the odor labels key is set #, value={subset odor labels dict}
        :return: 
        '''
        active_odors = {}
        preferred_grid_size = self.plot_grid[0] * self.plot_grid[1]
        n_odor_variants = self.dm.n_unique_odor_combos
        if kwargs:
            self.active_odors = kwargs
            return
        
        active_labels = copy.deepcopy(self.dm.odor_labels)

        if n_odor_variants > preferred_grid_size:
            common_odor = active_labels.pop(1)
            n_to_pick = preferred_grid_size-1
            n_picked = 0
            i_set = 0
            active_odors[i_set] = {}
            active_odors[i_set][1] = common_odor
            for key, odor in active_labels.items():
                active_odors[i_set][key] = odor
                n_picked += 1
                if n_picked == n_to_pick:
                    i_set += 1
                    active_odors[i_set] = {}
                    active_odors[i_set][1] = common_odor
                    n_picked = 0
        else:
            active_odors[0] = active_labels
        return active_odors

    def extract_from_dm(self, dm):
        # Close matplotlib figures that may still be open from earlier interactive sessions
        self.close_preexisting_figs = dm.close_preexisting_figs
        # Pen colors to use for plots.  Automatically steps colors if plotting with overlap. Slice symbol required.
        self.colors_id = dm.colors_id
        # Default or reset pen color
        self.default_color = dm.default_color
        self.title_size = dm.title_size
        self.y_range = dm.y_range
        # Default plot shape with 3 rows and 4 columns of plots
        # need space between [ ] and values to slice correctly
        self.plot_grid = dm.plot_grid
        # if Single clear old plots.  If overlap, plots on same axis with pen color iterated
        self.plot_mode = dm.plot_mode
        # Size of marker for odor delivery
        self.marker  = dm.marker
        # Whether to plot a polar summary in upper left
        self.polar_summary = dm.polar_summary
        # Whether to plot response data with baseline subtraction
        self.plt_baselinesubtracted = dm.plt_baselinesubtracted
        #plot_zero_line adds a dashed y=0 line to certain plots
        self.plot_zero_line = dm.plot_zero_line
        self.y_range_integral = dm.y_range_integral
        # Whether to include a sorted plot of the integral responses to see distribution
        self.sort_integrals = dm.sort_integrals
        # Criterion for ID of positive integrals
        self.sig_prob = dm.sig_prob
        # Color to use for sorted plotting of integrals in integrals plot
        self.sorted_color = dm.sorted_color
        self.stim_line_width = dm.stim_line_width
        # position is ymin times stim_position if negative or ymax * stim_position if positive
        self.stim_position = dm.stim_position
        # If pickling DataPlot what to do with open figs, options:
        #     pickle: open plots preserved by pickling.  unpickle plots names in old_save_dict["old_figs"]
        #     replot: keeps list of what was in active_plots.  Redoes plots for figs 1, 2, 3
        #     both: option for both pickle and replot.  Reactivation of figs default is pickle unless overrridden
        #     anything else: don't save plots unless overrridden
        self.save_figs = dm.save_figs
        self.time_unit = dm.time_unit
        self.concat_slice = dm.response_mask

        needed_colors = self.dm.n_unique_odor_combos - len(self.dm.colors_id)
        if needed_colors > 0:
            self.add_random_colors(needed_colors)

    def pg(self, g=None):
        '''
        Figs[2] Function to Plot by giving a glomerulus ID. Useful if you want to view a specific glomerulus
        :param g: glomerulus ID
        :return: plot
        '''
        if g  is not None and g in self.dm.glomeruli:
            # index is typically g - 1, due to arrays starting at 0
            self.plot_id = self.dm.glomeruli.index(g)
        elif g  is not None and g not in self.dm.glomeruli:
            print("Requested glomerulus not in data set: ", g)
            return
        else:
            # prevous plotting indexed self.plot_id to +=1 so just use this if None
            pass


        fig_id = 2
        # Check to see if Fig[3] already exists or create it
        if isinstance(self.figs[fig_id], (plt.Figure)):
            self.figs[fig_id]
        else:
            self.figs[fig_id] = plt.figure(fig_id + self.fig_shift)
            self.figs[fig_id].canvas.mpl_connect('close_event', self.handle_close)
            self.active_figs[fig_id] = True

        if self.plot_mode == "single":
            plt.clf()
            self.color[2] = self.default_color
            self.active_plots[2] = [self.dm.glomeruli[self.plot_id]]
        if self.plot_mode == "overlap":
            if isinstance(self.active_plots[2], (list)):
                self.active_plots[2].append(self.dm.glomeruli[self.plot_id])
            else:
                self.active_plots[2] = [(self.dm.glomeruli[self.plot_id])]

        print("plotting glomerulus: ", self.dm.glomeruli[self.plot_id])
        self.plot_resp(self.dm.glomeruli[self.plot_id])
        if self.plot_mode == "overlap":
            self.color[2] += 1
            if self.color[2] >= len(self.colors_id):
                self.color[2] = 0
        self.plot_id += 1
        self.clear_plot[2] = False


    def pn(self, n=None):
        '''
        Figs[2] Function for Plotting using the array index for the data.  Useful particularly with 0 to reset to 1st glomerulus
        or if you forget numbers of a subset of glomeruli
        :param n: index for glomerulus to plot.
        :return: plot
        '''
        if n is not None and n < len(self.dm.glomeruli):
            self.plot_id = n
        elif n is not None and n >= len(self.dm.glomeruli):
            print("Requested glomerulus n not in data set: ", n)
            return
        else:    # prevous plotting indexed self.plot_id to +=1 so just use this if None
            pass

        fig_id = 2
        # Check to see if Fig[3] already exists or create it
        if isinstance(self.figs[fig_id], (plt.Figure)):
            self.figs[fig_id]
        else:
            self.figs[fig_id] = plt.figure(fig_id + self.fig_shift)
            self.figs[fig_id].canvas.mpl_connect('close_event', self.handle_close)
            self.active_figs[fig_id] = True

        if self.plot_mode == "single":
            plt.clf()
            self.color[2] = self.default_color
            self.active_plots[2] = [self.dm.glomeruli[self.plot_id]]
        if self.plot_mode == "overlap":
            if isinstance(self.active_plots[2], (list)):
                self.active_plots[2].append(self.dm.glomeruli[self.plot_id])
            else:
                self.active_plots[2] = [(self.dm.glomeruli[self.plot_id])]

        print("plotting glomerulus: ", self.dm.glomeruli[self.plot_id])
        self.plot_resp(self.dm.glomeruli[self.plot_id])
        if self.plot_mode == "overlap":
            self.color[2] += 1
            if self.color[2] >= len(self.colors_id):
                self.color[2] = 0

        self.clear_plot[2] = False
        self.plot_id += 1


    def plot_all(self, glom=None, *args, **kwargs):
        '''
        Plots all data for a glomerulus in Figs[3] with stimuli plotted beneath
        :param glom: glomerulus to use
        :return:
        '''
        if "fig_id" in kwargs.keys():
            fig_id = kwargs["fig_id"]
        else:
            fig_id = 3
        # Check to see if Fig[3] already exists or create it
        if isinstance(self.figs[fig_id], (plt.Figure)):
            self.figs[fig_id]
        else:
            self.figs[fig_id] = plt.figure(fig_id + self.fig_shift)
            self.figs[fig_id].canvas.mpl_connect('close_event', self.handle_close)
            self.active_figs[fig_id] = True

        # checks to see if there are plots in Fig[3] or if not create a new axis for plotting
        if len(self.figs[fig_id].axes) > 0:
            ax = self.figs[fig_id].axes[0]
            plt.sca(ax)
        else:
            if args:
                ax = plt.subplot(args[0])
            elif "sub_fig" in kwargs.keys():
                ax = plt.subplot(kwargs["sub_fig"])
            else:
                ax = plt.subplot(111)
            plt.sca(ax)
        #print("AP2", self.active_plots)

        # if no glomerulus info is given create a default to plot
        if glom is None:
            if isinstance(self.active_plots[fig_id], (list)) and len(self.active_plots[fig_id]) > 0:
                glom = self.active_plots[fig_id][-1]
            else:
                glom = self.dm.glomeruli[0]

        # Check whether plotting overlap or single plots and handle accordingly
        if self.plot_mode == "single":
            ax.clear()
            self.color[fig_id] = self.default_color
        if len(ax.get_lines()) < 1:
            self.clear_plot[fig_id] = True
            self.color[fig_id] = self.default_color
        else:
            self.clear_plot[fig_id] = False
            self.color[fig_id] += 1

        # Gather the data needed to make the plot and set the y-scaling limits appropriately
        glom_id = self.dm.glomeruli.index(glom)
        # print(odor[0], odor[1], glom)

        # this is a time vs. response plot, so x = times
        xplt = self.dm.exp_dict["times"]

        # Typically plotting baseline_subtracted data unless plt_baselinesubtracted is False
        if self.plt_baselinesubtracted:
            yplt = self.dm.exp_dict["sdata-baseline"][:, glom_id]
        else:
            yplt = self.dm.exp_dict["sdata"][:, glom_id]

        # extract out the y_unit so all y values plot with the same unit skip units if not a pq.Quantity False
        if isinstance(yplt, (pq.Quantity)):
            use_units = True
            y_unit = yplt.units
        else:
            use_units = False
        if isinstance(self.y_range,(list, tuple)):
            if use_units:
                y_range = make_quantity(self.y_range, y_unit)
            else:
                y_range = self.y_range
            ax.set_ylim(y_range[0],y_range[1])  # Graph scaling to use
        else:
            ax.autoscale(enable=True, axis='y')

        ax.autoscale(enable=True, axis='x')

        # Autoscale x.  The plot can later be changed interactively
        if self.x_range_all == "auto":
            ax.autoscale(enable=True, axis='x')
        else:
            if use_units:
                ax.set_xlim(pq.Quantity(self.x_range_all[0], self.time_unit), pq.Quantity(self.x_range_all[1], self.time_unit))
            else:
                ax.set_xlim(self.x_range_all[0], self.x_range_all[1])

        if self.plot_zero_line and self.clear_plot[3]:
            if use_units:
                y_zero = make_quantity(np.zeros_like(xplt), y_unit)
            else:
                y_zero = np.zeros(xplt.shape, dtype=np.float32)
            ax.plot(xplt, y_zero, color="gray", linestyle="--")

        ax.plot(xplt, yplt, color=self.colors_id[self.color[3]], linewidth=3.0)
        if self.stim_position < 0:
            if isinstance(self.y_range,(list, tuple)):
                stim_position = self.y_range[0] * abs(self.stim_position)
            else:
                stim_position = np.min(yplt) * abs(self.stim_position)
        else:
            if isinstance(self.y_range,(list, tuple)):
                stim_position = self.y_range[1] * abs(self.stim_position)
            else:
                stim_position = np.max(yplt) * abs(self.stim_position)

        if use_units and not isinstance(stim_position, (pq.Quantity)):
            sp = make_quantity(stim_position, y_unit)
        else:
            sp = stim_position

        sw = self.stim_line_width

        # plots all the odor stimuli as distinct color line segments
        if self.clear_plot[3]:
            for n, odor_id in sorted(self.dm.odor_labels.items()):
                active_trials = self.dm.pd_odor[self.dm.pd_odor["odorant"] == odor_id]
                start = active_trials["trial_start_time"]
                stop = active_trials["trial_end_time"]
                plt_pts = [[(make_quantity(x1, self.time_unit),sp), (make_quantity(x2, self.time_unit),sp)] for x1, x2 in zip(start, stop)]
                lc = mc.LineCollection(plt_pts, color=self.colors_id[n-1], linewidths=sw)
                ax.add_collection(lc)
            if self.add_vbars:
                for n, odor_id in sorted(self.dm.odor_labels.items()):
                    active_trials = self.dm.pd_odor[self.dm.pd_odor["odorant"] == odor_id]
                    start = active_trials["trial_start_time"]
                    stop = active_trials["trial_end_time"]
                    plt_pts = [(make_quantity(x1, self.time_unit), make_quantity(x2, self.time_unit)) for x1, x2 in zip(start, stop)]
                    for begin, end in plt_pts:
                        ax.axvspan(begin, end, facecolor= self.colors_id[n-1], alpha=self.vbar_trans)

        # Plotting as single plots or overlapping
        #print("AP41", self.active_plots)

        if self.plot_mode == "single":
            self.color[fig_id] = self.default_color
            self.active_plots[fig_id] = [self.dm.glomeruli[glom_id]]
        if self.plot_mode == "overlap":
            if isinstance(self.active_plots[3], (list)):
                self.active_plots[fig_id].append(self.dm.glomeruli[glom_id])
            else:
                self.active_plots[fig_id] = [(self.dm.glomeruli[glom_id])]

        print("plotting glomerulus: ", self.dm.glomeruli[glom_id], self.active_plots)
        if self.plot_mode == "overlap":
            self.color[fig_id] += 1
            if self.color[fig_id] >= len(self.colors_id):
                self.color[fig_id] = 0
        title = "Glomerulus" + str(self.active_plots[fig_id])
        # Adding legend for all lines except dashed 0 line
        lines = [lines for i, lines in enumerate(ax.lines) if i != 0]
        ax.legend(lines, self.active_plots[fig_id], fontsize=self.title_size,loc=1)

        plt.suptitle(title, fontsize=self.title_size)

    def handle_close(self, event):
        '''
        Cleans up a closed figure in matplotlib
        :param event:
        :return: Removes figure from figures and active_figures and resets plotting defaults
        '''
        to_close = None
        event.canvas.figure.axes[0].has_been_closed = True
        print('Closed Figure!')
        self.before_cleanup = [type(a) for a in self.figs]
        for i, fig in enumerate(self.figs):
            if not isinstance(fig, (int)):
                try:      # Find closed figure. Only this figure will have property .has_been_closed
                    if fig.axes[0].has_been_closed:
                        to_close = i
                        break
                except:
                    pass

        print("closing fig ", to_close)
        if to_close is not None:
            plt.close(self.figs[to_close])
            self.figs[to_close] = 0
            self.active_figs[to_close] = False
            self.color[i] = self.default_color
            self.active_plots[i] = 0

    def close_plots(self, close_info=False):
        '''
         used to get rid of plots cleanly.  Any plots that don't close are not owned by this object and can
         be closed using plt.close()
        :param close_info:  specific plots to clear or a list, "all", etc.
        :return:
        '''
        if close_info == False:
            print("No existing plots cleared")
            return
        elif close_info == True or isinstance(close_info, (str)) and close_info.lower() in ("all", "yes"):
            # shuts down any previous fig versions
            close_list = self.active_figs[:]
        elif isinstance(close_info, (list, tuple)):
            close_list = [False] * self.n_figures
            for fig in close_info:
                if self.active_figs[fig]:
                    close_list[fig] = True
        elif isinstance(close_info, (int)):
            close_list = [False] * self.n_figures
            if self.active_figs[close_info]:
                    close_list[close_info] = True
        elif isinstance(close_info, (str)) and close_info.isnumeric():
            close_list = [False] * self.n_figures
            if self.active_figs[int(close_info)]:
                    close_list[int(close_info)] = True
        else:
            print("unknown clear plot command", close_info)
            close_list = []

        # if value in close_list is True, then close this plot and clear it from the lists
        for to_close, truth in enumerate(close_list):
            if truth:
                plt.close(self.figs[to_close])
                self.figs[to_close] = 0
                self.active_figs[to_close] = False

    def plot_resp(self, glom):
        '''
        Plots block of odor responses to specific odors.  Main plotting subroutine called by various
        operations of the object to plot odor responses
        for a particular glomerulus or to overlap responses if flag self.plot_mode = "overlap"
        :param glom:  Which glomerulus ID to use in the plot.  (Note: not the index, the actual ID)
        :return: a matplotlib plot
        '''
        plot_slice = self.concat_slice
        ax = [0]*self.dm.n_unique_odor_combos
        if isinstance(self.figs[2], (plt.Figure)):
            self.figs[2]
        else:
            self.figs[2] = plt.figure(2+ self.fig_shift)
            self.figs[2]._label = "2"
            self.figs[2].canvas.mpl_connect('close_event', self.handle_close)
            self.active_figs[2] = True
        glom_id = self.dm.glomeruli.index(glom)
        if plot_slice != self.dm.exp_dict["response_data"]["region"]:
            self.dm.make_response_data(region=plot_slice, force=True)

        title = "Glomerulus" + str(self.active_plots[2])
        plt.suptitle(title, fontsize=self.title_size)
        plt_dict = self.active_odors[self.active_odor_set]

        for n, n_odor in enumerate(plt_dict.values()):
            n = n+1
            concat_data = self.dm.exp_dict["response_data"][n_odor]["concat"]
            xplt = concat_data["time"]
            yplt = concat_data["resp"]
            yplt_min = yplt.min()
            if "auto" in self.y_range:
                if self.stim_position < 0.8 * yplt_min:
                    self.stim_position = 0.8 * yplt_min
            xstim = concat_data["stim_time"]
            # Ignore upper left plot if puttng a polar summary in that spot
            if self.polar_summary and n==1:
                # TODO: Correct axis creation only when need
                if len(self.figs[2].axes) > n-1:
                    if self.polar_summary and not self.figs[2].axes[0].__class__ == self.ispolar:
                        ax.remove()
                        ax = plt.subplot(self.plot_grid[0], self.plot_grid[1], n, projection='polar')
                        self.polar_response(ax, glom)
                        self.polar1 = "active"
                    else:
                        ax = self.figs[2].axes[n-1]
                        self.polar_response(ax, glom)
                        self.polar1 = "active"
                else:
                    ax = plt.subplot(self.plot_grid[0], self.plot_grid[1], n, projection='polar')
                    self.polar1 = "active"
                    self.polar_response(ax, glom)
            else:
                if len(self.figs[2].axes) > n-1:
                    ax = self.figs[2].axes[n-1]
                    plt.sca(ax)
                else:
                    ax = plt.subplot(self.plot_grid[0], self.plot_grid[1], n)
                    plt.sca(ax)
                    #self.fig2_axes[n] = ax[n]
                self._plot_response(glom_id, xplt, yplt, xstim, ax, n_odor)
        print("Lines = ",len(ax.lines))
        # If desired plot polar summary graph in position 1

    def _plot_response(self, glom_id, xplt, yplt, xstim, ax, n_odor):
        if isinstance(yplt, (pq.Quantity)):
            use_units = True
            unit = yplt.units
        else:
            use_units = False

        #self.fig2_axes[n] = ax[n]
        if self.plot_mode == "single":
            ax.clear()
        if len(ax.get_lines()) < 1:
            empty_plot = True
        else:
            empty_plot = False
        #print("testing axes", ax)
        if isinstance(self.y_range,(list, tuple)):
            if use_units:
                y_range = make_quantity(self.y_range, unit)
            else:
                y_range = self.y_range
            ax.set_ylim(y_range[0],y_range[1])  # Graph scaling to use
        elif self.y_range == "common_autoscale":
            if use_units:
                y_range = make_quantity(self.dm.autoscale_responses, unit)
            else:
                y_range = self.dm.autoscale_responses
            ax.set_ylim(y_range[0], y_range[1])
        else:
            ax.autoscale(enable=True, axis='y')
        ax.autoscale(enable=True, axis='x')
        # print(odor[0], odor[1], glom)
        plt.title = n_odor
        if use_units:
            ystim = make_quantity(np.zeros_like(xstim)+ self.stim_position, unit)

        else:
            ystim = np.zeros(xstim.shape, dtype=np.float32) + self.stim_position

        ax.plot(xplt, yplt[:, glom_id], color=self.colors_id[self.color[2]])
        if self.plot_zero_line and empty_plot:
            if use_units:
                y_zero = make_quantity(np.zeros_like(xplt), unit)
            else:
                y_zero = np.zeros(xplt.shape, dtype=np.float32)
            ax.plot(xplt, y_zero, color="gray", linestyle="--")
        if empty_plot:
            ax.plot(xstim, ystim, marker="s", color='black', markersize=self.marker)
            ax.legend(self.active_plots[2], title=n_odor, loc=1)
        else:
            lines = [lines for i, lines in enumerate(ax.lines) if i not in (1,2)]
            ax.legend(lines, self.active_plots[2], title=n_odor, loc=1)

    def polar_response(self, ax, g):
        '''
        Subroutine for creating Receptive field polar plot using hours on a clock to map odor responses
        :param g:  which glomerulus to plot
        :param plot_at: location of polar plot.  Typically upper left.
        :return:
        '''
        enable_y_autoscale = False
        # this is to create complete circle of plotting
        active_keys = sorted([key for key in self.active_odors[self.active_odor_set].keys()])
        active_indicies = [key-1 for key in active_keys]

        response_integrals = np.vstack((self.dm.exp_dict["response_integral"][active_indicies,:], self.dm.exp_dict["response_integral"][active_indicies[0],:]))

        g_id = self.dm.glomeruli.index(g)
        pi = np.pi
        #angles = np.asarray([pi/3.0, pi/6.0, 0.0, 11.*pi/6.0, 5.*pi/3., 3.*pi/2., 4.*pi/3., 7.*pi/6., pi, 5.*pi/6., 2.*pi/3., pi/2., pi/3.0,])
        angles = clock_angles(len(active_keys))
        ticks = [str(key) for key in active_keys]
        # get the integrated responses for the glomerulus
        glom_responses = response_integrals[:,g_id]

        #Map negatives to 0 for polar plotting
        #glom_responses = np.where(gr >0, gr, np.zeros_like(gr))

        # Make Scaling reasonable for the plot
        gmax = np.max(glom_responses)
        y_range_integral = self.dm.y_range_integral
        if isinstance(self.y_range_integral, (list, tuple)):
            yrange = self.y_range_integral

        elif self.y_range_integral == "common_autoscale":
            yrange = self.dm.autoscale_integrals
        else:
            enable_y_autoscale = True
            yrange = self.dm.autoscale_integrals
        if enable_y_autoscale:
            if gmax > y_range_integral[1]:
                yrange = [y_range_integral[0], gmax]
            elif gmax > y_range_integral[1] / 2.0:
                yrange = y_range_integral
            elif gmax > y_range_integral[1] / 5.0:
                yrange = [y_range_integral[0] / 2.0, y_range_integral[1] / 2.0]
            else:
                yrange = [y_range_integral[0] / 5.0, y_range_integral[1] / 5.0]

        plt.sca(ax)

        if self.plot_mode == "single":
            use_color = self.colors_id[self.default_color]
            self.polar_max = yrange[1]

        else:
            use_color = self.colors_id[self.color[2]]
            if yrange[1] < self.polar_max:
                yrange[1] = self.polar_max
            else:
                self.polar_max = yrange[1]    # For overlapping plots use the largest y range.
        if len(ax.get_lines()) < 1:
            empty_plot = True
            self.polar_max = yrange[1]
        else:
            empty_plot = False

        ax.plot(angles, glom_responses, color=use_color)
        if empty_plot:
            ax.set_xticks(angles)
            ax.set_xticklabels(ticks)  # 60, 30, 0, 330, 300, 270, 240, 210, 180, 150, 120, 90),

        ax.set_rlim(yrange[0], yrange[1])

    def plot_integrals(self):
        '''
        Fig[1] Subrouting for Plotting integrals for all glomeruli for each odor.
        :return: a matplotlib plot
        '''

        enable_y_autoscale = False

        # Integrate the responses after baseline subtraction
        if "response_integral" not in self.dm.exp_dict.keys():
            self.dm.exp_dict["response_integral"] = self.simple_integral()
        if isinstance(self.figs[1], (plt.Figure)):
            self.figs[1]
        else:
            self.figs[1] = plt.figure(1 + self.fig_shift)
            self.figs[1].canvas.mpl_connect('close_event', self.handle_close)
            self.active_figs[1] = True

        if isinstance(self.y_range_integral, (list, tuple)):
            y_range_integral = self.y_range_integral

        elif self.y_range_integral == "common_autoscale":
            y_range_integral = self.dm.autoscale_integrals
        else:
            enable_y_autoscale = True


        title = "Odor Response Integrals"
        plt.suptitle(title, fontsize=self.title_size)
        plt_dict = self.active_odors[self.active_odor_set]
        for n, odor_info in enumerate(plt_dict.keys()):
            n = n + 1
            if len(self.figs[1].axes) > n-1:
                ax = self.figs[1].axes[n-1]
                plt.sca(ax)
            else:
                ax = plt.subplot(self.plot_grid[0], self.plot_grid[1], n)
                plt.sca(ax)

            plt_id = odor_info
            plt_title = self.dm.odor_labels[odor_info]

            if enable_y_autoscale:
                ax.autoscale(enable=True, axis='y')
            else:
                ax.set_ylim(y_range_integral[0],y_range_integral[1])  # Graph scaling to use
            ax.set_xlim(0, max(self.dm.glomeruli) + 2)
            # create a zero baseline
            xplt = np.asarray(self.dm.glomeruli)
            yzero = np.zeros_like(xplt)

            # extract the correct integrated odor responses
            yplt = self.dm.exp_dict["response_integral"][plt_id-1, :]
            if self.sort_integrals:
                syplt = -np.sort(-yplt)
            ax.plot(xplt, yplt, color="black")
            if self.sort_integrals:
                ax.plot(xplt, syplt, color=self.sorted_color)
            cyplt = np.ones_like(xplt) 
            ax.plot(xplt, cyplt * self.dm.exp_dict["response_cand"]["stats"]["pos_criterion"], linestyle="-", color="blue", linewidth=1.0)
            ax.plot(xplt, cyplt * self.dm.exp_dict["response_cand"]["stats"]["neg_criterion"], linestyle="-", color="blue", linewidth=1.0)
            ax.plot(xplt, yzero, linestyle=":", color="red", linewidth=1.0)
            ax.legend(['Glomerulus', 'Sorted'], title=plt_title, loc=1)
        #plt.pause(.5)

    def make_ef(self, *args, **kwargs):
        self.ef = Explore_f0(self, *args, **kwargs)

    def add_random_colors(self, n=1):
        start_c = max(self.dm.colors_id.keys()) + 1
        for i in range(n):
            self.dm.colors_id[start_c + i] = np.random.rand(3)

    def plot_overlap(self, glom, odor):
        '''
        Divides up the data based on duration of stimuli from 0.5, 1.0, 3.0, 5.0
        :param glom: Which glomerulus to plot
        :param s_dur: stimulus durations
        :param r_data: responses to the stimuli
        :return: figure plot showing the responses
        '''
        if isinstance(odor, (int)):
            odor = self.dm.odor_labels[odor]
        rd = self.dm.exp_dict["response_data"][odor]["overlap"]
        r_data = rd["resp"]
        s_dur = rd["time"]
        divisions = [0.8, 1.5, 3.5]
        s_dur = np.asarray(s_dur)
        r_05 = r_data[s_dur < divisions[0],:]
        r_1 = r_data[np.logical_and(s_dur>divisions[0], s_dur < divisions[1]), :]
        r_3 = r_data[np.logical_and(s_dur>divisions[1], s_dur < divisions[2]), :]
        r_5 = r_data[s_dur > divisions[2],:]
        g_slice = glom-1
        t = [0.5,1.0, 3.0, 5.0]
        r_05m = np.mean(r_05, axis=0)
        r_1m = np.mean(r_1, axis=0)
        r_3m = np.mean(r_3, axis=0)
        r_5m = np.mean(r_5, axis=0)
        m_y = [r_05m[g_slice], r_1m[g_slice], r_3m[g_slice],r_5m[g_slice]]
        r_05s = np.std(r_05, axis=0)
        r_1s = np.std(r_1, axis=0)
        r_3s = np.std(r_3, axis=0)
        r_5s = np.std(r_5, axis=0)
        x_05 = np.ones(r_05.shape[0]) *0.5
        s_y = [r_05s[g_slice], r_1s[g_slice], r_3s[g_slice],r_5s[g_slice]]
        x_1 = np.ones(r_1.shape[0]) * 1.
        x_3 = np.ones(r_3.shape[0]) * 3.
        x_5 = np.ones(r_5.shape[0]) * 5.

        if isinstance(self.figs[5], (plt.Figure)):
            self.figs[5]
        else:
            self.figs[5] = plt.figure(5 + self.fig_shift)
            self.figs[5]._label = "5"
            self.figs[5].canvas.mpl_connect('close_event', self.handle_close)
            self.active_figs[5] = True
        glom_id = self.dm.glomeruli.index(glom)

        ax = plt.subplot(111)

        ax.scatter(x_1, r_1[:,g_slice])
        ax.scatter(x_05, r_05[:,g_slice])
        ax.scatter(x_3, r_3[:,g_slice])
        ax.scatter(x_5, r_5[:,g_slice])
        ax.scatter(t, m_y, marker="o", s=70, color="black")
        ax.errorbar(t, m_y, yerr=s_y, uplims=True, lolims=True, capsize=5, color="black")
        ax.plot(t, m_y, color="black")
        ax.set_title("Glomerulus: " + str(glom)+ "  Odor: 5", fontsize=24)




class Explore_f0():
    _what_i_am = "Explore_f0"
    def __init__(self, dmp, corner_f=None, accelerate=False, *args, **kwargs):
        try:
            id = dmp.what_am_i
        except:
            print("Need to initialized with a Data carrying object not", type(dmp))
            return
        if id == "Plot_dm":
            self.dmp = dmp
            self.dm = self.dmp.dm
        elif id == "DataMan":
            self.dmp = plot_dm(dmp)
            self.dm = self.dmp.dm
        elif id == "DataManFactory":
            self.dmf = dmp
            self.dm_gen = next(self.dmf)
            self.dm = next(self.dm_gen)
            self.dmp = self.plot_dm(self.dm)

        # Default filter for f0 and df/f data
        self.corner_f = corner_f
        # filter raw if True causes raw data to be filtered for the plot
        self.filter_raw = True
        # Whether to decimate points if > 1e6.
        self.accelerate = accelerate
        # How to scale and shift resp and tread data 
        self.scale_f0 = 1.0
        self.scale_df_f = 1.0
        self.slide_f0 = 20.0
        self.slide_df_f = -0.1
        self.explore_f0_xrange = "auto"
        # What color to use for response_data
        self.rd_color = "red"
        # Figure Ids to use for F0 and Df/f plots
        self.f0_fig_id = 6
        self.df_f_fig_id = 7
        self.visible = {"data": True, "resp": True, "tread": True}
        if self.dm.f0_method in self.dm.known_f0_methods:
            true_plts = [self.dm.f0_method]
        else:
            true_plts = [self.dm.known_f0_method]
        for true_p in true_plts:
            self.visible[true_p] = True


        # how to filter respiration and tread data if empty dict no filtering
        self.resp_filter = {"btype": "band_pass", "corner_f": [0.5, 10.]}
        self.tread_filter = {"btype": "low_pass", "corner_f": 0.5}
        f0_methods = copy.deepcopy(self.dm.known_f0_methods)
        if "combo" in f0_methods:
            f0_methods.remove("combo")
        if len(self.dm.f0s) < len(f0_methods):
            self.dm.f0s = self.dm.make_f0(f0_methods)
            self.dm.df_fs = self.dm.make_df_f("all")
        #self.compare_f0_methods(corner_f=corner_f)
        self.r_avg = float(self.dm.raw_fluorescence.mean())
        self.r_amp = float(self.dm.raw_fluorescence.max()) - float(self.dm.raw_fluorescence.min())
        self.prep_resp(**kwargs)
        self.prep_tread(**kwargs)

    @property
    def what_am_i(self):
        return self._what_i_am

    def prep_resp(self, **kwargs):
        '''
                Method to prepare respiratory data for plotting.  Will allow refiltering by passing filter.
        Prep involves filtering and normalizing data and shifting to center at 0.
        :param filter: Filter info: {btype:x, corner_f:y} if {} no filtering.  
        :return: cleaned up respiratory data
        '''
        if "resp_time" not in self.dm.available_data:
            self.dm.prepare_respiration(**kwargs)
        if "resp_time" not in self.dm.available_data:
            return
        resp_trace = self.dm.exp_dict["respiration"]["trace"]
        resp_time = self.dm.exp_dict["respiration"]["time"]
        resp_amp = float(resp_trace.max()) - float(resp_trace.min())
        resp_avg = float(resp_trace.mean())

        self.b = (copy.deepcopy(resp_trace) - resp_avg) / resp_amp
        self.bt = copy.deepcopy(resp_time)


    def prep_tread(self, **kwargs):
        '''
        Method to prepare treadmill data for plotting.  Will allow refiltering by passing filter.
        Prep involves filtering and normalizing data and shifting to center at 0.
        :param filter: Filter info: {btype:x, corner_f:y} if {} no filtering.  
        :return: cleaned up treadmill data
        '''
        if "tread_time" not in self.dm.available_data:
            self.dm.prepare_treadmill(**kwargs)
        if "tread_time" not in self.dm.available_data:
            return
        tread_time = self.dm.exp_dict["treadmill"]["time"]
        tread_velocity = self.dm.exp_dict["treadmill"]["velocity"]
        self.ct = tread_time
        self.cdt = (tread_time[-1] - tread_time[0])/ float(tread_time.shape[0])
        tread_amp = float(tread_velocity.max()) - float(tread_velocity.min())
        tread_avg = float(tread_velocity.mean())
        self.c =(copy.deepcopy(tread_velocity) - tread_avg) / tread_amp
            


    def compare_f0_methods(self, test_glom=5, corner_f=None, return_early=False):
        '''
        used to compare different methods to generate an f0 trace.  Also plots treadmill and respiratory
        data to see if artifacts are removed by method
        :param test_glo: which glom data to plot
        :param corner_f: how to filter f0 data. if None will use self.corner_f
        :param return_early: If True will stop plot before adding checkboxes
        :return: loads DataPlot df_f_raw and F0 arrays and tread and resp data
        '''
        scale = self.scale_f0
        fig_id = self.f0_fig_id

        if isinstance(self.dmp.figs[fig_id], (plt.Figure)):
            fig_f0 = self.dmp.figs[fig_id]
            ax_f0 = fig_f0.axes[0]
            self.rax_f0 = fig_f0.axes[1]
        else:
            fig_f0, ax_f0 = plt.subplots(nrows=1, ncols=1, num=fig_id)
            self.dmp.figs[fig_id] = fig_f0   #plt.figure(fig_id)
            self.dmp.figs[fig_id]._label = str(fig_id)
            self.dmp.figs[fig_id].canvas.mpl_connect('close_event', self.dmp.handle_close)
            self.dmp.active_figs[fig_id] = True
            self.rax_f0 = plt.axes([0.05, 0.5, 0.1, 0.15])
            title = "Compare F0 Methods"
            plt.suptitle(title, fontsize=self.dmp.title_size)
        
        self.rax_f0.clear()
        self.dmp.clear_plot[fig_id] = True
        self.dmp.color[fig_id] = self.dmp.default_color
        ax_f0.clear()
        self.dmp.lines[fig_id] = {}
        self.dmp.active_plots[fig_id] = test_glom
        label = "Glomerulus" + str(test_glom)
        #if self.y_range == "auto":
        #    ax.autoscale(enable=True, axis='y')
        #else:
        #    ax.set_ylim(self.y_range[0],self.y_range[1])  # Graph scaling to use
        # Autoscale x.  The plot can later be changed interactively
        if self.explore_f0_xrange == "auto":
            ax_f0.autoscale(enable=True, axis='x')
        else:
            ax_f0.set_xlim(pq.Quantity(self.dmp.x_range_all[0], self.dmp.time_unit), pq.Quantity(self.dmp.x_range_all[1], self.dmp.time_unit))

        #ax = plt.subplot(111)
        if corner_f is None:
            corner_f = self.dm.corner_f
        else:
            self.dm.corner_f = corner_f
        self.t = self.dm.times
        if self.filter_raw:
            self.r = smooth(copy.deepcopy(self.dm.raw_fluorescence), corner_f=corner_f, dt=self.dm.avg_dt)
        else:
            self.r = copy.deepcopy(self.dm.raw_fluorescence)
        avg = float(self.r.mean())
        r_amp = float(self.r.max()) - float(self.r.min())
        if "data" in self.visible.keys():
            visibility = True
        else:
            visibility = False
        if "resp_time" in self.dm.available_data:
            if "resp" in self.visible.keys():
                visibility = True
            else:
                visibility = False
            l2 = ax_f0.plot(self.bt, (self.b * r_amp * scale) + avg + self.slide_f0, color="darkgoldenrod", label="respiration", visible=visibility)
            self.dmp.lines[fig_id]["respiration"] = l2

        if "tread_time" in self.dm.available_data:
            if "tread" in self.visible.keys():
                visibility = True
            else:
                visibility = False

            l3 = ax_f0.plot(self.ct, (self.c * r_amp * scale) + avg + self.slide_f0, color="black", linewidth=3, label="velocity", visible=visibility)
            self.dmp.lines[fig_id]["velocity"] = l3
        l1 = ax_f0.plot(self.t, self.r[:,test_glom-1], color=self.rd_color, label=label, linewidth=3, visible=visibility)
        self.dmp.lines[fig_id][label] = l1
        self.dmp.color[fig_id] = 3
        for key, value in self.dm.f0s.items():
            print("color", key, self.dmp.colors_id[self.dmp.color[fig_id]])
            if key in self.visible.keys():
                visibility = True
            else:
                visibility = False
            self.dmp.lines[fig_id][key] = ax_f0.plot(self.t, value[:,test_glom-1], color=self.dmp.colors_id[self.dmp.color[fig_id]], label=key, visible=visibility)
            self.dmp.color[fig_id] += 1
        stim_position = avg *0.9 - 0.1

        sp = stim_position
        print("stim_position_f0", sp)
        sw = self.dmp.stim_line_width

        # plots all the odor stimuli as distinct color line segments
        if self.dmp.clear_plot[3]:
            for n, odor_id in sorted(self.dm.odor_labels.items()):
                active_trials = self.dm.pd_odor[self.dm.pd_odor["odorant"] == odor_id]
                start = active_trials["trial_start_time"]
                stop = active_trials["trial_end_time"]
                plt_pts = [[(make_quantity(x1, self.dmp.time_unit),sp), (make_quantity(x2, self.dmp.time_unit),sp)] for x1, x2 in zip(start, stop)]
                lc = mc.LineCollection(plt_pts, color=self.dmp.colors_id[n-1], linewidths=sw)
                ax_f0.add_collection(lc)
            if self.dmp.add_vbars:
                for n, odor_id in sorted(self.dm.odor_labels.items()):
                    active_trials = self.dm.pd_odor[self.dm.pd_odor["odorant"] == odor_id]
                    start = active_trials["trial_start_time"]
                    stop = active_trials["trial_end_time"]
                    plt_pts = [(make_quantity(x1, self.dmp.time_unit), make_quantity(x2, self.dmp.time_unit)) for x1, x2 in zip(start, stop)]
                    for begin, end in plt_pts:
                        ax_f0.axvspan(begin, end, facecolor= self.dmp.colors_id[n-1], alpha=self.dmp.vbar_trans)

        ax_f0.legend(loc=1)
        if return_early:
            return
        plt.subplots_adjust(left=0.2)
        self.plines = [l[0] for l in self.dmp.lines[fig_id].values()]
        self.labels = [str(pline.get_label()) for pline in self.plines]
        self.visibility = [pline.get_visible() for pline in self.plines]
        #self.rax_f0 = plt.axes([0.05, 0.5, 0.1, 0.15])
        self.check = CheckButtons(self.rax_f0, self.labels, self.visibility)

        self.check.on_clicked(self.func1)

    def compare_df_f_methods(self, test_glom=5, corner_f=None, smooth_df_fs=True, return_early=False):
        '''
        used to select an f0 to generate df/f
        :param test: which tests to try
        :param plot_it: Whether to generate plots to better characterize the f0 method
        :return: loads DataPlot df_f_raw and F0 arrays
        '''
        scale = self.scale_df_f
        fig_id = self.df_f_fig_id

        if isinstance(self.dmp.figs[fig_id], (plt.Figure)):
            print("fig already exists, reset figure")
            fig_df_f = self.dmp.figs[fig_id]
            ax_df_f = fig_df_f.axes[0]
            self.rax_df_f = fig_df_f.axes[1]
        
        else:
            print("fig new, setup figure")
            fig_df_f, ax_df_f = plt.subplots(nrows=1, ncols=1, num=fig_id)
            self.dmp.figs[fig_id] = fig_df_f   #plt.figure(fig_id)
            self.dmp.figs[fig_id]._label = str(fig_id)
            self.dmp.figs[fig_id].canvas.mpl_connect('close_event', self.dmp.handle_close)
            self.dmp.active_figs[fig_id] = True
            self.rax_df_f = plt.axes([0.05, 0.5, 0.1, 0.15])
            title = "Compare df/f Methods"
            plt.suptitle(title, fontsize=self.dmp.title_size)
        
        self.dmp.clear_plot[fig_id] = True
        self.dmp.color[fig_id] = self.dmp.default_color
        ax_df_f.clear()
        self.rax_df_f.clear()
        self.dmp.lines[fig_id] = {}
        self.dmp.active_plots[fig_id] = test_glom
        label = "Glomerulus" + str(test_glom)
        #if self.y_range == "auto":
        #    ax.autoscale(enable=True, axis='y')
        #else:
        #    ax.set_ylim(self.y_range[0],self.y_range[1])  # Graph scaling to use
        # Autoscale x.  The plot can later be changed interactively
        if self.explore_f0_xrange == "auto":
            ax_df_f.autoscale(enable=True, axis='x')
        else:
            ax_df_f.set_xlim(pq.Quantity(self.dmp.x_range_all[0], self.dmp.time_unit), pq.Quantity(self.dmp.x_range_all[1], self.dmp.time_unit))

        if corner_f is None:
            corner_f = self.dm.corner_f
        else:
            self.dm.corner_f = corner_f
        df_fs = {}
        if smooth_df_fs:
            for key, value in self.dm.df_fs.items():
                df_fs[key] = smooth(copy.deepcopy(value), corner_f=corner_f, dt=self.dm.avg_dt)
        else:
            df_fs = copy.deepcopy(self.dm.df_fs)
        self.t = self.dm.times

        if self.filter_raw:
            self.r = smooth(copy.deepcopy(self.dm.df_f_raw), corner_f=corner_f, dt=self.dm.avg_dt)
        else:
            self.r = copy.deepcopy(self.dm.df_f_raw)
        avg = float(self.r.mean())
        r_amp = float(self.r.max()) - float(self.r.min())
        if "data" in self.visible.keys():
            visibility = True
        else:
            visibility = False


        if "resp_time" in self.dm.available_data:
            if "resp" in self.visible.keys():
                visibility = True
            else:
                visibility = False

            l2 = ax_df_f.plot(self.bt, (self.b * r_amp * scale) + avg + self.slide_df_f, color="darkgoldenrod", label="respiration", visible=visibility)
            self.dmp.lines[fig_id]["respiration"] = l2

        if "tread_time" in self.dm.available_data:
            if "tread" in self.visible.keys():
                visibility = True
            else:
                visibility = False
            l3 = ax_df_f.plot(self.ct, (self.c * r_amp * scale) + avg + self.slide_df_f, color="black", linewidth=3, label="velocity", visible=visibility)
            self.dmp.lines[fig_id]["velocity"] = l3
        l1 = ax_df_f.plot(self.t, self.r[:,test_glom-1], color=self.rd_color, label=label, linewidth=3, visible=visibility)
        self.dmp.lines[fig_id][label] = l1

        self.dmp.color[fig_id] = 3
        for key, value in df_fs.items():
            print("color", key, self.dmp.colors_id[self.dmp.color[fig_id]] )
            if key in self.visible.keys():
                visibility = True
            else:
                visibility = False
            self.dmp.lines[fig_id][key] = ax_df_f.plot(self.t, value[:,test_glom-1], color=self.dmp.colors_id[self.dmp.color[fig_id]], label=key, visible=visibility)
            self.dmp.color[fig_id] += 1

        # This block is to show where stimuli are present
        if self.dmp.stim_position < 0:
            if isinstance(self.dmp.y_range,(list, tuple)):
                stim_position = self.dmp.y_range[0] * abs(self.dmp.stim_position)
            else:
                stim_position = abs(self.dmp.stim_position)    # np.min(yplt) *
        else:
            if isinstance(self.dmp.y_range,(list, tuple)):
                stim_position = self.dmp.y_range[1] * abs(self.dmp.stim_position)
            else:
                stim_position = abs(self.dmp.stim_position)    # np.max(yplt) *

        sp = stim_position

        sw = self.dmp.stim_line_width
        # plots all the odor stimuli as distinct color line segments
        if self.dmp.clear_plot[fig_id]:
            for n, odor_id in sorted(self.dm.odor_labels.items()):
                active_trials = self.dm.pd_odor[self.dm.pd_odor["odorant"] == odor_id]
                start = active_trials["trial_start_time"]
                stop = active_trials["trial_end_time"]
                plt_pts = [[(make_quantity(x1, self.dmp.time_unit),sp), (make_quantity(x2, self.dmp.time_unit),sp)] for x1, x2 in zip(start, stop)]
                lc = mc.LineCollection(plt_pts, color=self.dmp.colors_id[n-1], linewidths=sw)
                ax_df_f.add_collection(lc)
            if self.dmp.add_vbars:
                for n, odor_id in sorted(self.dm.odor_labels.items()):
                    active_trials = self.dm.pd_odor[self.dm.pd_odor["odorant"] == odor_id]
                    start = active_trials["trial_start_time"]
                    stop = active_trials["trial_end_time"]
                    plt_pts = [(make_quantity(x1, self.dmp.time_unit), make_quantity(x2, self.dmp.time_unit)) for x1, x2 in zip(start, stop)]
                    for begin, end in plt_pts:
                        ax_df_f.axvspan(begin, end, facecolor= self.dmp.colors_id[n-1], alpha=self.dmp.vbar_trans)

        ax_df_f.legend(loc=1)
        if return_early:
            return
        plt.subplots_adjust(left=0.2)
        self.plines = [l[0] for l in self.dmp.lines[fig_id].values()]
        self.labels = [str(pline.get_label()) for pline in self.plines]
        self.visibility = [pline.get_visible() for pline in self.plines]
        self.check = CheckButtons(self.rax_df_f, self.labels, self.visibility)

        self.check.on_clicked(self.func2)

    def func1(self, label):
        #index = self.labels.index(label)
        print("Click", label)
        self.dmp.lines[self.f0_fig_id][label][0].set_visible(not self.dmp.lines[self.f0_fig_id][label][0].get_visible())
        self.dmp.figs[self.f0_fig_id].axes[0].legend(loc=1)
        plt.draw()

    def func2(self, label):
        #index = self.labels.index(label)
        print("Click", label)
        self.dmp.lines[self.df_f_fig_id][label][0].set_visible(not self.dmp.lines[self.df_f_fig_id][label][0].get_visible())
        self.dmp.figs[self.df_f_fig_id].axes[0].legend(loc=1)
        plt.draw()

def clock_angles(n):
    tick = 360// n
    angles = []
    for i in range(0, n+1):
        degree = 90 - tick * i
        if degree >= 0:
            angles.append(math.radians(float(degree)))
        else:
            angles.append(2.0 * np.pi + math.radians(float(degree)))
    return angles

def add_vbars(dmp, fig_num=6, axis_num=0):
    fig = dmp.figs[fig_num]
    ax = fig.axes[0]
    for n, odor_id in sorted(dmp.dm.odor_labels.items()):
        active_trials = dmp.dm.pd_odor[dmp.dm.pd_odor["odorant"] == odor_id]
        start = active_trials["trial_start_time"]
        stop = active_trials["trial_end_time"]
        plt_pts = [(make_quantity(x1, dmp.time_unit), make_quantity(x2, dmp.time_unit)) for x1, x2 in zip(start, stop)]
        for begin, end in plt_pts:
            plt.axvspan(begin, end, facecolor= dmp.colors_id[n-1], alpha=0.2)

def select_grid(n_axes):
    '''
    Not currently implemented but can change plot grid based on number of odor deliveries needed to be plotted
    :param n_axes:
    :return:
    '''
    if n_axes <= 6:
        plot_grid = [2,3]
    elif n_axes <= 9:
        plot_grid = [3,3]
    elif n_axes <= 12:
        plot_grid = [3,4]
    elif n_axes <= 16:
        plot_grid = [4,4]
    else:
        plot_grid = False
    return plot_grid
