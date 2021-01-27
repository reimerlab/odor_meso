import os
from matplotlib import pyplot
from pipeline import odor, experiment
from DataMan import helpers
from plot_dm import plot_dm

# Default values of args located in uinit.init
def data_plot(key,
              close_preexisting_figs=True,
              sig_prob=0.99, 
              colors_id={ 0: "black", 1: "red" , 2: "orange", 3: "green", 4: "blue", 5: "purple", 6: "rosybrown", 7: "teal", 8: "saddlebrown", 9: "magenta", 10: "cornflowerblue", 11: "gold", 12: "lime"},
              default_color=0,
              title_size=24,
              y_range='common_autoscale',
              plot_grid=[ 3, 4 ],
              plot_mode='single',
              marker=3,
              polar_summary=True,
              plt_baselinesubtracted=True,
              plot_zero_line=True,
              y_range_integral='common_autoscale',
              sort_integrals=True,
              sorted_color='darkorange',
              stim_line_width=5,
              stim_position=-0.09,
              save_figs='replot',
              response_mask='odor_post',
              time_unit='sec'):

     data_plotintegral = (odor.OdorAnalysis.PlotIntegral & key).fetch1()
     paths_init = os.path.join('/data/odor_meso/external/dataplot_storage','paths.init')#os.environ.get('DATAPLOT_STORAGE')

     # Reconstruct DataMan object
     dm                                                     = type('', (), {})()
     dm.what_am_i                                           = "DataMan"
     dm.expt_id                                             = key['experiment_id'] - 1
     dm.my_paths                                            = helpers.My_paths(paths_init=paths_init)
     dm.n_unique_odor_combos                                = data_plotintegral['n_unique_odor_combos']
     dm.autoscale_integrals                                 = data_plotintegral['autoscale_integrals']
     dm.odor_labels                                         = data_plotintegral['odor_labels']
     dm.glomeruli                                           = data_plotintegral['glomeruli']
     dm.exp_dict                                            = {}
     dm.exp_dict["response_integral"]                       = data_plotintegral['response_integral']
     dm.exp_dict['response_cand']                           = {}
     dm.exp_dict['response_cand']['stats']                  = {}
     dm.exp_dict["response_cand"]["stats"]["pos_criterion"] = data_plotintegral['pos_criterion']
     dm.exp_dict["response_cand"]["stats"]["neg_criterion"] = data_plotintegral['neg_criterion']
     dm.close_preexisting_figs                              = close_preexisting_figs
     dm.colors_id                                           = colors_id
     dm.default_color                                       = default_color
     dm.title_size                                          = title_size
     dm.y_range                                             = y_range
     dm.plot_grid                                           = plot_grid
     dm.plot_mode                                           = plot_mode
     dm.marker                                              = marker
     dm.polar_summary                                       = polar_summary
     dm.plt_baselinesubtracted                              = plt_baselinesubtracted
     dm.plot_zero_line                                      = plot_zero_line
     dm.y_range_integral                                    = y_range_integral
     dm.sort_integrals                                      = sort_integrals
     dm.sig_prob                                            = sig_prob
     dm.sorted_color                                        = sorted_color
     dm.stim_line_width                                     = stim_line_width
     dm.stim_position                                       = stim_position
     dm.save_figs                                           = save_figs
     dm.time_unit                                           = time_unit
     dm.response_mask                                       = response_mask

     op = plot_dm(dm)
     op.plot_integrals()
     figure = pyplot.gcf()
     figure.set_size_inches(14, 10)

