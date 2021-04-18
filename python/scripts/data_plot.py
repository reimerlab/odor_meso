import os
import numpy as np
from matplotlib import pyplot, animation
from matplotlib_scalebar.scalebar import ScaleBar
from IPython.display import HTML
from scipy.interpolate import griddata
from tqdm import tqdm

from pipeline import experiment, meso, odor
from DataMan import helpers
from plot_dm import plot_dm


# Wrapper for plot_dm
def plot_dataman(key,
              close_preexisting_figs=True,
              sig_prob=0.99, 
              colors_id={0: "black", 
                         1: "red" , 
                         2: "orange", 
                         3: "green", 
                         4: "blue", 
                         5: "purple", 
                         6: "rosybrown", 
                         7: "teal", 
                         8: "saddlebrown", 
                         9: "magenta", 
                         10: "cornflowerblue", 
                         11: "gold", 
                         12: "lime"},
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
    
    # Default values of args located in uinit.init
    data_plotintegral = (odor.OdorAnalysis.PlotIntegral & key).fetch1()

    # Reconstruct DataMan object
    dm                                                     = type('', (), {})()
    dm.what_am_i                                           = "DataMan"
    dm.expt_id                                             = key['experiment_id'] - 1
    dm.my_paths                                            = helpers.My_paths(paths_init='/data/odor_meso/paths.init')
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


def plot_fields(fields, experiment_id, fps, vmin=1, vmax=1, 
                segmentation_methods=[], threshold_fraction=0.1, mask_sig=None):

    # Determine order of fields according to x position
    experiment_id_query = experiment.ExperimentalIdentifier() & f'experiment_id={experiment_id}'
    fields_sorted = np.argsort((meso.ScanInfo.Field & experiment_id_query).fetch('x'))

    print('Image min intensity: ', [fields[f].min() for f in fields_sorted])
    print('Image max intensity: ', [fields[f].max() for f in fields_sorted])

    # figure, axes = pyplot.subplots(1, len(fields), figsize=(8,8))
    figure, axes = pyplot.subplots(1, 1, figsize=(8,8))

    # Plot fields according to x position
    image = []
    # for a, f in enumerate(fields_sorted):
    f=1
    a=0
    axes.set_title(f'Field {f+1}', color='black')
    axes.axis('off')

    # Calculate aspect ratio of pixels
    px_height, px_width, um_height, um_width = (meso.ScanInfo.Field & \
                                                experiment_id_query & \
                                                f'field={f+1}').fetch1('px_height','px_width','um_height','um_width')
    dx_height = um_height/px_height
    dx_width = um_width/px_width

    image.append(axes.imshow(fields[f][:, :, 0], 
                                vmin=vmin, 
                                vmax=vmax, 
                                aspect=dx_height/dx_width))#axes[a]

    pyplot.close()

    def animate(i):
        # for a, f in enumerate(fields_sorted):
        #     image[a].set_data(fields[f][:, :, i])
        image[0].set_data(fields[f][:, :, i])
    video = animation.FuncAnimation(fig = figure, 
                                    func = animate, 
                                    frames = fields[0].shape[2], 
                                    interval = 1000 / fps)

    if segmentation_methods:
        f=1
        a=0
        # for a, f in enumerate(fields_sorted):
        for segmentation_method in segmentation_methods:
            mask_pixels, mask_weights = (meso.Segmentation.Mask & \
                                            experiment_id_query & \
                                            f'field={f+1}' & \
                                            f'segmentation_method={segmentation_method}' \
                                        ).fetch('pixels', 'weights')

            # Calculate mask weights threshold
            weight_min=1
            weight_max=0
            for i in range(0, np.shape(mask_weights)[0]):
                if np.min(mask_weights[i]) < weight_min:
                    weight_min = np.min(mask_weights[i])
                if np.max(mask_weights[i]) > weight_max:
                    weight_max = np.max(mask_weights[i])

            weight_threshold = threshold_fraction * (weight_max-weight_min) + weight_min
            
            # Plot mask contours
            for i in tqdm(range(0, np.shape(mask_pixels)[0])):
                mask_pixels_threshold = np.delete(mask_pixels[i], 
                                                    np.where(mask_weights[i] < weight_threshold))

                mask_pixels_unravel = np.unravel_index(np.ndarray.astype(mask_pixels_threshold,'int64'),
                                                        np.shape(fields[f])[0:2], 
                                                        order='F')
                mask_image = np.zeros(np.shape(fields[f])[0:2], dtype=np.bool)
                mask_image[mask_pixels_unravel] = True

                #TODO fix color for segmentation 1 vs 2 or significant masks
                if (mask_sig is not None) and (i in mask_sig[f]):
                    axes.contour(mask_image, colors='red', linewidths=0.5)#axes[a]
                elif segmentation_method == 1:
                    axes.contour(mask_image, colors='white', linewidths=0.5)
                elif segmentation_method > 1:
                    axes.contour(mask_image, colors='red', linewidths=0.5)

            # Plot scale bar
            # px_width, um_width = (meso.ScanInfo.Field & experiment_id_query & \
            #                       f'field={f+1}').fetch1('px_width','um_width')
            # dx_width = um_width/px_width

            # scalebar = ScaleBar(dx=dx_width, 
            #                     units='um', 
            #                     dimension='si-length', 
            #                     scale_loc='bottom',
            #                     fixed_value = 100, 
            #                     color='w', 
            #                     box_alpha=0, 
            #                     pad=0.5,
            #                     location='lower right')
            # axes[a].add_artist(scalebar)

    return HTML(video.to_html5_video(embed_limit=100))


def plot_stitched_movie(experiment_id_query, image, fps, vmin=1, vmax=1, masks=False, segmentation_method=1, threshold_fraction=0.1):
    figure = pyplot.figure(figsize=(8,8))
    
    print('Image min intensity: ', np.nanmin(image))
    print('Image max intensity: ', np.nanmax(image))
    
    im1 = pyplot.imshow(image[:, :, 0], 
                        vmin=np.nanmin(image)*vmin, 
                        vmax=np.nanmax(image)*vmax, 
                        aspect="equal")

    figure.patch.set_facecolor('xkcd:black')
    axis = figure.axes[0]
    axis.tick_params(axis='both', colors='white')

    pyplot.close()

    def animate(i):
        im1.set_data(image[:, :, i])

    video = animation.FuncAnimation(fig = figure, 
                                    func = animate, 
                                    frames = image.shape[2], 
                                    interval = 1000 / fps)

    if masks:
        #     print('Weight min:       ', weight_min)
        #     print('Weight max:       ', weight_max)
        #     print('Weight threshold: ', weight_threshold)

            # Move masks to stitched image space
        #     key = {'animal_id': 1571, 'session': 1, 'scan_idx': 1, 'pipeline_version': 1}
        #         all_heights = (meso.ScanInfo.Field & key).fetch('px_height')
        #         all_widths = (meso.ScanInfo.Field & key).fetch('px_width')
        #         length_total = sum(np.multiply(all_heights, all_widths))
        #         all_coordinates = np.empty((length_total, 2))
        #         all_values = np.empty((length_total))
        #         counter = 0

        for field in [1,2,3]:
            scan_field = (meso.ScanInfo.Field & experiment_id_query & f'field={field}').fetch1()
            image_size = [scan_field['px_height'], scan_field['px_width']]

            um_per_px_width = scan_field['um_width'] / scan_field['px_width']
            um_per_px_height = scan_field['um_height'] / scan_field['px_height']

            mask_pixels, mask_weights = (meso.Segmentation.Mask 
                                        & f'field={field}' 
                                        & f'segmentation_method={segmentation_method}'
                                        ).fetch('pixels', 'weights')

            # Calculate mask weights threshold
            weight_min=1
            weight_max=0
            for i in range(0, np.shape(mask_weights)[0]):
                if np.min(mask_weights[i]) < weight_min:
                    weight_min = np.min(mask_weights[i])

                if np.max(mask_weights[i]) > weight_max:
                    weight_max = np.max(mask_weights[i])

            weight_threshold = threshold_fraction * (weight_max-weight_min) + weight_min

    #             mask_image = np.zeros(image_size, dtype=int)

            # Plot each mask
            for i in tqdm(range(0, np.shape(mask_pixels)[0])):
                mask_threshold = np.delete(mask_pixels[i], np.where(mask_weights[i] < weight_threshold))
                mask_unravel = np.unravel_index(np.ndarray.astype(mask_threshold,'int64'), image_size, order='F')
    #                 mask_image[mask_pixels_unravel] = 1

        #     x = um_per_px_width*(np.array(range(scan_field['px_width'])) - scan_field['px_width']/2)
        #     y = um_per_px_height*(np.array(range(scan_field['px_height'])) - scan_field['px_height']/2)

            # Motion correct
        #     x_shifted = np.repeat(x[...,np.newaxis], frame+1-frame, axis=1) + x_shifts[frame:frame+1]
        #     y_shifted = np.repeat(y[...,np.newaxis], frame+1-frame, axis=1) + y_shifts[frame:frame+1]

                mask_transform = np.empty((len(mask_threshold), 2))
                mask_values = np.ones(len(mask_threshold))

                # Stitch masks
    #             for j in range(len(mask_threshold)):
                mask_transform[:,0] = um_per_px_width*(mask_unravel[1][:]  - (scan_field['px_width']/2)) + scan_field['x']
                mask_transform[:,1] = um_per_px_height*(mask_unravel[0][:] - (scan_field['px_height']/2)) + scan_field['y']
                # print(np.shape(mask_pixels[i]), np.shape(mask_threshold), np.shape(mask_unravel), np.shape(mask_transform))
    #             print(i, mask_unravel, mask_transform)
    #                 print(i, mask_unravel[0][i], mask_unravel[1][i], mask_transform[i,0], mask_transform[i,1])

    #                 all_coordinates = np.empty((len(mask_pixels_threshold), 2))
    #                 all_values = np.empty((len(mask_pixels_threshold)))
    #                 counter = 0

    #                 for i in range(scan_field['px_width']):
    #                     for j in range(scan_field['px_height']):
    #                         all_coordinates[counter,0] = um_per_px_width*(i - (scan_field['px_width']/2)) + scan_field['x']
    #                         all_coordinates[counter,1] = um_per_px_height*(j - (scan_field['px_height']/2)) + scan_field['y']
    #                         all_values[counter] = mask_image[j, i]*1

    #                         counter+=1

                # Interpolate fields onto a uniform grid
                width_min = 108.12#np.min(all_coordinates[:,0]) #108.12
                width_max = 1903.12#np.max(all_coordinates[:,0]) #3108.12
                height_min = -1070.208#np.min(all_coordinates[:,1]) #-1070.208
                height_max = 1725.1253333333334#np.max(all_coordinates[:,1]) #1729.792
                # print(width_min,width_max,height_min,height_max)
                # print(np.min(mask_transform[:,0]),
                #       np.max(mask_transform[:,0]),
                #       np.min(mask_transform[:,1]),
                #       np.max(mask_transform[:,1]))
                print(np.shape(mask_transform), np.shape(mask_values))
                grid_x, grid_y = np.mgrid[height_min:height_max:560j, width_min:width_max:360j]
                mask_image = griddata(mask_transform, mask_values, (grid_y, grid_x), method='linear')

    #                 mask_image = griddata(all_coordinates, all_values, (grid_y, grid_x), method='nearest')
                
                axis.contour(mask_image, colors='white', linewidths=0.5)
    #         mask_image_threshold = np.ma.masked_where(mask_image < 1, mask_image)
    #         axis.imshow(mask_image_threshold, aspect="equal", interpolation='none', alpha=0.5)
    
    scalebar = ScaleBar(dx=5, units='mm', dimension='si-length', length_fraction=0.111, color='w', box_alpha=0, location='lower right')
    axis.add_artist(scalebar)

    return HTML(video.to_html5_video(embed_limit=100))


# # Determine complete image extent
# for field in [1,2,3]:
#     scan_field = (meso.ScanInfo.Field & key & f'field={field}').fetch1()
#     um_per_px_width = scan_field['um_width'] / scan_field['px_width']
#     um_per_px_height = scan_field['um_height'] / scan_field['px_height']

#     width_min = um_per_px_width*(0 - scan_field['px_width']/2) + scan_field['x']
#     width_max = um_per_px_width*(360 - scan_field['px_width']/2) + scan_field['x']
#     height_min = um_per_px_height*(0 - scan_field['px_height']/2) + scan_field['y']
#     height_max = um_per_px_height*(560 - scan_field['px_height']/2) + scan_field['y']
#     print('Complete image coordinates: ', width_min, width_max, height_min, height_max)
    
# Complete image coordinates:  708.12 2508.12 -1070.208 1729.792
# Complete image coordinates:  108.12 1908.12 -1070.208 1543.1253333333334 # why is this different?
# Complete image coordinates:  1308.12 3108.12 -1070.208 1729.792

def plot_traces(odor_index, plot_integral, fluorescence, ylim_f=[0,1000], ylim_df_f=[-1,1]):
    for i in plot_integral['glomeruli']:
        if (plot_integral['response_integral'][odor_index-1][i-1] < plot_integral['neg_criterion']) or \
        (plot_integral['response_integral'][odor_index-1][i-1] > plot_integral['pos_criterion']):
            
            fig = pyplot.figure(figsize=(12,10))
            fig.suptitle(f'Mask {i} with signficant repsonse \
                ({plot_integral["response_integral"][odor_index-1][i-1]:.2f}) \
                    to odor ({plot_integral["odor_labels"][odor_index]})')

            ax1 = pyplot.subplot(2, 3, 1)
            pyplot.plot(fluorescence['raw_fluorescence'][:,i-1])
            ax1.title.set_text('Raw Fluorescence')
            ax1.set_ylim(ylim_f)
            ax1.set_xlabel('Frame (#)')
            
            ax2 = pyplot.subplot(2, 3, 2)
            pyplot.plot(fluorescence['f0s']['constant_med'][:,i-1])
            ax2.title.set_text('F0 - constant_med')
            ax2.set_ylim(ylim_f)
            ax2.set_xlabel('Frame (#)')

            ax3 = pyplot.subplot(2, 3, 3)
            pyplot.plot(fluorescence['f0s']['bkg_smooth'][:,i-1])
            ax3.title.set_text('F0 - bkg_smooth')
            ax3.set_ylim(ylim_f)
            ax3.set_xlabel('Frame (#)')
            
            ax4 = pyplot.subplot(2, 3, 4)
            pyplot.plot(fluorescence['df_fs']['constant_med'][:,i-1])
            ax4.title.set_text('dF/F - constant_med')
            ax4.set_ylim(ylim_df_f)
            ax4.set_xlabel('Frame (#)')

            ax5 = pyplot.subplot(2, 3, 5)
            pyplot.plot(fluorescence['df_fs']['bkg_smooth'][:,i-1])
            ax5.title.set_text('dF/F - bkg_smooth')
            ax5.set_ylim(ylim_df_f)
            ax5.set_xlabel('Frame (#)')

            ax6 = pyplot.subplot(2, 3, 6)
            pyplot.plot(fluorescence['df_fs']['combo'][:,i-1])
            ax6.title.set_text('dF/F - combo')
            ax6.set_ylim(ylim_df_f)
            ax6.set_xlabel('Frame (#)')

