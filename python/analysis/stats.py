import numpy as np
#from filters import smooth
#from scipy.stats import ttest_ind as tt
#from scipy.stats import normaltest as nt
from scipy.stats import mannwhitneyu as mwu
from scipy.stats import wilcoxon as wsr
from scipy.stats import norm, zscore
import copy
import matplotlib.pyplot as plt
#import matplotlib as mpl
import quantities as pq

def mann_whitney_u(exp_dict, glom, odors="all", n_compare=1, anal=("odor", "bkg"), use_continuity=True, tails=2, testing=False, sig_p=1.0e-4):
    '''
    Runs mann whitney u test for one glomerulus for a specific odors vs a background set of data of similar construction
    from the same glomerulus.  Uses an average dt value for all comparisions rather than exact experimental time.
    :param exp_dict:  source of data and masks
    :param odors: which odors to test, or all.  Takes lists, tuples, int, or "all"
    :param glom: which glomerulus to test
    :param n_compare: how many random background signals to compare to
    :param bkg_anal:  whether to analyze the background consistency or the odor vs background
    :param use_continuity: a scipy.stats.mann_whitney_u parameter
    :param tails: chooses value for scipy.stats.mann_whitney_u alternative parameter
    :param testing: a debug flag.  can also trigger early returns
    :return: dictionary of results
    '''


    if tails == 2:
        alternative = 'two-sided'
    elif tails == 1:
        alternative = 'greater'
    if odors == "all":
        odors = [n for n, od in exp_dict["odors"]]
    elif isinstance(odors, (int)):
        odors = [odors]

    mask_dict = exp_dict["mask_dict"]
    glom_data = exp_dict["sdata-leak"][:,glom-1]
    #if testing: print("loading", mask_dict["pre_maxodor_post_mask"][mask_dict["pre_maxodor_post_mask"]==-1].shape, glom_data.shape)

    p_dict = {}
    stat_dict = {}
    cp_dict = {}
    cstat_dict = {}
    for odor in odors:
        stat_list = []
        p_list = []
        stat_np = np.ones([n_compare], dtype=np.float32)
        p_np = np.ones([n_compare], dtype=np.float32)
        c_p_np = np.ones([n_compare,n_compare], dtype=np.float32)
        c_stat_np = np.ones([n_compare,n_compare], dtype=np.float32)
        # Generate control comparisons to determine significance
        c_stat_list = []
        c_p_list = []
        odor_mask = mask_dict["odor_post_mask_by_odor"][odor]
        response_data = glom_data[odor_mask<0]
        # Grab control data from region outside [odor_start, odor_end+buffer]
        ctrl_masks = match_stimuli(mask_in=odor_mask, stimuli=exp_dict["stimuli"], odor=odor, dt=exp_dict["avg_dt"], n_masks=n_compare,
                                   pre_buffer=0.0, buffer=exp_dict["buffer"])
        if "odor" in anal:
            for i, mask in enumerate(ctrl_masks):
                y = glom_data[mask==-(odor+1)]
                if testing: print("inputs", response_data.shape, y.shape, response_data.mean(), response_data.std(), y.mean(), y.std(), mask.shape)
                stat, p = mwu(response_data, y, use_continuity=use_continuity, alternative=alternative)
                stat_list.append([i, stat])
                stat_np[i] = stat
                p_np[i] = p
                p_list.append([i,p])
            p_dict[odor] = p_np       #p_list
            stat_dict[odor] = stat_np     #stat_list
        if n_compare > 1 and "bkg" in anal:
            for i in range(n_compare-1):
                for j in range(i, n_compare):
                    cx = glom_data[ctrl_masks[i]==-(odor+1)]
                    cy = glom_data[ctrl_masks[j]==-(odor+1)]
                    if testing: print("inputs", cx.shape, cy.shape, cx.mean(), cx.std(), cy.mean(), cy.std(), mask.shape)
                    cstat, cp = mwu(cx, cy, use_continuity=True, alternative=alternative)
                    c_stat_list.append([i,j,cstat])
                    c_stat_np[i,j] = cstat
                    c_stat_np[j,i] = cstat
                    c_p_np[i,j] = cp
                    c_p_np[j,i] = cp
                    c_p_list.append((i,j, cp))
            cstat_dict[odor] = c_stat_np      #c_stat_list
            cp_dict[odor] = c_p_np
            #print("ctrl shapes ", cp_dict[odor].shape, cstat_dict[odor].shape )    #c_p_list
    mwu_dict = {"odor": odor, "glom": glom, "stat": stat_dict, "p": p_dict, "n_comp": n_compare, "ctrl_stat": cstat_dict,
                "ctrl_p": cp_dict, "odor_mask":odor_mask, "ctrl_masks": ctrl_masks}

    return mwu_dict

def match_stimuli(dm, odor, n_masks=1, mask_region="odor_post", match_region="odor_post"):
    '''
    if fed a mask with 1 values in where stimuli are present, will find places for new mask.slices that is non-overlapping
    values inserted will be odor+1, so odor 1 will have mask values 2.  If match region is "common", then all stimuli are
    treated the same and a common slice will be used regardless of the true odor duration.

    :param dp: at DataPlot object
    :param odor: Which odor to create matched stimuli for as [key#, odor_name].
    :param n_masks: how many matched stimuli to create
    :param mask_region: what part of the data to exclude as likely containing odor responses
    :param match_region: what region to match.  Can be different from mask_region
    :return: a mask dictionary and a slices dictionary containing the matched periods from outside the mask_region
    '''

    cslices = {}
    cmasks = {}
    dm_mask, __ = dm.mask_slice_manager(how_to=mask_region, by_odor=odor[1])
    mask = np.abs(copy.deepcopy(dm_mask))   # Need to copy mask in otherwise original mask is modified.
    n_times = mask.shape[0]
    __, slices = dm.mask_slice_manager(how_to=match_region, by_odor=odor[1])
    if "common" in match_region:
        __, odor_slices = dm.mask_slice_manager(how_to="common", by_odor=odor[1])

    else:
        __, odor_slices = dm.mask_slice_manager(how_to=mask_region, by_odor=odor[1])

    for n in range(n_masks):
        cmask = np.zeros_like(mask)
        total_duration = 0
        total_odor_duration = 0
        cresp_slices = []
        codor_slices = []
        for (start, stop), (odor_start, odor_stop) in zip(slices, odor_slices):
            not_free = True
            duration = stop - start
            relative_start_odor = odor_start - start
            shift_needed = max(0, relative_start_odor)
            odor_duration = odor_stop - odor_start
            n_tries = 0
            while not_free:
                if n_tries >= 5000:
                    print("couldn't find a good location")
                    return False
                n_start = int(np.random.randint(shift_needed, n_times - duration, size=1))
                if np.abs(mask[slice(n_start, n_start + duration)]).max() < 0.5:
                    not_free = False
                n_tries += 1
            cmask[slice(n_start, n_start + duration)] = odor[0] + 1  #
            mask[slice(n_start, n_start + duration)] =  odor[0] + 1  # 1 is masking real odor stimuli so add 1 to odor key# to be unique
            cresp_slices.append([n_start, n_start + duration])
            codor_slices.append([n_start + relative_start_odor, n_start + relative_start_odor + odor_duration])
            total_duration += duration
            total_odor_duration += odor_duration
        #print('check_mask', total_duration, np.where(mask_in==-(odor + 1))[0].shape[0], np.where(mask==-(odor + 1))[0].shape[0])
        #masks[n] = mask
        cmasks[n] = cmask
        cslices[n] = {"response": cresp_slices, "odor": codor_slices, "duration": total_duration, "odor_duration": total_odor_duration}
    return cmasks, cslices

def anal_mwu_dict(mwu_dict):
    c_p = mwu_dict["ctrl_p"]
    ctrl_ps = {}
    odor_ps = {}
    odor_ps_z = {}
    odor_ps_p = {}

    for pkey, pvalue in c_p.items():
        nrow, ncol = pvalue.shape
        i_mask = np.eye(ncol, dtype=np.int)
        ps = []
        for i in range(nrow):
            pst = pvalue[i_mask[:,i] !=1,i]
            pstc = np.where(pst>0.0, pst, 1e-200)
            ps.append(np.median(-np.log10(pstc)))
        ps_np = np.array(ps)
        ps_mean = ps_np.mean()
        ps_std = ps_np.std()
        ctrl_ps[pkey] = [ps_mean, ps_std]
    odor_p = mwu_dict["p"]
    for pkey, pvalue in odor_p.items():
        p_set = np.where(odor_p[pkey]>1e-100, odor_p[pkey], 1e-100)
        odor_ps[pkey] = np.median(np.log10(p_set))
        odor_ps[pkey] = np.where(odor_ps[pkey] is np.inf, -200.0, -odor_ps[pkey])
        print("min, max", pkey, odor_ps[pkey].min(), odor_ps[pkey].max())

        odor_ps_z[pkey] = (odor_ps[pkey] - ctrl_ps[pkey][0])/ctrl_ps[pkey][1]
        odor_ps_p[pkey] = norm.cdf(odor_ps_z[pkey])
        odor_ps_p[pkey] = -np.log10(np.where(odor_ps_p[pkey] > 1.0, 1.0, odor_ps_p[pkey]))
        #print("odor_ps", pkey, odor_ps[pkey].min(), odor_ps[pkey].max(), np.median(odor_ps[pkey]))

    mwu_dict["odor_ps_to_ctrl-ps"] = odor_ps_p
    mwu_dict["odor_median-ps"] = odor_ps
    mwu_dict["ctrl_ps"] =  ctrl_ps

    return mwu_dict

def analyze_odor_specificity(expt_dict, h_map=None, gloms="all", compare=10, test="mwu"):
    '''
    This method performs a repeated analysis of significance using mann_whitney_u method
    "mwu" or "t" test, then outputs 2d numpy arrays with median p for odors and background

    :param expt_dict:
    :param gloms:
    :param compare:
    :return:
    '''
    resp_dict={}
    p_summary = []
    cp_summary = []

    if gloms == "all":
        gloms = expt_dict["glomeruli"]
    elif isinstance(gloms, (int)):
        gloms = [gloms]

    n_odors = len(expt_dict["odors"])
    print("Analyzed Glomerulus ", end="", flush=True)

    for glom in gloms:
        resp_dict[glom] = mann_whitney_u(expt_dict, glom, odors="all", n_compare=compare)
        resp_dict[glom] = anal_mwu_dict(resp_dict[glom])
        p_temp = []
        cp_temp = []
        for odor in range(n_odors):
            p_temp.append(resp_dict[glom]["odor_median-ps"][odor+1])
            cp_temp.append(resp_dict[glom]["odor_ps_to_ctrl-ps"][odor+1])
        p_summary.append(np.asarray(p_temp))
        cp_summary.append(np.asarray(cp_temp))
        print(str(glom), end=" ", flush=True)
        if not glom % 20: print("")
    pt = np.vstack(p_summary)
    cpt = np.vstack(cp_summary)
    p_all = -np.log10(np.where(pt>0.0, pt, 1e-200))
    cp_all = -np.log10(np.where(cpt>0.0, cpt, 1e-200))
    h_map = plt_heatmaps(p_all, cp_all, expt_dict["odors"], h_map, scale=100., fsize=6)
    return resp_dict, p_all, cp_all, h_map

def analyze_glom_specificity(expt_dict, h_map=None, gloms="all", compare=10, test="mwu"):
    '''
    This method performs a repeated analysis of significance using mann_whitney_u method
    "mwu" or "t" test, then outputs 2d numpy arrays with median p for odors compared to bkg

    :param expt_dict:
    :param gloms:
    :param compare:
    :return:
    '''
    resp_dict={}
    p_summary = []
    p_to_c_summary = []

    if gloms == "all":
        gloms = expt_dict["glomeruli"]
    elif isinstance(gloms, (int)):
        gloms = [gloms]

    n_odors = len(expt_dict["odors"])
    print("Analyzed Glomerulus ", end="", flush=True)

    for glom in gloms:
        resp_dict[glom] = mann_whitney_u(expt_dict, glom, odors="all", n_compare=compare)
        resp_dict[glom] = anal_mwu_dict(resp_dict[glom])
        p_temp = []
        p_to_c_temp = []
        # Stacking the results for all glomeruli into two arrays mwu_dict["odor_ps_to_ctrl-ps"] = odor_ps_p
        #     mwu_dict["odor_median-ps"] = odor_ps
        for odor in range(n_odors):
            p_temp.append(resp_dict[glom]["odor_ps_to_ctrl-ps"][odor+1])
            p_to_c_temp.append(resp_dict[glom]["odor_median-ps"][odor+1])
        p_summary.append(np.asarray(p_temp))
        p_to_c_summary.append(np.asarray(p_to_c_temp))
        print(str(glom), end=" ", flush=True)
        if not glom % 20: print("")
    pt = np.vstack(p_summary)
    p_to_ct = np.vstack(p_to_c_summary)
    #p_all = -np.log10(np.where(pt>0.0, pt, 1e-200))
    #cp_all = -np.log10(np.where(cpt>0.0, cpt, 1e-200))
    h_map = plt_heatmaps(pt, p_to_ct, expt_dict["odors"], h_map, scale=100., fsize=6)
    return resp_dict, pt, p_to_ct, h_map


def plt_heatmaps(pdata, opdata, odors, h_map=None, scale=100., cmap="OrRd", fsize=16, aspect=0.5):
    '''
    Plots a heatmap for the statistical analyses.  provides a reference for the figure that
    can be pickled
    :param pdata:
    :param opdata:
    :param odors:
    :param scale:
    :param cmap:
    :param fsize:
    :return:
    '''
    if h_map is None or not isinstance(h_fig, (plt.Figure)):
        h_map = plt.figure(4)
        ax1 = plt.subplot(121)
        ax2 = plt.subplot(122)
    else:
        ax1 = h_map.axes[0]
        ax1 = h_map.axes[1]
    rows, cols = pdata.shape
    print("rows, cols", rows, cols)
    ax1.imshow(pdata / scale, cmap=cmap, vmin=0.0, vmax=1.0)
    rows, cols = pdata.shape
    y_labels = []
    y_labels.extend([str(i) for i in range(1, rows + 1)])
    x_labels = [lab for i, lab in odors]
    print(y_labels, x_labels)
    ax1.set_ylim(float(rows) - 0.5, -0.5)
    ax1.set_yticks(np.arange(0.0,rows+1))

    ax1.set_yticklabels(y_labels, fontsize=fsize)
    ax1.set_xlim(-0.5,12.5)
    ax1.set_xticks(np.arange(0,13))
    ax1.set_xticklabels(x_labels, fontsize=fsize, rotation="vertical")
    forceAspect(ax1, aspect)
    ax2.imshow(opdata / scale, cmap=cmap, vmin=0.0, vmax=1.0)
    rows, cols = opdata.shape
    y_labels = [""]
    y_labels.extend([str(i) for i in range(1, rows + 1)])
    x_labels = [lab for i, lab in odors]
    ax2.set_yticklabels(y_labels,fontsize=fsize)
    ax2.set_xlim(-0.5,11.5)
    ax2.set_xticks(np.arange(0,12))
    ax2.set_xticklabels(x_labels, fontsize=fsize, rotation="vertical")
    forceAspect(ax2, aspect)

    plt.show()
    return h_map



def test_print(n):
    print("Analyzed Glomerulus ", end="")

    for i in range(n):

        print(str(i), end=" ", flush=True)

def forceAspect(ax,aspect):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

def adj_fonts(ax, which_ax="x", size=14):
    if "x" in which_ax:
        for tick in ax.get_xticklabels():
            tick.set_fontsize(size)
    if "y" in which_ax:
        for tick in ax.get_yticklabels():
            tick.set_fontsize(size)
    plt.pause(0.5)

def vcorrcoef(row_vectors,vector, mask=None):
    '''
    Computes correlation coefs between a row vector y with n elements, and a matrix of row vectors n x m to return m corr coef
    :param X:
    :param y:
    :return:
    '''
    if mask is not None:
        t_row_vectors = row_vectors[mask>0.5,:].T
        vector = np.expand_dims(vector[mask>0.5], axis=0)
    else:
        t_row_vectors = row_vectors.T
        vector = np.expand_dims(vector, axis=0)
    row_vector_means = np.reshape(np.mean(t_row_vectors, axis=1), (t_row_vectors.shape[0], 1))
    vector_mean = np.mean(vector)
    r_numerator = np.sum((t_row_vectors - row_vector_means)*(vector - vector_mean), axis=1)
    r_denominator = np.sqrt(np.sum((t_row_vectors - row_vector_means)**2, axis=1)*np.sum((vector - vector_mean)**2))
    r_coef = r_numerator/r_denominator
    return r_coef

def twod_vcorr(vin, masks=None, exclusion_mask=True, criterion=0.5):
    '''
    Takes a two_d array, typically [datapoints, glomeruli] and does cross correlation across the data_point vectors.
    Returns a dictionary with glomeruli with significant correlations with other glomeruli,  Analysis can be restricted
    by providing a [1,0] mask for the regions to restrict the analysis to
    :param vin:  a 2d array of vectors to calculate cross-correlations.
    :param masks:  masking arrays for data set to restrict for the analysis
    :param exclusion_mask:  whether to calculate correlation for data in region outside the masked regions
    :param criterion:  A value the x-correlation needs to exceed to be added to the return_dict
    :return:  All glomeruli that are above criterion for x-correlation with this glomerulus.  Does not include itself.
    '''
    glom_corr = {}
    return_dict = {}
    # mask is a 1d array of 0,1 where data is only used if mask ==1.  Not a true mask, but whatever.
    # Running deepcopy to make sure original vectors are not still linked
    if masks is None:
        masks = {"no_mask": np.ones_like(vin.shape[0], dtype=np.int)}
    elif exclusion_mask:   #Creates a mask that excludes points that may be part of any odor response
        a_mask = next(iter(masks.values()))  # Takes the first mask in the masks_dict for construction of exclusion masks
        exclusion_mask = np.zeros_like(a_mask)
        for mask in masks.values():  # Marks all areas where any mask is non-zero
            exclusion_mask += mask
        exclusion_mask = np.logical_not(exclusion_mask)  # Inverts the exclusion mask to mark all areas that were zero
        masks["exclusion_mask"] = exclusion_mask
    for key, mask in masks.items():
        if mask.max() < 0.5:  #Make sure it is a [0,1] mask, if not invert
            mask = -1 * mask
        print("inputsin v mask", vin.shape, mask.shape)
        v = copy.deepcopy(np.array(vin[mask>0.5,:]))
        print("input v", v.shape)
        # The means of the datapoint vectors
        vm = np.mean(v, axis=0)
        # The difference from the mean for each vector
        v_diff = (v - vm)
        # Expand dims so broadcasting during matmul works correctly.  Otherwise this is the square of v_diff
        denom1 = np.expand_dims(np.sum(v_diff * v_diff, axis=0), axis=0)
        print("denom1", denom1.shape)
        # This is the formation of all the cross multiplied v_diff vectors
        neum = np.matmul(v_diff.T, v_diff)
        print("neum", neum.shape)
        # Creates a square denom matrix that will be divided into the neum matrix after taking the sq root.
        denom = np.matmul(denom1.T, denom1)
        print("denom", denom.shape)
        rcoef = neum / np.sqrt(denom)
        # Subtracting 1 from the diagonal to remove self-self correlations.
        no_self = rcoef - np.eye(rcoef.shape[0], rcoef.shape[1])
        # Finding glomeruli that correlate above criterion and storing in a dict with the correlation coef.
        # +1 is used because the np arrays start at 0 but glomeruli index starts at 1.
        for i in range(no_self.shape[0]):
            temp = no_self[:,i]
            g_xcor = np.where(temp>criterion)[0]
            gc_coef = temp[g_xcor]
            glom_corr[i+1] = [[g_x+1, coef] for g_x, coef in zip(g_xcor, gc_coef)]
        # Removes dictionary elements where there are no significant correlations.
        sig_gloms = {key: value for key, value in glom_corr.items() if len(value) > 0}
        return_dict[key] = {"r_gs": rcoef, "sig_gs": sig_gloms, "sig_level": criterion}
    return return_dict

def twod_xcorr(odr, criterion=0.5):
    '''
    Takes a two_d array, typically [datapoints, glomeruli] and does cross correlation across the data_point vectors.
    Returns a dictionary with glomeruli with significant correlations with other glomeruli,  Analysis can be restricted
    by providing a [1,0] mask for the regions to restrict the analysis to
    :param vin:  a 2d array of vectors to calculate cross-correlations.
    :param masks:  masking arrays for data set to restrict for the analysis
    :param exclusion_mask:  whether to calculate correlation for data in region outside the masked regions
    :param criterion:  A value the x-correlation needs to exceed to be added to the return_dict
    :return:  All glomeruli that are above criterion for x-correlation with this glomerulus.  Does not include itself.
    '''
    glom_corr = {}
    return_dict = {}
    comp_set = {}
    # mask is a 1d array of 0,1 where data is only used if mask ==1.  Not a true mask, but whatever.
    # Running deepcopy to make sure original vectors are not still linked
    for odor in odr.dp.odors:
        return_dict[odor[1]] = {}
        comp_sets = []
        responses = odr.integrated_vectors[odor[0]]["responses"]
        comp_set["resp" + str(odor[1])] = [responses, responses]
        if "shams" in odr.integrated_vectors[odor[0]].keys():
            for i, value in enumerate(odr.integrated_vectors[odor[0]]["shams"]):
                comp_set["resp" + str(odor[1]) + "_sham" + str(i)] = [responses, value]

        for key, value in comp_set.items():
            v1 = copy.deepcopy(np.asarray(value[0]))   # Need to remove quantities for matmul
            v2 = copy.deepcopy(np.asarray(value[1]))
            print("input v", v1.shape, v2.shape)
            # The means of the datapoint vectors
            vm1 = np.mean(v1, axis=0)
            vm2 = np.mean(v2, axis=0)
            # The difference from the mean for each vector
            v1_diff = (v1 - vm1)
            v2_diff = (v2 - vm2)
            # Expand dims so broadcasting during matmul works correctly.  Otherwise this is the square of v_diff
            denom1 = np.expand_dims(np.sqrt(np.sum(v1_diff * v1_diff, axis=0)), axis=0)
            denom2 = np.expand_dims(np.sqrt(np.sum(v2_diff * v2_diff, axis=0)), axis=0)

            print("denom", denom1.shape, denom2.shape)
            # This is the formation of all the cross multiplied v_diff vectors
            neum = np.matmul(v1_diff.T, v2_diff)
            print("neum", neum.shape)
            # Creates a square denom matrix that will be divided into the neum matrix after taking the sq root.
            denom = np.matmul(denom1.T, denom2)
            print("denom", denom.shape)
            rcoef = neum / denom
            # Subtracting 1 from the diagonal to remove self-self correlations.
            no_self = rcoef - np.eye(rcoef.shape[0], rcoef.shape[1])
            # Finding glomeruli that correlate above criterion and storing in a dict with the correlation coef.
            # +1 is used because the np arrays start at 0 but glomeruli index starts at 1.
            for i in range(no_self.shape[0]):
                temp = no_self[:,i]
                g_xcor = np.where(temp>criterion)[0]
                gc_coef = temp[g_xcor]
                glom_corr[i+1] = [[g_x+1, coef] for g_x, coef in zip(g_xcor, gc_coef)]
            # Removes dictionary elements where there are no significant correlations.
            sig_gloms = {key1: value1 for key1, value1 in glom_corr.items() if len(value) > 0}
            return_dict[odor[1]][key] = {"r_gs": rcoef, "sig_gs": sig_gloms, "sig_level": criterion}
    return return_dict


def corr_vs_ctrl(dp, odor_id, criterion=0.55):
    cmasks, __ = match_stimuli(dp, odor_id, n_masks=1, region="odor_post")
    vin = dp.exp_dict["sdata-leak"]
    active_masks = dp.exp_dict["mask_dict"]["odor_post_mask_by_odor"]
    ctrl_masks = {}
    ctrl_masks[odor_id] = cmasks
    combined_masks = {}
    for key, mask in active_masks.items():
        if key == odor_id:
            if np.min(mask) < -0.1:
                combined_masks["a" + str(key)] = (-1* mask).astype(np.bool)  # Some masks are -1,0
            else:
                combined_masks["a" + str(key)] = mask.astype(np.bool)

    for key, m_set in ctrl_masks.items():
        if key == odor_id:
            for vkey, vmask in m_set.items():
                combined_masks["c" + str(vkey)] = vmask.astype(np.bool)

    corr_dict = twod_vcorr(vin, masks=combined_masks, exclusion_mask=False, criterion=criterion)
    return corr_dict

def odr_corr_vs_ctrl(odr, odor_id, criterion=0.55):
    vin = odr.integrated_vectors[odor_id]["responses"]
    ctrls = odr.integrated_vectors[odor_id]["shams"]
    active_masks = None
    ctrl_masks = vin.shape[0]
    ctrl_masks[odor_id] = cmasks
    combined_masks = {}
    for key, mask in active_masks.items():
        if key == odor_id:
            if np.min(mask) < -0.1:
                combined_masks["a" + str(key)] = (-1* mask).astype(np.bool)  # Some masks are -1,0
            else:
                combined_masks["a" + str(key)] = mask.astype(np.bool)

    for key, m_set in ctrl_masks.items():
        if key == odor_id:
            for vkey, vmask in m_set.items():
                combined_masks["c" + str(vkey)] = vmask.astype(np.bool)

    corr_dict = twod_vcorr(vin, masks=combined_masks, exclusion_mask=False, criterion=criterion)
    return corr_dict


def sig_responders(odr, sig_prob=0.975, criterion=0.35):
    '''
    Takes a two_d array, typically [datapoints, glomeruli] and does cross correlation across the data_point vectors.
    Returns a dictionary with glomeruli with significant correlations with other glomeruli,  Analysis can be restricted
    by providing a [1,0] mask for the regions to restrict the analysis to
    :param vin:  a 2d array of vectors to calculate cross-correlations.
    :param masks:  masking arrays for data set to restrict for the analysis
    :param exclusion_mask:  whether to calculate correlation for data in region outside the masked regions
    :param criterion:  A value the x-correlation needs to exceed to be added to the return_dict
    :return:  All glomeruli that are above criterion for x-correlation with this glomerulus.  Does not include itself.
    '''
    return_dict = {"p_level": sig_prob}
    norm_d = norm(loc=0.0, scale=1.0)  #stabdard normal distribution

    for odor in odr.dp.odors:
        responses = odr.integrated_vectors[odor[0]]["responses"]
        n_shams = len(odr.integrated_vectors[odor[0]]["shams"])
        shape_r = responses.shape
        comp_set = np.zeros([n_shams, shape_r[0], shape_r[1]])
        for i, value in enumerate(odr.integrated_vectors[odor[0]]["shams"]):
            comp_set[i,:,:] = value
        comp_avg = np.average(comp_set, axis=0)
        comp_sd = np.std(comp_set, axis=0)
        resp_z = (responses - comp_avg) / comp_sd
        sig_resp = np.where(resp_z > norm_d.ppf(sig_prob), 1, 0)
        p_sig = np.sum(sig_resp, axis=0) / float(shape_r[0])
        #z_sig = zscore(p_sig)
        #sigs = np.where(np.abs(z_sig) > norm_d.ppf(sig_prob))[0]
        sigs = np.where(p_sig > criterion)[0]
        sig_responders =  sigs + 1
        sig_levels = p_sig[sigs]
        return_dict[odor[0]] = {"responders": sig_responders, "respond_p": sig_levels, "all_probs": p_sig}

    return return_dict


def two_way_anova(odr):
    pass
def one_way_anova(ord):
    pass
