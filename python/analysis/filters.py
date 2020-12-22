import numpy as np
# import math
from scipy.signal import butter, filtfilt, sosfiltfilt   # savgol_filter
import quantities as pq

def smooth(data, corner_f=1.0, dt=0.078, axis=0, silence=False):
    '''
    Performs a 5 pole butterworth filter on data
    :param data: data set to analyze
    :param corner_f: corner frequency to filter in Hz
    :param dt: average interframe time interval in sec
    :param axis:  which direction the data run
    :return:  filtered data
    '''
    add_unit = False
    transpose_it = False
    if axis == 1:
        transpose_it = True
        data = data.T
        axis = 0
    if corner_f == "raw":
        return data
    if isinstance(data, (pq.Quantity)):
        return_unit = data.units
        add_unit = True
    buffer = int(2. / (corner_f * dt))  # points to add to beginning and end of data to avoid end effects
    if len(data.shape) == 1:
        end_buff = data[-buffer:]
        start_buff = data[:buffer]
    else:
        end_buff = data[-buffer:, :]
        start_buff = data[:buffer, :]
    data = np.append(start_buff, data, axis=0)
    data = np.append(data, end_buff, axis=0)
    if not silence: print("endpad_filter ", buffer, data.shape)
    sample_freq = 1.0 / dt
    #norm_corner = corner_f /(sample_freq / 2.0)
    sos = butter(4, float(corner_f), btype="lowpass", output="sos", fs=float(sample_freq))
    sdata = sosfiltfilt(sos, data, axis=axis)
    if len(sdata.shape) == 1:
        sdata = sdata[buffer:-buffer]
    else:
        sdata = sdata[buffer:-buffer,:]
    if transpose_it:
        sdata = sdata.T
    if add_unit:
        sdata = pq.Quantity(sdata, return_unit)
    return sdata

def slider(array, window=500, axis=0, **kwargs):
    '''
    Slides a window across the data putting the median of the selected points into the value at the window center.
    :param array:  2 d array of data to filter in 1 direction
    :param window: filtering window
    :param axis: which axis to filter
    :param kwargs: any parameters or overrides
    :return: lower filtered array which should correspond to the background.
    '''
    if kwargs is None:
        kwargs = {}
    # Checks for a filter function override in kwargs
    if "filter" in kwargs.keys():
        filter = kwargs.pop("filter")
    else:
        filter = np.median
    test_array = np.asarray([i for i in range(100)], dtype=np.float32)
    if not silence: print("testing filter", filter(test_array))
    # Checks for an axis override in kwargs and transposes
    if axis == 1 and "axis" not in kwargs.keys():
        array = array.T
        kwargs["axis"] = 0
    elif "axis" in kwargs.keys() and kwargs["axis"] == 1:
        array = array.T
        kwargs["axis"] = 0
        axis = 1
    else:
        kwargs["axis"] = 0
    ashape = array.shape
    range0 = ashape[0]
    range1 = ashape[1]
    midpoint = window // 2
    # cycling is set to provide data for arry elements
    cycles = range0
    #clean up N-term of array which often has an artifact but not significant data
    array[0:20,:] = array[20:40,:]
    # Buffer ends of array to remove edge effects.
    farray = np.zeros([range0 + window, range1])
    if isinstance(array, (pq.Quantity)):
        farray = pq.Quantity(farray, array.units)
    out_array = np.zeros_like(array)
    farray[0:midpoint,:] = array[0:midpoint]
    farray[midpoint:midpoint+range0, :] = array
    farray[midpoint+range0:,:] = array[-midpoint:,:]
    #print("kwargs", kwargs, cycles, midpoint, range0)
    for i in range(cycles):
        # Take a window sized slice of the data
        minset = farray[i:i+window, :]
        # Place median value into output array
        out_array[i,:] = filter(minset, **kwargs)
        #if not i % 500:
        #    print("point", i, out_array[i,4])
    # If the array had been transposed, return it to its original shape
    if axis == 1:
        out_array = out_array.T

    return out_array

def high_pass(data, corner_f=1.0, dt=0.078, axis=0, silence=False):
    '''
    Performs a 5 pole butterworth filter on data
    :param data: data set to analyze
    :param corner_f: corner frequency to filter in Hz
    :param dt: average interframe time interval in sec
    :param axis:  which direction the data run
    :return:  filtered data
    '''
    add_unit = False
    transpose_it = False
    if axis == 1:
        transpose_it = True
        data = data.T
        axis = 0
    if corner_f == "raw":
        return data
    if isinstance(data, (pq.Quantity)):
        return_unit = data.units
        add_unit = True
    buffer = int(2. / (corner_f * dt))  # points to add to beginning and end of data to avoid end effects
    if len(data.shape) == 1:
        end_buff = data[-buffer:]
        start_buff = data[:buffer]
    else:
        end_buff = data[-buffer:, :]
        start_buff = data[:buffer, :]
    data = np.append(start_buff, data, axis=0)
    data = np.append(data, end_buff, axis=0)
    if not silence: print("endpad_filter ", buffer, data.shape)
    sample_freq = 1.0 / dt
    #norm_corner = corner_f /(sample_freq / 2.0)
    sos = butter(4, float(corner_f), btype="highpass", output="sos", fs=float(sample_freq))
    sdata = sosfiltfilt(sos, data, axis=axis)
    if len(sdata.shape) == 1:
        sdata = sdata[buffer:-buffer]
    else:
        sdata = sdata[buffer:-buffer,:]
    if transpose_it:
        sdata = sdata.T
    if add_unit:
        sdata = pq.Quantity(sdata, return_unit)
    return sdata

def band_pass(data, corner_f=[50.0,70.0], dt=0.078, axis=0, silence=False):
    '''
    Performs a 5 pole butterworth filter on data
    :param data: data set to analyze
    :param corner_f: corner frequency to filter in Hz
    :param dt: average interframe time interval in sec
    :param axis:  which direction the data run
    :return:  filtered data
    '''
    add_unit = False
    transpose_it = False
    if axis == 1:
        transpose_it = True
        data = data.T
        axis = 0
    if corner_f == "raw":
        return data
    if isinstance(data, (pq.Quantity)):
        return_unit = data.units
        add_unit = True
    buffer = int(2. / (corner_f[0] * dt))  # points to add to beginning and end of data to avoid end effects
    if len(data.shape) == 1:
        end_buff = data[-buffer:]
        start_buff = data[:buffer]
    else:
        end_buff = data[-buffer:, :]
        start_buff = data[:buffer, :]
    data = np.append(start_buff, data, axis=0)
    data = np.append(data, end_buff, axis=0)
    if not silence: print("endpad_filter ", buffer, data.shape)
    sample_freq = 1.0 / dt
    #norm_corner = corner_f /(sample_freq / 2.0)
    sos = butter(4, Wn=(float(corner_f[0]), float(corner_f[1])), btype="bandpass", output="sos", fs=float(sample_freq))
    sdata = sosfiltfilt(sos, data, axis=axis)
    if len(sdata.shape) == 1:
        sdata = sdata[buffer:-buffer]
    else:
        sdata = sdata[buffer:-buffer,:]
    if transpose_it:
        sdata = sdata.T
    if add_unit:
        sdata = pq.Quantity(sdata, return_unit)
    return sdata

def band_stop(data, corner_f=[50.0,70.0], dt=0.078, axis=0, silence=False):
    '''
    Performs a 5 pole butterworth filter on data
    :param data: data set to analyze
    :param corner_f: corner frequency to filter in Hz
    :param dt: average interframe time interval in sec
    :param axis:  which direction the data run
    :return:  filtered data
    '''
    add_unit = False
    transpose_it = False
    if axis == 1:
        transpose_it = True
        data = data.T
        axis = 0
    if corner_f == "raw":
        return data
    if isinstance(data, (pq.Quantity)):
        return_unit = data.units
        add_unit = True
    buffer = int(2. / (corner_f[0] * dt))  # points to add to beginning and end of data to avoid end effects
    if len(data.shape) == 1:
        end_buff = data[-buffer:]
        start_buff = data[:buffer]
    else:
        end_buff = data[-buffer:, :]
        start_buff = data[:buffer, :]
    data = np.append(start_buff, data, axis=0)
    data = np.append(data, end_buff, axis=0)
    if not silence: print("endpad_filter ", buffer, data.shape)
    sample_freq = 1.0 / dt
    #norm_corner = corner_f /(sample_freq / 2.0)
    sos = butter(4, Wn=(float(corner_f[0]), float(corner_f[1])), btype="bandstop", output="sos", fs=float(sample_freq))
    sdata = sosfiltfilt(sos, data, axis=axis)
    if len(sdata.shape) == 1:
        sdata = sdata[buffer:-buffer]
    else:
        sdata = sdata[buffer:-buffer,:]
    if transpose_it:
        sdata = sdata.T
    if add_unit:
        sdata = pq.Quantity(sdata, return_unit)
    return sdata


def data_reduce(t, y, by=11, y_axis=0, y_avg=True):
    '''
    Assumes that t is one dimensional and y is one or two dimensional and that y_axis dimension
    is same length as t
    :param t:  time series
    :param y:  y data can be 1 or 2 dimensional
    :param by: reduction factor
    :param y_axis:  axis matching time series
    :param y_ave:  whether to decimate y or average y for each point
    :return: time series and y data reduced by the by factor
    '''
    # n_extra tracks whether by divides evenly into len(t).  If not add one more data point
    n_extra = 0
    if y_axis == 1:
        y = y.T
    len_t = t.shape[0]
    #print(len_t, type(len_t))
    # Short is the remainder to indicate that one extra data point needs to be added to get everything
    short = len_t % by
    l_ax = len_t // by
    # pad- How many points need to be added to the end to get match to by with no remainder
    pad = by - short
    y_dim = len(y.shape)
    # if there is padding needed, then pad < by
    if pad < by:
        n_extra = 1
        t = np.append(t, t[-pad:], axis=0)
        if y_dim == 1:
            y = np.append(y, y[-pad:], axis=0)
        else:
            y = np.append(y, y[-pad:,:], axis=0)
    # Condition that average of the removed points is used as a replacement value
    if y_avg:
        # To make the average center abound the replacement point add half the by number to the front of the array
        # And the rest to the end of the array
        y_pad = by // 2
        if y_dim == 1:
            y = np.append(y[:y_pad], y, axis=0)
            y = np.append(y, y[-(by-y_pad):], axis=0)
        else:
            y = np.append(y[:y_pad,:], y, axis=0)
            y = np.append(y, y[-(by-y_pad):,:], axis=0)
        # ax0 should now be an even division of the shape
        ax0 = y.shape[0] // by
        #print("shapes y,ax, by", y.shape, ax0, by, pad, y_pad, l_ax)
        # Slightly different syntax if 1d or 2d
        if y_dim == 1:
            y_2d = np.reshape(y, [ax0, by])
            y_1d = np.average(y_2d, axis=1)
            y_red = y_1d[:-1]
        else:
            y_3d = np.reshape(y, [ax0, by, y.shape[1]])
            y_2d = np.average(y_3d, axis=1)
            y_red = y_2d[:-1, :]

    else:
        if y_dim == 1:
            y_3d = np.reshape(y, [l_ax+n_extra, by])
            y_red = y_3d[:,0]
        else:
            y_3d = np.reshape(y, [l_ax+n_extra, by, y.shape[1]])
            y_red = y_3d[:,0,:]

    t_2d = np.reshape(t, [l_ax+n_extra, by])
    #print("tshapes", t.shape, t_2d.shape, l_ax, n_extra, by)
    t_red = t_2d[:,0]

    if y_axis == 1:
        y_red = y_red.T

    return t_red, y_red
