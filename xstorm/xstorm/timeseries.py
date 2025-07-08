import numpy as np
import numpy.ma as ma
from scipy.ndimage import label

import xarray as xr
from xarray import apply_ufunc

from .utils import _update_name, _update_units, get_dim_from_coord


def conditional_average(da, coord=None, level=2.5, delta=None):
    """
    Conditional average along a single dimension.

    Assumes input DataArray has already had the mean subtracted
    and been divided by its standard deviation.

    Algorithm first finds all contiguous regions where signal
    exceeds the value specified by `level`. Then it finds the
    position of the maximum value of each region. Then, working
    from the start of the array, it drops any maxima which are not
    separated from the last maxima by more than `delta` along `coord`.
    Then samples a window of width `delta` either side of each
    of the remaining peaks, and returns the average of those windows.

    Parameters
    ----------
    da : xarray.DataArray
    coord : str
    level : float

    delta : float

    Returns
    -------
    cond_avg : xarray.DataArray
        Averaged signal, of length 2*delta along coord.
    """

    if coord not in da.coords:
        raise ValueError
    dim = get_dim_from_coord(coord, da)

    # Check for non-constant sampling rates
    diffs = da[coord].diff(dim).values
    if not np.allclose(diffs, diffs[0]):
        raise ValueError(f"Input data has non-constant sampling rate along" "{dim}")

    step = (da[coord][1] - da[coord][0]).values
    window = int(np.floor(delta / step))

    result = apply_ufunc(
        _cond_avg_ufunc,
        da,
        input_core_dims=[[dim]],
        output_core_dims=[["window_dim"]],
        kwargs={"level": level, "window": window},
        dask="parallelized",
        output_sizes={"window_dim": 2 * window + 1},
        output_dtypes=[np.float64],
        keep_attrs=True,
    )

    # Construct new coord along windowed axis with name of original dimension
    result = result.rename(window_dim=dim)
    new_coord_data = np.arange(
        start=-window * step, stop=(window + 1) * step, step=step
    )
    result[coord] = xr.DataArray(new_coord_data, dims=[dim], name=coord)

    return result


def _cond_avg_1d_ufunc(arr, level=2.5, window=50):
    """
    Conditional average along 1D data.

    numpy-style ufunc so expects a numpy-style array where the final
    axis is the axis to be conditionally-averaged along.

    Returns
    -------
    cond_avg : array
        numpy-style array of length 2*window+1 along the final axis.
    """

    peak_inds = _peak_finder_ufunc(arr, level, window)

    windows = [arr[slice(ind - window, ind + window + 1)] for ind in peak_inds]
    a = np.stack(windows, axis=-1)
    cond_avg = np.mean(a, axis=-1)
    return cond_avg


def _cond_avg_ufunc(arr, level=2.5, window=50):
    return np.apply_along_axis(
        _cond_avg_1d_ufunc, axis=-1, arr=arr, level=level, window=window
    )


def _peak_finder_ufunc(arr, level=2.5, window=50):

    # Find  contiguous regions above level
    blob_labels, n_blobs = label(arr > level)

    peak_inds = []
    previous_indmax = 0
    for i in range(1, n_blobs + 1):
        blob = ma.masked_where(blob_labels != i, arr)
        indmax = np.argmax(blob, axis=-1)

        # Only keep maxima which are separated by more than window num of points
        if (previous_indmax + window) < indmax < (len(arr) - window):
            peak_inds.append(indmax)
            previous_indmax = indmax

    # Report number of peaks found
    print(f"Found {len(peak_inds)} peaks")

    return np.asarray(peak_inds)


def waiting_times(da, coord="time", level=2.5, delta=None):
    """
    Computes the waiting times between level-crossing events.

    Parameters
    ----------
    da
    coord
    level
    delta

    Returns
    -------

    """

    if coord not in da.coords:
        raise ValueError
    dim = get_dim_from_coord(coord, da)

    # TODO check for non-constant sampling rates?
    step = (da[coord][1] - da[coord][0]).values
    window = np.fix(delta / step)

    # TODO would this fail on 1D data that returns ragged result?
    raise NotImplementedError
    crossings = apply_ufunc(
        _level_crossings_ufunc,
        da,
        input_core_dims=[dim],
        kwargs={"level": level, "window": window},
        dask="parallelized",
        output_dtypes=[np.int],
        keep_attrs=True,
    )

    times = da[coord].isel(dim=crossings)
    waiting_times = times.diff(dim)

    return _update_name(waiting_times, suffix="_waiting_times")
