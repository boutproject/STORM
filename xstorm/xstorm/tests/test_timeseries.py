import pytest
import numpy as np

import xarray as xr
import xarray.testing as xrt

from xstorm.timeseries import (
    _peak_finder_ufunc,
    _cond_avg_1d_ufunc,
    _cond_avg_ufunc,
    conditional_average,
)
from xstorm.accessors import StormDataArrayAccessor


class TestPeakFinder:
    def test_find_peaks_1d_ufunc(self):
        arr = np.random.rand(20)
        print(arr.shape)
        arr[5] = 1.7
        arr[14] = 1.8
        arr[15] = 1.9

        peaks = _peak_finder_ufunc(arr, level=1.5, window=3)
        expected = np.array([5, 15])
        np.testing.assert_array_equal(peaks, expected)

    @pytest.mark.skip
    def test_find_peaks_2d_ufunc(self):
        arr = np.random.rand(20, 2)
        arr[5, 0] = 1.7
        arr[14, 0] = 1.8
        arr[15, 0] = 1.9
        arr[5, 0] = 1.7
        arr[14, 0] = 1.9
        arr[15, 0] = 1.8

        crossings = _peak_finder_ufunc(arr, level=1.5, window=3)
        expected = np.array([[5, 15], [7, 12]])
        np.testing.assert_array_equal(peaks, expected)


class TestConditionalAverage:
    def test_1d_cond_avg_ufunc(self):
        arr = np.random.rand(20)
        arr[5] = 1.7
        arr[14] = 1.8
        arr[15] = 1.9

        cond_avg = _cond_avg_1d_ufunc(arr, level=1.5, window=2)

        window1 = np.take(arr, indices=[3, 4, 5, 6, 7])
        window2 = np.take(arr, indices=[13, 14, 15, 16, 17])
        expected = 0.5 * (window1 + window2)
        np.testing.assert_array_equal(cond_avg, expected)

    def test_2d_cond_avg_ufunc(self):
        arr = np.random.rand(2, 20)
        arr[0, 5] = 1.7
        arr[0, 14] = 1.8
        arr[0, 15] = 1.9
        arr[1, 10] = 1.7
        arr[1, 17] = 1.8

        cond_avg = _cond_avg_ufunc(arr, level=1.5, window=2)

        window1 = np.take(arr[0, :], indices=[3, 4, 5, 6, 7])
        window2 = np.take(arr[0, :], indices=[13, 14, 15, 16, 17])
        window3 = np.take(arr[1, :], indices=[8, 9, 10, 11, 12])
        window4 = np.take(arr[1, :], indices=[15, 16, 17, 18, 19])
        avg1 = 0.5 * (window1 + window2)
        avg2 = 0.5 * (window3 + window4)
        expected = np.stack([avg1, avg2], axis=0)
        np.testing.assert_array_equal(cond_avg, expected)

    def test_1d_cond_avg(self):
        arr = np.random.rand(21)
        arr[5] = 1.7
        arr[14] = 1.8
        arr[15] = 1.9
        times = np.linspace(5, 10, num=21)
        da = xr.DataArray(data=arr, dims=["t"], coords=[times])

        cond_avg = conditional_average(da, coord="t", level=1.5, delta=0.6)

        window1 = np.take(arr, indices=[3, 4, 5, 6, 7])
        window2 = np.take(arr, indices=[13, 14, 15, 16, 17])
        expected_data = 0.5 * (window1 + window2)
        expected_coord = xr.DataArray([-0.5, -0.25, 0.0, 0.25, 0.5], dims=["t"])
        expected = xr.DataArray(expected_data, dims=["t"], coords={"t": expected_coord})
        xrt.assert_equal(cond_avg, expected)

    def test_2d_cond_avg(self):
        ...

    def test_uneven_sampling(self):
        arr = np.random.rand(20)
        times = np.linspace(5, 10, num=20)
        times[7] = 32.1
        da = xr.DataArray(data=arr, dims=["t"], coords=[times])

        with pytest.raises(ValueError):
            conditional_average(da, coord="t", level=1.5)


class TestWaitingTimes:
    ...


def power(signal, dim):
    """Compute total power in signal"""

    if dim in signal.coords:
        T = signal.coords[dim].max() - signal.coords[dim].min()
    else:
        T = len(signal.sizes[dim])

    return (signal**2).integrate(coord=dim) / T


class TestPowerSpectralDensity:
    def test_integrated_power(self):
        # Create signal
        size = 100
        rng = np.random.default_rng(seed=54328)
        data = rng.random(size=size)
        time = xr.DataArray(np.linspace(0.2, 2.2, size), dims=["time"])
        signal = xr.DataArray(data=data, dims=["time"], coords={"time": time})

        # Compute and integrate PSD
        psd = signal.storm.power_spectra(dim="time")
        integrated_power = psd.integrate(coord="freq_time")

        # Compute total power in signal directly
        actual = power(signal, dim="time")

        # Compare
        xrt.assert_allclose(integrated_power, actual, atol=1e-3)
