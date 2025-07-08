import pytest

from pathlib import Path
import numpy as np
import numpy.testing as npt
from xbout.tests.test_load import bout_xyt_example_files
import xarray.testing as xrt

from xstorm import open_stormdataset
from xstorm.load import apply_unnormalise


def create_BOUTinp(path, *, realistic_geometry="doublenull"):
    inputpath = str(Path(path).parent) + "/BOUT.inp"

    with open(inputpath, "w") as f:
        f.write("[filaments]\n")
        f.write(f"realistic_geometry = {realistic_geometry}\n")
        f.write("Z = 1\n")
        f.write("B_0 = 1\n")
        f.write("m_i = 1\n")
        f.write("T_i0 = 1\n")
        f.write("T_e0 = 1\n")
        f.write("n_0 = 1\n")
        f.write("Lx = 1\n")
        f.write("Ly = 1\n")
        f.write("Lz = 1\n")

        # create empty [mesh] section to avoid error when trying to read from it
        f.write("[mesh]")

    return inputpath


class TestSave:
    @pytest.mark.parametrize("realistic_geometry", ["doublenull", "slab"])
    def test_reload_all(
        self, tmp_path_factory, bout_xyt_example_files, realistic_geometry
    ):
        # Create data
        path = bout_xyt_example_files(
            tmp_path_factory, nxpe=4, nype=5, nt=1, grid="grid", write_to_disk=True
        )

        gridpath = str(Path(path).parent) + "/grid.nc"

        inputpath = create_BOUTinp(path, realistic_geometry=realistic_geometry)

        # Load it as a stormdataset
        original = open_stormdataset(
            datapath=path,
            gridfilepath=gridpath,
            inputfilepath=inputpath,
        )

        # Save it to a netCDF file
        savepath = str(Path(path).parent) + "/temp_boutdata.nc"
        original.storm.save(savepath=savepath)

        # recreate regions after loading to avoid machine-precision rounding errors in
        # dask computations
        from xbout.region import _create_regions_toroidal

        original.load()
        del original.attrs["regions"]
        original = _create_regions_toroidal(original)

        # Load it again
        recovered = open_stormdataset(savepath, inputfilepath=inputpath)

        # recreate regions after loading to avoid machine-precision rounding errors in
        # dask computations
        recovered.load()
        del recovered.attrs["regions"]
        recovered = _create_regions_toroidal(recovered)

        # Compare
        xrt.assert_identical(original.load(), recovered.load())

    def test_reload_separate_variables(self, tmp_path_factory, bout_xyt_example_files):

        path = bout_xyt_example_files(
            tmp_path_factory, nxpe=4, nype=1, nt=1, grid="grid", write_to_disk=True
        )

        gridpath = str(Path(path).parent) + "/grid.nc"

        inputpath = create_BOUTinp(path)

        # Load it as a stormdataset
        original = open_stormdataset(
            datapath=path,
            gridfilepath=gridpath,
            inputfilepath=inputpath,
        )

        # Save it to a netCDF file
        savepath = str(Path(path).parent) + "/temp_boutdata.nc"
        original.storm.save(savepath=savepath, separate_vars=True)

        # Load it again
        savepath = str(Path(path).parent) + "/temp_boutdata_*.nc"
        recovered = open_stormdataset(savepath, inputfilepath=inputpath)

        # Compare
        xrt.assert_identical(recovered, original)

    @pytest.mark.parametrize("realistic_geometry", ["doublenull", "slab"])
    def test_reapply_unnormalise(
        self, tmp_path_factory, bout_xyt_example_files, realistic_geometry
    ):
        # Create data
        path = bout_xyt_example_files(
            tmp_path_factory, nxpe=4, nype=5, nt=1, grid="grid", write_to_disk=True
        )

        gridpath = str(Path(path).parent) + "/grid.nc"

        inputpath = create_BOUTinp(path, realistic_geometry=realistic_geometry)

        # Load it as a stormdataset
        ds = open_stormdataset(
            datapath=path,
            gridfilepath=gridpath,
            inputfilepath=inputpath,
        )

        with pytest.warns(UserWarning):
            ds2 = apply_unnormalise(ds)

        xrt.assert_identical(ds, ds2)


class TestMethods:
    def test_integrate_midpoints(self, tmp_path_factory, bout_xyt_example_files):
        # Create data
        dataset_list = bout_xyt_example_files(
            None, lengths=(4, 100, 110, 120), nxpe=1, nype=1, nt=1, syn_data_type=1
        )

        temp_dir = tmp_path_factory.mktemp("filaments")

        inputpath = create_BOUTinp(temp_dir, realistic_geometry="slab")

        ds = open_stormdataset(dataset_list, inputfilepath=inputpath)

        t = np.linspace(0.0, 8.0, 4)[:, np.newaxis, np.newaxis, np.newaxis]
        x = np.linspace(0.05, 9.95, 100)[np.newaxis, :, np.newaxis, np.newaxis]
        y = np.linspace(0.1, 21.9, 110)[np.newaxis, np.newaxis, :, np.newaxis]
        z = np.linspace(0.15, 35.85, 120)[np.newaxis, np.newaxis, np.newaxis, :]
        ds["time"].data[...] = t.squeeze()
        ds["dx"].data[...] = 0.1
        ds["dy"].data[...] = 0.2
        ds["dz"].data[...] = 0.3

        tfunc = 1.5 * t
        xfunc = x**2
        yfunc = 10.0 * y - y**2
        zfunc = 2.0 * z**2 - 30.0 * z
        ds["n"].data[...] = tfunc * xfunc * yfunc * zfunc

        tintegral = 48.0
        xintegral = 1000.0 / 3.0
        yintegral = 5.0 * 22.0**2 - 22.0**3 / 3.0
        zintegral = 2.0 * 36.0**3 / 3.0 - 15.0 * 36.0**2
        npt.assert_allclose(
            ds.storm.integrate_midpoints("n", dims="time"),
            (tintegral * xfunc * yfunc * zfunc).squeeze(),
        )
        npt.assert_allclose(
            ds.storm.integrate_midpoints("n", dims="radial"),
            (tfunc * xintegral * yfunc * zfunc).squeeze(),
            rtol=1.0e-4,
        )
        npt.assert_allclose(
            ds.storm.integrate_midpoints("n", dims="parallel"),
            (tfunc * xfunc * yintegral * zfunc).squeeze(),
            rtol=1.0e-4,
        )
        npt.assert_allclose(
            ds.storm.integrate_midpoints("n", dims="binormal"),
            (tfunc * xfunc * yfunc * zintegral).squeeze(),
            rtol=1.0e-4,
        )
        npt.assert_allclose(
            ds.storm.integrate_midpoints("n", dims=["time", "radial"]),
            (tintegral * xintegral * yfunc * zfunc).squeeze(),
            rtol=1.0e-4,
        )
        npt.assert_allclose(
            ds.storm.integrate_midpoints("n", dims=["time", "parallel"]),
            (tintegral * xfunc * yintegral * zfunc).squeeze(),
            rtol=1.0e-4,
        )
        npt.assert_allclose(
            ds.storm.integrate_midpoints("n", dims=["time", "binormal"]),
            (tintegral * xfunc * yfunc * zintegral).squeeze(),
            rtol=1.0e-4,
        )
        npt.assert_allclose(
            ds.storm.integrate_midpoints("n", dims=["radial", "parallel"]),
            (tfunc * xintegral * yintegral * zfunc).squeeze(),
            rtol=1.0e-4,
        )
        npt.assert_allclose(
            ds.storm.integrate_midpoints("n", dims=["radial", "binormal"]),
            (tfunc * xintegral * yfunc * zintegral).squeeze(),
            rtol=1.0e-4,
        )
        npt.assert_allclose(
            ds.storm.integrate_midpoints("n", dims=["parallel", "binormal"]),
            (tfunc * xfunc * yintegral * zintegral).squeeze(),
            rtol=1.2e-4,
        )
        npt.assert_allclose(
            ds.storm.integrate_midpoints("n", dims=["time", "radial", "parallel"]),
            (tintegral * xintegral * yintegral * zfunc).squeeze(),
            rtol=1.0e-4,
        )
        npt.assert_allclose(
            ds.storm.integrate_midpoints("n", dims=["time", "radial", "binormal"]),
            (tintegral * xintegral * yfunc * zintegral).squeeze(),
            rtol=1.0e-4,
        )
        npt.assert_allclose(
            ds.storm.integrate_midpoints("n", dims=["time", "parallel", "binormal"]),
            (tintegral * xfunc * yintegral * zintegral).squeeze(),
            rtol=1.2e-4,
        )
        # default dims
        npt.assert_allclose(
            ds.storm.integrate_midpoints("n"),
            (tfunc * xintegral * yintegral * zintegral).squeeze(),
            rtol=1.4e-4,
        )
        npt.assert_allclose(
            ds.storm.integrate_midpoints(
                "n", dims=["time", "radial", "parallel", "binormal"]
            ),
            (tintegral * xintegral * yintegral * zintegral),
            rtol=1.4e-4,
        )
