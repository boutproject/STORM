import pytest

from xbout.tests.test_load import bout_xyt_example_files
import numpy as np
import numpy.testing as npt
import xarray.testing as xrt

from xstorm import open_stormdataset

from .test_stormdataset import create_BOUTinp


class TestMethods:
    def test_spatial_dims(self, tmpdir_factory, bout_xyt_example_files):
        # Create data
        dataset_list = bout_xyt_example_files(
            None, lengths=(2, 3, 5, 2), nxpe=1, nype=1, nt=1
        )

        temp_dir = tmpdir_factory.mktemp("filaments")

        inputpath = create_BOUTinp(temp_dir, realistic_geometry="slab")

        ds = open_stormdataset(dataset_list, inputfilepath=inputpath)

        assert set(ds.filament._spatial_dims) == set(["radial", "parallel", "binormal"])

    def test_get_background_subtracted_fluctuation(
        self, tmpdir_factory, bout_xyt_example_files
    ):
        # Create data
        dataset_list = bout_xyt_example_files(
            None, lengths=(2, 3, 5, 2), nxpe=1, nype=1, nt=1, syn_data_type=1
        )

        temp_dir = tmpdir_factory.mktemp("filaments")

        inputpath = create_BOUTinp(temp_dir, realistic_geometry="slab")

        ds = open_stormdataset(dataset_list, inputfilepath=inputpath)

        ds["n_eq"] = ds["n"].isel(time=0) / 2.0

        deltan, total_mass = ds.filament._get_background_subtracted_fluctuation(
            "n", return_mass=True
        )

        npt.assert_allclose(deltan, 0.5e6 * np.ones([2, 3, 5, 2]))

        npt.assert_allclose(total_mass, 30.0 * 0.5e6 * np.ones(2))

    @pytest.mark.parametrize(
        "dims", [None, pytest.param(["radial", "binormal"], marks=pytest.mark.long)]
    )
    def test_CoM(self, tmpdir_factory, bout_xyt_example_files, dims):
        # Create data
        dataset_list = bout_xyt_example_files(
            None, lengths=(2, 10, 11, 12), nxpe=1, nype=1, nt=1, syn_data_type=1
        )

        temp_dir = tmpdir_factory.mktemp("filaments")

        inputpath = create_BOUTinp(temp_dir, realistic_geometry="slab")

        ds = open_stormdataset(dataset_list, inputfilepath=inputpath)

        ds["n_eq"] = ds["n"].isel(time=0).copy(data=np.ones([10, 11, 12]))

        x0 = np.array([0.2, 0.3])[:, np.newaxis, np.newaxis, np.newaxis]
        y0 = np.array([0.4, 0.5])[:, np.newaxis, np.newaxis, np.newaxis]
        z0 = np.array([0.6, 0.7])[:, np.newaxis, np.newaxis, np.newaxis]
        x = np.linspace(0.0, 1.0, 10)[np.newaxis, :, np.newaxis, np.newaxis] - x0
        y = np.linspace(0.0, 1.0, 11)[np.newaxis, np.newaxis, :, np.newaxis] - y0
        z = np.linspace(0.0, 1.0, 12)[np.newaxis, np.newaxis, np.newaxis, :] - z0
        ds["radial"].data[:] = np.linspace(0.0, 1.0, 10)
        ds["parallel"].data[:] = np.linspace(0.0, 1.0, 11)
        ds["binormal"].data[:] = np.linspace(0.0, 1.0, 12)
        ds["n"] = ds["n"].copy(data=(1.0 + np.exp(-73.0 * (x**2 + y**2 + z**2))))

        CoM = ds.filament.CoM("n", dims=dims)
        if dims is None:
            npt.assert_allclose(CoM["radial"], x0.squeeze(), rtol=1.0e-3)
            npt.assert_allclose(CoM["parallel"], y0.squeeze(), rtol=1.0e-3)
            npt.assert_allclose(CoM["binormal"], z0.squeeze(), rtol=1.0e-3)
        elif dims == ["radial", "binormal"]:
            npt.assert_allclose(
                CoM["radial"], np.repeat(x0[:, 0, :, 0], 11, axis=1), rtol=1.0e-3
            )
            npt.assert_allclose(
                CoM["parallel"].transpose("time", "parallel"),
                np.repeat(np.linspace(0.0, 1.0, 11)[np.newaxis, :], 2, axis=0),
                rtol=1.0e-3,
            )
            npt.assert_allclose(
                CoM["binormal"], np.repeat(z0[:, 0, :, 0], 11, axis=1), rtol=1.0e-3
            )
        else:
            raise ValueError(f"Unrecognised value dims={dims}")

    @pytest.mark.parametrize(
        "dims", [None, pytest.param(["radial", "binormal"], marks=pytest.mark.long)]
    )
    def test_CoM_velocity_finitedifference(
        self, tmpdir_factory, bout_xyt_example_files, dims
    ):
        # Create data
        dataset_list = bout_xyt_example_files(
            None, lengths=(5, 10, 11, 12), nxpe=1, nype=1, nt=1, syn_data_type=1
        )

        temp_dir = tmpdir_factory.mktemp("filaments")

        inputpath = create_BOUTinp(temp_dir, realistic_geometry="slab")

        ds = open_stormdataset(dataset_list, inputfilepath=inputpath)

        ds["n_eq"] = ds["n"].isel(time=0).copy(data=np.ones([10, 11, 12]))

        x0 = np.linspace(0.2, 0.3, 5)[:, np.newaxis, np.newaxis, np.newaxis]
        y0 = np.linspace(0.5, 0.4, 5)[:, np.newaxis, np.newaxis, np.newaxis]
        z0 = np.linspace(0.5, 0.7, 5)[:, np.newaxis, np.newaxis, np.newaxis]
        x = np.linspace(0.0, 1.0, 10)[np.newaxis, :, np.newaxis, np.newaxis] - x0
        y = np.linspace(0.0, 1.0, 11)[np.newaxis, np.newaxis, :, np.newaxis] - y0
        z = np.linspace(0.0, 1.0, 12)[np.newaxis, np.newaxis, np.newaxis, :] - z0
        ds["time"].data[:] = np.linspace(0.0, 2.0, 5)
        ds["radial"].data[:] = np.linspace(0.0, 1.0, 10)
        ds["parallel"].data[:] = np.linspace(0.0, 1.0, 11)
        ds["binormal"].data[:] = np.linspace(0.0, 1.0, 12)
        ds["n"] = ds["n"].copy(data=(1.0 + np.exp(-73.0 * (x**2 + y**2 + z**2))))

        vr = ds.filament.CoM_velocity(
            "n", "radial", method="finite-difference", dims=dims
        )
        vp = ds.filament.CoM_velocity(
            "n", "parallel", method="finite-difference", dims=dims
        )
        vb = ds.filament.CoM_velocity(
            "n", "binormal", method="finite-difference", dims=dims
        )

        if dims is None:
            npt.assert_allclose(vr, 0.05 * np.ones(5), rtol=6.0e-3)
            npt.assert_allclose(vp, -0.05 * np.ones(5), rtol=6.0e-3)
            npt.assert_allclose(vb, 0.1 * np.ones(5), rtol=6.0e-3)
        elif dims == ["radial", "binormal"]:
            npt.assert_allclose(vr, 0.05 * np.ones([5, 11]), rtol=6.0e-3)
            npt.assert_allclose(vb, 0.1 * np.ones([5, 11]), rtol=6.0e-3)
        else:
            raise ValueError(f"Unrecognised value dims={dims}")

    @pytest.mark.parametrize(
        "dims", [None, pytest.param(["radial", "binormal"], marks=pytest.mark.long)]
    )
    def test_CoM_velocity_ExB(self, tmpdir_factory, bout_xyt_example_files, dims):
        # Create data
        dataset_list = bout_xyt_example_files(
            None, lengths=(5, 10, 11, 12), nxpe=1, nype=1, nt=1, syn_data_type=1
        )

        temp_dir = tmpdir_factory.mktemp("filaments")

        inputpath = create_BOUTinp(temp_dir, realistic_geometry="slab")

        ds = open_stormdataset(dataset_list, inputfilepath=inputpath)

        ds["n_eq"] = ds["n"].isel(time=0).copy(data=np.ones([10, 11, 12]))

        x0 = np.linspace(0.2, 0.3, 5)[:, np.newaxis, np.newaxis, np.newaxis]
        y0 = np.linspace(0.5, 0.4, 5)[:, np.newaxis, np.newaxis, np.newaxis]
        z0 = np.linspace(0.5, 0.7, 5)[:, np.newaxis, np.newaxis, np.newaxis]
        x = np.linspace(0.0, 1.0, 10)[np.newaxis, :, np.newaxis, np.newaxis] - x0
        y = np.linspace(0.0, 1.0, 11)[np.newaxis, np.newaxis, :, np.newaxis] - y0
        z = np.linspace(0.0, 1.0, 12)[np.newaxis, np.newaxis, np.newaxis, :] - z0
        ds["time"].data[:] = np.linspace(0.0, 2.0, 5)
        ds["radial"].data[:] = np.linspace(0.0, 1.0, 10)
        ds["parallel"].data[:] = np.linspace(0.0, 1.0, 11)
        ds["binormal"].data[:] = np.linspace(0.0, 1.0, 12)
        ds["n"] = ds["n"].copy(data=(1.0 + np.exp(-73.0 * (x**2 + y**2 + z**2))))
        ds["phi"] = ds["n"].copy(data=(4.0 * x + 0.0 * y - 3.0 * z))
        ds["V"] = ds["n"].copy()
        ds["V"].data[...] = 5.0

        vr = ds.filament.CoM_velocity("n", "radial", dims=dims)
        vp = ds.filament.CoM_velocity("n", "parallel", dims=dims)
        vb = ds.filament.CoM_velocity("n", "binormal", dims=dims)

        if dims is None:
            npt.assert_allclose(vr, -3.0 * np.ones(5))
            npt.assert_allclose(vp, 5.0 * np.ones(5))
            npt.assert_allclose(vb, -4.0 * np.ones(5))
        elif dims == ["radial", "binormal"]:
            npt.assert_allclose(vr, -3.0 * np.ones([5, 11]))
            npt.assert_allclose(vb, -4.0 * np.ones([5, 11]))
        else:
            raise ValueError(f"Unrecognised value dims={dims}")

    def test_with_origin_at_initial_filament_position_defaults(
        self, tmpdir_factory, bout_xyt_example_files
    ):
        # Create data
        dataset_list = bout_xyt_example_files(
            None, lengths=(2, 10, 10, 10), nxpe=1, nype=1, nt=1, syn_data_type=1
        )

        temp_dir = tmpdir_factory.mktemp("filaments")

        inputpath = create_BOUTinp(temp_dir, realistic_geometry="slab")

        ds = open_stormdataset(dataset_list, inputfilepath=inputpath)

        ds.options["filaments"]["Lx"] = 3.0
        ds.options["filaments"]["Ly"] = 4.0
        ds.options["filaments"]["Lz"] = 5.0
        ds.options.getSection("blob")
        ds.params["rho_s0"] = 2.0

        ds["radial"].data[...] = np.linspace(0.3, 5.7, 10)
        ds["parallel"].data[...] = np.linspace(-3.6, 3.6, 10)
        ds["binormal"].data[...] = np.linspace(0.0, 10.0, 10, endpoint=False)

        ds_out = ds.filament.with_origin_at_initial_filament_position()

        # check original dataset was not modified
        npt.assert_equal(ds["radial"], np.linspace(0.3, 5.7, 10))
        npt.assert_equal(ds["parallel"], np.linspace(-3.6, 3.6, 10))
        npt.assert_equal(ds["binormal"], np.linspace(0.0, 10.0, 10, endpoint=False))

        # Check result
        npt.assert_allclose(ds_out["radial"], np.linspace(-1.2, 4.2, 10), rtol=2.0e-15)
        npt.assert_allclose(
            ds_out["parallel"], np.linspace(-3.6, 3.6, 10), rtol=2.0e-15
        )
        npt.assert_allclose(
            ds_out["binormal"], np.linspace(-5.0, 4.0, 10), rtol=2.0e-15
        )

    def test_with_origin_at_initial_filament_position_offsets(
        self, tmpdir_factory, bout_xyt_example_files
    ):
        # Create data
        dataset_list = bout_xyt_example_files(
            None, lengths=(2, 10, 10, 10), nxpe=1, nype=1, nt=1, syn_data_type=1
        )

        temp_dir = tmpdir_factory.mktemp("filaments")

        inputpath = create_BOUTinp(temp_dir, realistic_geometry="slab")

        ds = open_stormdataset(dataset_list, inputfilepath=inputpath)

        ds.options["filaments"]["Lx"] = 3.0
        ds.options["filaments"]["Ly"] = 4.0
        ds.options["filaments"]["Lz"] = 5.0
        ds.options.getSection("blob")
        ds.options["blob"]["xoffset"] = 0.4
        ds.options["blob"]["yoffset"] = 0.6
        ds.options["blob"]["zoffset"] = 0.7
        ds.params["rho_s0"] = 2.0

        ds["radial"].data[...] = np.linspace(0.3, 5.7, 10)
        ds["parallel"].data[...] = np.linspace(-3.6, 3.6, 10)
        ds["binormal"].data[...] = np.linspace(0.0, 10.0, 10, endpoint=False)

        ds_out = ds.filament.with_origin_at_initial_filament_position()

        # check original dataset was not modified
        npt.assert_equal(ds["radial"], np.linspace(0.3, 5.7, 10))
        npt.assert_equal(ds["parallel"], np.linspace(-3.6, 3.6, 10))
        npt.assert_equal(ds["binormal"], np.linspace(0.0, 10.0, 10, endpoint=False))

        # Check result
        npt.assert_allclose(ds_out["radial"], np.linspace(-2.1, 3.3, 10), rtol=2.0e-15)
        npt.assert_allclose(
            ds_out["parallel"], np.linspace(-4.4, 2.8, 10), rtol=2.0e-15
        )
        npt.assert_allclose(
            ds_out["binormal"], np.linspace(-7.0, 3.0, 10, endpoint=False), rtol=2.0e-15
        )

    def test_time_index_of_max_CoM_velocity(
        self, tmpdir_factory, bout_xyt_example_files
    ):
        # Create data
        dataset_list = bout_xyt_example_files(
            None, lengths=(5, 10, 11, 12), nxpe=1, nype=1, nt=1, syn_data_type=1
        )

        temp_dir = tmpdir_factory.mktemp("filaments")

        inputpath = create_BOUTinp(temp_dir, realistic_geometry="slab")

        ds = open_stormdataset(dataset_list, inputfilepath=inputpath)

        ds["n_eq"] = ds["n"].isel(time=0).copy(data=np.ones([10, 11, 12]))

        n = ds["time"]

        x0 = np.zeros(5)[:, np.newaxis, np.newaxis, np.newaxis]
        y0 = np.zeros(5)[:, np.newaxis, np.newaxis, np.newaxis]
        z0 = np.zeros(5)[:, np.newaxis, np.newaxis, np.newaxis]
        x = np.linspace(0.0, 1.0, 10)[np.newaxis, :, np.newaxis, np.newaxis] - x0
        y = np.linspace(0.0, 1.0, 11)[np.newaxis, np.newaxis, :, np.newaxis] - y0
        z = np.linspace(0.0, 1.0, 12)[np.newaxis, np.newaxis, np.newaxis, :] - z0
        ds["time"].data[:] = np.linspace(0.0, 2.0, 5)
        ds["radial"].data[:] = np.linspace(0.0, 1.0, 10)
        ds["parallel"].data[:] = np.linspace(0.0, 1.0, 11)
        ds["binormal"].data[:] = np.linspace(0.0, 1.0, 12)
        ds["n"] = ds["n"].copy(data=(1.0 + np.exp(-73.0 * (x**2 + y**2 + z**2))))
        ds["phi"] = ds["n"].copy(data=(-4.0 * x + 0.0 * y + 3.0 * z))
        ds["V"] = ds["n"].copy()
        ds["V"].data[...] = 5.0

        # Change phi, V to put maximum in CoM velocity
        x = x[0]
        y = y[0]
        z = z[0]
        ds["phi"][2, ...] = -4.0 * x + 0.0 * y + 4.0 * z
        ds["V"][3, ...] = 6.0
        ds["phi"][4, ...] = -5.0 * x + 0.0 * y + 3.0 * z

        assert ds.filament.time_index_of_max_CoM_velocity("n")["time"].values == 2
        assert (
            ds.filament.time_index_of_max_CoM_velocity("n", direction="parallel")[
                "time"
            ].values
            == 3
        )
        assert (
            ds.filament.time_index_of_max_CoM_velocity("n", direction="binormal")[
                "time"
            ].values
            == 4
        )

    def test_select_at_max_CoM_velocity(self, tmpdir_factory, bout_xyt_example_files):
        # Create data
        dataset_list = bout_xyt_example_files(
            None, lengths=(5, 10, 11, 12), nxpe=1, nype=1, nt=1, syn_data_type=1
        )

        temp_dir = tmpdir_factory.mktemp("filaments")

        inputpath = create_BOUTinp(temp_dir, realistic_geometry="slab")

        ds = open_stormdataset(dataset_list, inputfilepath=inputpath)

        ds["n_eq"] = ds["n"].isel(time=0).copy(data=np.ones([10, 11, 12]))

        n = ds["time"]

        x0 = np.zeros(5)[:, np.newaxis, np.newaxis, np.newaxis]
        y0 = np.zeros(5)[:, np.newaxis, np.newaxis, np.newaxis]
        z0 = np.zeros(5)[:, np.newaxis, np.newaxis, np.newaxis]
        x = np.linspace(0.0, 1.0, 10)[np.newaxis, :, np.newaxis, np.newaxis] - x0
        y = np.linspace(0.0, 1.0, 11)[np.newaxis, np.newaxis, :, np.newaxis] - y0
        z = np.linspace(0.0, 1.0, 12)[np.newaxis, np.newaxis, np.newaxis, :] - z0
        ds["time"].data[:] = np.linspace(0.0, 2.0, 5)
        ds["radial"].data[:] = np.linspace(0.0, 1.0, 10)
        ds["parallel"].data[:] = np.linspace(0.0, 1.0, 11)
        ds["binormal"].data[:] = np.linspace(0.0, 1.0, 12)
        ds["n"] = ds["n"].copy(data=(1.0 + np.exp(-73.0 * (x**2 + y**2 + z**2))))
        ds["phi"] = ds["n"].copy(data=(-4.0 * x + 0.0 * y + 3.0 * z))
        ds["V"] = ds["n"].copy()
        ds["V"].data[...] = 5.0

        # Change phi, V to put maximum in CoM velocity
        x = x[0]
        y = y[0]
        z = z[0]
        ds["phi"][2, ...] = -4.0 * x + 0.0 * y + 4.0 * z
        ds["V"][3, ...] = 6.0
        ds["phi"][4, ...] = -5.0 * x + 0.0 * y + 3.0 * z

        xrt.assert_identical(
            ds.filament.select_at_max_CoM_velocity("n"), ds.isel(time=2)
        )
        xrt.assert_identical(
            ds.filament.select_at_max_CoM_velocity("n", direction="parallel"),
            ds.isel(time=3),
        )
        xrt.assert_identical(
            ds.filament.select_at_max_CoM_velocity("n", direction="binormal"),
            ds.isel(time=4),
        )

    @pytest.mark.parametrize(
        "radial_width", [0.04, pytest.param(None, marks=pytest.mark.long)]
    )
    @pytest.mark.parametrize(
        "parallel_width", [None, pytest.param(0.04, marks=pytest.mark.long)]
    )
    @pytest.mark.parametrize(
        "binormal_width", [0.04, pytest.param(None, marks=pytest.mark.long)]
    )
    def test_interpolate_relative_to_CoM(
        self,
        tmpdir_factory,
        bout_xyt_example_files,
        radial_width,
        parallel_width,
        binormal_width,
    ):
        # Create data
        dataset_list = bout_xyt_example_files(
            None, lengths=(5, 40, 41, 42), nxpe=1, nype=1, nt=1, syn_data_type=1
        )

        temp_dir = tmpdir_factory.mktemp("filaments")

        inputpath = create_BOUTinp(temp_dir, realistic_geometry="slab")

        ds = open_stormdataset(dataset_list, inputfilepath=inputpath)

        ds["n_eq"] = ds["n"].isel(time=0).copy(data=np.ones([40, 41, 42]))

        x0 = np.linspace(0.2, 0.3, 5)[:, np.newaxis, np.newaxis, np.newaxis]
        y0 = np.linspace(0.5, 0.4, 5)[:, np.newaxis, np.newaxis, np.newaxis]
        z0 = np.linspace(0.5, 0.7, 5)[:, np.newaxis, np.newaxis, np.newaxis]
        x = np.linspace(0.0, 1.0, 40)[np.newaxis, :, np.newaxis, np.newaxis] - x0
        y = np.linspace(0.0, 1.0, 41)[np.newaxis, np.newaxis, :, np.newaxis] - y0
        z = np.linspace(0.0, 1.0, 42)[np.newaxis, np.newaxis, np.newaxis, :] - z0
        ds["time"].data[:] = np.linspace(0.0, 2.0, 5)
        ds["radial"].data[:] = np.linspace(0.0, 1.0, 40)
        ds["parallel"].data[:] = np.linspace(0.0, 1.0, 41)
        ds["binormal"].data[:] = np.linspace(0.0, 1.0, 42)
        ds["n"] = ds["n"].copy(data=(1.0 + np.exp(-73.0 * (x**2 + y**2 + z**2))))

        ds_relative = ds.filament.interpolate_relative_to_CoM(
            "n",
            radial_width=radial_width,
            parallel_width=parallel_width,
            binormal_width=binormal_width,
            n_radial=20,
            n_parallel=40,
            n_binormal=10,
        )

        if radial_width is None:
            x_rel = x
        else:
            x_rel = np.linspace(-0.018, 0.018, 20)[
                np.newaxis, :, np.newaxis, np.newaxis
            ]
        if parallel_width is None:
            y_rel = y
        else:
            y_rel = np.linspace(-0.019, 0.019, 40)[
                np.newaxis, np.newaxis, :, np.newaxis
            ]
        if binormal_width is None:
            z_rel = z
        else:
            z_rel = np.linspace(-0.016, 0.016, 10)[
                np.newaxis, np.newaxis, np.newaxis, :
            ]

        test_data = (
            1.0
            + np.exp(-73.0 * (x_rel**2 + y_rel**2 + z_rel**2))
            + 0.0 * ds["time"].values[:, np.newaxis, np.newaxis, np.newaxis]
        )

        npt.assert_allclose(ds_relative["n"], test_data, rtol=0.02)

    @pytest.mark.parametrize("direction", ["radial", "parallel", "binormal"])
    def test_integrate_over_subregion(
        self, tmpdir_factory, bout_xyt_example_files, direction
    ):
        # Create data
        dataset_list = bout_xyt_example_files(
            None, lengths=(2, 10, 12, 14), nxpe=1, nype=1, nt=1, syn_data_type=1
        )

        temp_dir = tmpdir_factory.mktemp("filaments")

        inputpath = create_BOUTinp(temp_dir, realistic_geometry="slab")

        ds = open_stormdataset(dataset_list, inputfilepath=inputpath)

        t = np.arange(2.0)[:, np.newaxis, np.newaxis, np.newaxis]
        x = np.linspace(-4.0 + 5.0 / 10, 6.0 - 5.0 / 10, 10)[
            np.newaxis, :, np.newaxis, np.newaxis
        ]
        y = np.linspace(-3.0 + 5.0 / 12, 7.0 - 5.0 / 12, 12)[
            np.newaxis, np.newaxis, :, np.newaxis
        ]
        z = np.linspace(-2.0 + 5.0 / 14, 8.0 - 5.0 / 14, 14)[
            np.newaxis, np.newaxis, np.newaxis, :
        ]
        ds["radial"].data[:] = x.squeeze()
        ds["parallel"].data[:] = y.squeeze()
        ds["binormal"].data[:] = z.squeeze()
        ds["dx"].data[...] = 1.0
        ds["dy"].data[...] = 10.0 / 12
        ds["dz"].data[...] = 10.0 / 14

        a = 0.0
        b = 0.0
        c = 0.0
        if direction == "radial":
            a = 1.0
        elif direction == "parallel":
            b = 1.0
        elif direction == "binormal":
            c = 1.0
        else:
            assert False

        ds["n"] = ds["n"].copy(
            data=(1.0 + a * (x + 4.0) + b * (y + 3.0) + c * (z + 2.0) + 0.0 * t)
        )

        ds["n_eq"] = ds["n"].isel(time=0).copy(data=np.ones([10, 12, 14]))

        if direction == "radial":
            outside = ds.filament.integrate_over_subregion("n", radial_sel=1.0)
            inside = ds.filament.integrate_over_subregion("n", radial_sel=slice(1.0))
            total = ds.filament.integrate_over_subregion("n", radial_sel=None)
        elif direction == "parallel":
            outside = ds.filament.integrate_over_subregion(
                "n", radial_sel=None, parallel_sel=2.0
            )
            inside = ds.filament.integrate_over_subregion(
                "n", radial_sel=None, parallel_sel=slice(2.0)
            )
            total = ds.filament.integrate_over_subregion(
                "n", radial_sel=None, parallel_sel=None
            )
        elif direction == "binormal":
            outside = ds.filament.integrate_over_subregion(
                "n", radial_sel=None, binormal_sel=3.0
            )
            inside = ds.filament.integrate_over_subregion(
                "n", radial_sel=None, binormal_sel=slice(3.0)
            )
            total = ds.filament.integrate_over_subregion(
                "n", radial_sel=None, binormal_sel=None
            )
        else:
            assert False

        npt.assert_allclose(outside, 7.5 * 500.0 * np.ones(2), rtol=2.0e-15)
        npt.assert_allclose(inside, 2.5 * 500.0 * np.ones(2), rtol=2.0e-15)
        npt.assert_allclose(total, 5.0 * 1000.0 * np.ones(2), rtol=2.0e-15)
