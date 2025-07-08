import pytest

from pathlib import Path
from xbout.tests.test_load import bout_xyt_example_files
import numpy as np
import numpy.testing as npt

from xstorm import open_stormdataset

from .test_stormdataset import create_BOUTinp


class TestMethods:
    @pytest.mark.parametrize(
        "keep_yboundaries", [pytest.param(True, marks=pytest.mark.long), False]
    )
    @pytest.mark.parametrize("myg", [pytest.param(2, marks=pytest.mark.long), 0])
    def test_extrapolate_to_sheath(
        self, tmp_path_factory, bout_xyt_example_files, keep_yboundaries, myg
    ):
        # in principle would be nicer to have several separate tests for different
        # methods, but merging into one here to save time creating the test Dataset

        # Create data
        path = bout_xyt_example_files(
            tmp_path_factory,
            lengths=(2, 3, 5, 2),
            guards={"y": myg},
            nxpe=3,
            nype=6,
            nt=1,
            grid="grid",
            topology="lower-disconnected-double-null",
            write_to_disk=True,
        )

        gridpath = str(Path(path).parent) + "/grid.nc"

        # inputpath = create_BOUTinp(path, realistic_geometry=realistic_geometry)
        inputpath = create_BOUTinp(path)

        # Load it as a stormdataset
        ds = open_stormdataset(
            datapath=path,
            gridfilepath=gridpath,
            inputfilepath=inputpath,
            keep_yboundaries=keep_yboundaries,
        )
        ds["zShift_CELL_YLOW"] = ds["zShift"].copy()
        ds["zShift_CELL_YLOW"].attrs["cell_location"] = "CELL_YLOW"
        ds = ds.set_coords("zShift_CELL_YLOW")

        if keep_yboundaries and myg > 0:
            new_data = np.ones([2, 9, 30 + 4 * myg, 2])

            new_data[:, :, myg, :] = 2.1
            new_data[:, :, myg + 1, :] = 3.01
            new_data[:, :, myg + 2, :] = 4.001

            new_data[:, :, -myg - 3, :] = 2.1
            new_data[:, :, -myg - 2, :] = 3.01
            new_data[:, :, -myg - 1, :] = 4.001
        else:
            new_data = np.ones([2, 9, 30, 2])

            new_data[:, :, 0, :] = 2.1
            new_data[:, :, 1, :] = 3.01
            new_data[:, :, 2, :] = 4.001

            new_data[:, :, -3, :] = 2.1
            new_data[:, :, -2, :] = 3.01
            new_data[:, :, -1, :] = 4.001

        ds["n"] = ds["n"].copy(data=new_data)

        n_extrap = ds["n"].sheath.extrapolate_to_lower()
        lower_bndry_lower_regions = [
            "lower_inner_SOL",
            "lower_inner_intersep",
            "lower_inner_PFR",
        ]
        lower_bndry_upper_regions = [
            "upper_outer_PFR",
            "upper_outer_intersep",
            "upper_outer_SOL",
        ]
        upper_bndry_lower_regions = [
            "lower_outer_PFR",
            "lower_outer_intersep",
            "lower_outer_SOL",
        ]
        upper_bndry_upper_regions = [
            "upper_inner_SOL",
            "upper_inner_intersep",
            "upper_inner_PFR",
        ]
        assert len(n_extrap) == (
            len(lower_bndry_lower_regions) + len(lower_bndry_upper_regions)
        )
        for region in lower_bndry_upper_regions:
            assert np.all(n_extrap[region].values == 1.0)
        for region in lower_bndry_lower_regions:
            # stencil of extrapolation is n_extrap = 3*n[y=0] - 3*n[y=1] + n[y=2]
            boundary_cell = 3.0 * 2.1 - 3.0 * 3.01 + 4.001
            sheath_value = 0.375 * boundary_cell + 0.75 * 2.1 - 0.125 * 3.01
            npt.assert_allclose(n_extrap[region].values, sheath_value, rtol=1.0e-15)

        if keep_yboundaries and myg > 0:
            n_extrap = ds["n"].sheath.extrapolate_to_upper()
            assert len(n_extrap) == (
                len(upper_bndry_lower_regions) + len(upper_bndry_upper_regions)
            )
            for region in upper_bndry_upper_regions:
                assert np.all(n_extrap[region].values == 1.0)
            for region in upper_bndry_lower_regions:
                # stencil of extrapolation is n_extrap = 3*n[y=-1] - 3*n[y=-2] + n[y=-3]
                boundary_cell = 3.0 * 4.001 - 3.0 * 3.01 + 2.1
                sheath_value = 0.375 * boundary_cell + 0.75 * 4.001 - 0.125 * 3.01
                npt.assert_allclose(n_extrap[region].values, sheath_value, rtol=1.0e-15)
        else:
            with pytest.raises(ValueError):
                n_extrap = ds["n"].sheath.extrapolate_to_upper()

        # test extrapolation of staggered field
        if keep_yboundaries and myg > 0:
            stag_data = np.ones([2, 9, 30 + 4 * myg, 2])

            stag_data[:, :, myg, :] = 2.0

            stag_data[:, :, -myg, :] = 3.0
        else:
            stag_data = np.ones([2, 9, 30, 2])

            stag_data[:, :, 0, :] = 2.0

        ds["U"] = ds["n"].copy(data=stag_data)
        ds["U"].attrs["cell_location"] = "CELL_YLOW"

        U_extrap = ds["U"].sheath.extrapolate_to_lower()
        assert len(U_extrap) == (
            len(lower_bndry_lower_regions) + len(lower_bndry_upper_regions)
        )
        for region in lower_bndry_upper_regions:
            assert np.all(U_extrap[region].values == 1.0)
        for region in lower_bndry_lower_regions:
            assert np.all(U_extrap[region].values == 2.0)

        if keep_yboundaries and myg > 0:
            U_extrap = ds["U"].sheath.extrapolate_to_upper()
            assert len(U_extrap) == (
                len(upper_bndry_lower_regions) + len(upper_bndry_upper_regions)
            )
            for region in upper_bndry_upper_regions:
                assert np.all(U_extrap[region].values == 1.0)
            for region in upper_bndry_lower_regions:
                assert np.all(U_extrap[region].values == 3.0)
        else:
            with pytest.raises(ValueError):
                U_extrap = ds["U"].sheath.extrapolate_to_upper()

        # test forced extrapolation of staggered field
        if keep_yboundaries and myg > 0:
            stag_data = np.ones([2, 9, 30 + 4 * myg, 2])

            stag_data[:, :, myg + 1, :] = 2.1
            stag_data[:, :, myg + 2, :] = 3.01
            stag_data[:, :, myg + 3, :] = 4.001

            stag_data[:, :, -myg - 3, :] = 2.1
            stag_data[:, :, -myg - 2, :] = 3.01
            stag_data[:, :, -myg - 1, :] = 4.001
        else:
            stag_data = np.ones([2, 9, 30, 2])

            stag_data[:, :, 1, :] = 2.1
            stag_data[:, :, 2, :] = 3.01
            stag_data[:, :, 3, :] = 4.001

            stag_data[:, :, -3, :] = 2.1
            stag_data[:, :, -2, :] = 3.01
            stag_data[:, :, -1, :] = 4.001

        ds["U"] = ds["n"].copy(data=stag_data)
        ds["U"].attrs["cell_location"] = "CELL_YLOW"

        U_extrap = ds["U"].sheath.extrapolate_to_lower(ylow_force_extrapolate=True)
        assert len(U_extrap) == (
            len(lower_bndry_lower_regions) + len(lower_bndry_upper_regions)
        )
        for region in lower_bndry_upper_regions:
            assert np.all(U_extrap[region].values == 1.0)
        for region in lower_bndry_lower_regions:
            # stencil of extrapolation is U_extrap = 3*U[y=1] - 3*U[y=2] + U[y=3]
            sheath_value = 3.0 * 2.1 - 3.0 * 3.01 + 4.001
            npt.assert_allclose(U_extrap[region].values, sheath_value, rtol=1.0e-15)

        if keep_yboundaries and myg > 0:
            U_extrap = ds["U"].sheath.extrapolate_to_upper(ylow_force_extrapolate=True)
            assert len(U_extrap) == (
                len(upper_bndry_lower_regions) + len(upper_bndry_upper_regions)
            )
            for region in upper_bndry_upper_regions:
                assert np.all(U_extrap[region].values == 1.0)
            for region in upper_bndry_lower_regions:
                # stencil of extrapolation is U_extrap = 3*U[y=-1] - 3*U[y=-2] + U[y=-3]
                sheath_value = 3.0 * 4.001 - 3.0 * 3.01 + 2.1
                npt.assert_allclose(U_extrap[region].values, sheath_value, rtol=1.0e-15)
        else:
            with pytest.raises(ValueError):
                U_extrap = ds["U"].sheath.extrapolate_to_upper(
                    ylow_force_extrapolate=True
                )

        ##########################################
        # test methods that extrapolate logarithms
        ##########################################
        ds["logn"] = ds["n"].copy(data=np.log(ds["n"]))
        ds["logT"] = ds["n"].copy(data=np.log(2.0 * ds["n"]))

        # Note test input file has n_0=1 (in cm^-3) so test dataset has
        # n_0=1e6 (in m^-3).
        # The result from extrapolating logn is multiplied by n_0 because the logarithmic
        # variables are not unnormalised.
        n_sheath = ds.sheath.n_lower()
        assert len(n_sheath) == (
            len(lower_bndry_lower_regions) + len(lower_bndry_upper_regions)
        )
        for region in lower_bndry_upper_regions:
            assert np.all(n_sheath[region].values == 1.0e6)
        for region in lower_bndry_lower_regions:
            # stencil of extrapolation is
            # logn_sheath = 3*logn[y=0] - 3*logn[y=1] + logn[y=2]
            boundary_cell = 2.1**3 / 3.01**3 * 4.001
            sheath_value = 1.0e6 * (
                boundary_cell**0.375 * 2.1**0.75 / 3.01**0.125
            )
            npt.assert_allclose(n_sheath[region].values, sheath_value, rtol=1.0e-15)

        if keep_yboundaries and myg > 0:
            n_sheath = ds.sheath.n_upper()
            assert len(n_sheath) == (
                len(upper_bndry_lower_regions) + len(upper_bndry_upper_regions)
            )
            for region in upper_bndry_upper_regions:
                assert np.all(n_sheath[region].values == 1.0e6)
            for region in upper_bndry_lower_regions:
                # stencil of extrapolation is
                # logn_sheath = 3*logn[y=-1] - 3*logn[y=-2] + logn[y=-3]
                boundary_cell = 4.001**3 / 3.01**3 * 2.1
                sheath_value = 1.0e6 * (
                    boundary_cell**0.375 * 4.001**0.75 / 3.01**0.125
                )
                npt.assert_allclose(n_sheath[region].values, sheath_value, rtol=1.0e-15)
        else:
            with pytest.raises(ValueError):
                n_sheath = ds.sheath.n_upper()

        T_e_sheath = ds.sheath.T_e_lower()
        assert len(T_e_sheath) == (
            len(lower_bndry_lower_regions) + len(lower_bndry_upper_regions)
        )
        for region in lower_bndry_upper_regions:
            npt.assert_allclose(T_e_sheath[region].values, 2.0)
        for region in lower_bndry_lower_regions:
            # stencil of extrapolation is
            # logT_e_sheath = 3*logT_e[y=0] - 3*logT_e[y=1] + logT_e[y=2]
            boundary_cell = 4.2**3 / 6.02**3 * 8.002
            sheath_value = boundary_cell**0.375 * 4.2**0.75 / 6.02**0.125
            npt.assert_allclose(T_e_sheath[region].values, sheath_value, rtol=1.0e-15)

        if keep_yboundaries and myg > 0:
            T_e_sheath = ds.sheath.T_e_upper()
            assert len(T_e_sheath) == (
                len(upper_bndry_lower_regions) + len(upper_bndry_upper_regions)
            )
            for region in upper_bndry_upper_regions:
                npt.assert_allclose(T_e_sheath[region].values, 2.0)
            for region in upper_bndry_lower_regions:
                # stencil of extrapolation is
                # logT_e_sheath = 3*logT_e[y=-1] - 3*logT_e[y=-2] + logT_e[y=-3]
                boundary_cell = 8.002**3 / 6.02**3 * 4.2
                sheath_value = boundary_cell**0.375 * 8.002**0.75 / 6.02**0.125
                npt.assert_allclose(
                    T_e_sheath[region].values, sheath_value, rtol=1.0e-15
                )
        else:
            with pytest.raises(ValueError):
                T_e_sheath = ds.sheath.T_e_upper()

        # first test T_i with no T_i in Dataset - should return 0.0
        T_i_sheath = ds.sheath.T_i_lower()
        assert len(T_i_sheath) == (
            len(lower_bndry_lower_regions) + len(lower_bndry_upper_regions)
        )
        for region in lower_bndry_upper_regions:
            assert T_i_sheath[region] == 0.0
        for region in lower_bndry_lower_regions:
            assert T_i_sheath[region] == 0.0

        T_i_sheath = ds.sheath.T_i_upper()
        assert len(T_i_sheath) == (
            len(upper_bndry_lower_regions) + len(upper_bndry_upper_regions)
        )
        for region in upper_bndry_upper_regions:
            assert T_i_sheath[region] == 0.0
        for region in upper_bndry_lower_regions:
            assert T_i_sheath[region] == 0.0

        # reset logn and logT to simplify tests for calculated sheath quantities
        ds["logn"] = ds["n"].copy(data=1.5 * np.ones(ds["n"].shape))
        ds["logT"] = ds["n"].copy(data=np.ones(ds["n"].shape))

        me_over_mi = 5.485799e-4  # (electron mass)/(1 amu)
        phi_ambi = ds.sheath.ambipolar_potential_lower()
        expected = -np.exp(1.0) * 0.5 * np.log(2.0 * np.pi * me_over_mi)
        for phi in phi_ambi.values():
            npt.assert_allclose(phi, expected)
        if keep_yboundaries and myg > 0:
            phi_ambi = ds.sheath.ambipolar_potential_upper()
            for phi in phi_ambi.values():
                npt.assert_allclose(phi, expected)
        else:
            with pytest.raises(ValueError):
                phi_ambi = ds.sheath.ambipolar_potential_upper()

        me_plus_mi = 1.66145000496e-27  # (electron mass) + (1 amu)
        echarge = 1.60217662e-19
        # Test extrapolated |U| less than sound speed first
        ds["U"] = ds["U"].copy(data=np.zeros(ds["U"].shape))
        U_sheath = ds.sheath.U_lower()
        expected = -np.sqrt(echarge * np.exp(1.0) / me_plus_mi)
        for U in U_sheath.values():
            npt.assert_allclose(U, expected)
        if keep_yboundaries and myg > 0:
            U_sheath = ds.sheath.U_upper()
            expected = np.sqrt(echarge * np.exp(1.0) / me_plus_mi)
            for U in U_sheath.values():
                npt.assert_allclose(U, expected)
        else:
            with pytest.raises(ValueError):
                U_sheath = ds.sheath.U_upper()

        # Test extrapolated |U| greater than sound speed
        ds["U"] = ds["U"].copy(data=-1.0e9 * np.ones(ds["U"].shape))
        U_sheath = ds.sheath.U_lower()
        for U in U_sheath.values():
            npt.assert_allclose(U, -1.0e9)
        ds["U"] = ds["U"].copy(data=1.0e9 * np.ones(ds["U"].shape))
        if keep_yboundaries and myg > 0:
            U_sheath = ds.sheath.U_upper()
            for U in U_sheath.values():
                npt.assert_allclose(U, 1.0e9)
        else:
            with pytest.raises(ValueError):
                U_sheath = ds.sheath.U_upper()

        e_over_mi = 96485332.9
        ds["phi"] = ds["n"].copy(data=3.0 * np.ones(ds["n"].shape))
        ds["V"] = ds["U"].copy(data=np.zeros(ds["U"].shape))
        V_sheath = ds.sheath.V_lower()
        V_prefactor = np.sqrt(
            1.0
            / 2.0
            / np.pi
            / me_over_mi
            / (1.0 + me_over_mi)
            * e_over_mi
            * np.exp(1.0)
        )
        exponential = np.exp(-3.0 / np.exp(1.0))
        for V in V_sheath.values():
            npt.assert_allclose(V, -V_prefactor * exponential)
        if keep_yboundaries and myg > 0:
            V_sheath = ds.sheath.V_upper()
            for V in V_sheath.values():
                npt.assert_allclose(V, V_prefactor * exponential)
        else:
            with pytest.raises(ValueError):
                V_sheath = ds.sheath.V_upper()

        # if phi is negative, exponential term is limited to 1.0
        ds["phi"] = ds["n"].copy(data=-0.1 * np.ones(ds["n"].shape))
        V_sheath = ds.sheath.V_lower()
        V_prefactor = np.sqrt(
            1.0
            / 2.0
            / np.pi
            / me_over_mi
            / (1.0 + me_over_mi)
            * e_over_mi
            * np.exp(1.0)
        )
        exponential = 1.0
        for V in V_sheath.values():
            npt.assert_allclose(V, -V_prefactor * exponential)
        if keep_yboundaries and myg > 0:
            V_sheath = ds.sheath.V_upper()
            for V in V_sheath.values():
                npt.assert_allclose(V, V_prefactor * exponential)
        else:
            with pytest.raises(ValueError):
                V_sheath = ds.sheath.V_upper()

        ds["phi"] = ds["n"].copy(data=np.ones(ds["n"].shape))
        ds["U"] = ds["U"].copy(data=np.zeros(ds["U"].shape))
        J_sheath = ds.sheath.current_lower()
        exponential = np.exp(-1.0 / np.exp(1.0))
        for J in J_sheath.values():
            npt.assert_allclose(
                J,
                -echarge
                * np.exp(1.5)
                * 1.0e6
                * (
                    np.sqrt(echarge * np.exp(1.0) / me_plus_mi)
                    - V_prefactor * exponential
                ),
            )
        if keep_yboundaries and myg > 0:
            J_sheath = ds.sheath.current_upper()
            for J in J_sheath.values():
                npt.assert_allclose(
                    J,
                    echarge
                    * np.exp(1.5)
                    * 1.0e6
                    * (
                        np.sqrt(echarge * np.exp(1.0) / me_plus_mi)
                        - V_prefactor * exponential
                    ),
                )
        else:
            with pytest.raises(ValueError):
                J_sheath = ds.sheath.current_upper()

        qpar_sheath = ds.sheath.qpar_lower()
        for qpar in qpar_sheath.values():
            npt.assert_allclose(
                qpar,
                -(
                    (-0.5 * np.log(2.0 * np.pi * me_over_mi) - 0.5)
                    * echarge
                    * np.exp(1.0)
                    - 0.5 * 9.10938356e-31 * (V_prefactor * exponential) ** 2
                )
                * 1.0e6
                * np.exp(1.5)
                * V_prefactor
                * exponential,
            )
        if keep_yboundaries and myg > 0:
            qpar_sheath = ds.sheath.qpar_upper()
            for qpar in qpar_sheath.values():
                npt.assert_allclose(
                    qpar,
                    (
                        (-0.5 * np.log(2.0 * np.pi * me_over_mi) - 0.5)
                        * echarge
                        * np.exp(1.0)
                        - 0.5 * 9.10938356e-31 * (V_prefactor * exponential) ** 2
                    )
                    * 1.0e6
                    * np.exp(1.5)
                    * V_prefactor
                    * exponential,
                )
        else:
            with pytest.raises(ValueError):
                qpar_sheath = ds.sheath.qpar_upper()

        Qpar_sheath = ds.sheath.electron_total_energy_flux_lower()
        for Qpar in Qpar_sheath.values():
            npt.assert_allclose(
                Qpar,
                -(
                    (-0.5 * np.log(2.0 * np.pi * me_over_mi) - 0.5)
                    * echarge
                    * np.exp(1.0)
                    - 0.5 * 9.10938356e-31 * (V_prefactor * exponential) ** 2
                )
                * 1.0e6
                * np.exp(1.5)
                * V_prefactor
                * exponential
                - 5.0
                / 2.0
                * 1.0e6
                * np.exp(1.5)
                * echarge
                * np.exp(1.0)
                * V_prefactor
                * exponential
                - 1.0
                / 2.0
                * 9.10938356e-31
                * 1.0e6
                * np.exp(1.5)
                * (V_prefactor * exponential) ** 3,
            )
        if keep_yboundaries and myg > 0:
            Qpar_sheath = ds.sheath.electron_total_energy_flux_upper()
            for Qpar in Qpar_sheath.values():
                npt.assert_allclose(
                    Qpar,
                    (
                        (-0.5 * np.log(2.0 * np.pi * me_over_mi) - 0.5)
                        * echarge
                        * np.exp(1.0)
                        - 0.5 * 9.10938356e-31 * (V_prefactor * exponential) ** 2
                    )
                    * 1.0e6
                    * np.exp(1.5)
                    * V_prefactor
                    * exponential
                    + 5.0
                    / 2.0
                    * 1.0e6
                    * np.exp(1.5)
                    * echarge
                    * np.exp(1.0)
                    * V_prefactor
                    * exponential
                    + 1.0
                    / 2.0
                    * 9.10938356e-31
                    * 1.0e6
                    * np.exp(1.5)
                    * (V_prefactor * exponential) ** 3,
                )
        else:
            with pytest.raises(ValueError):
                Qpar_sheath = ds.sheath.electron_total_energy_flux_upper()

        ds["T_i"] = ds["T"]
        ds["logTi"] = ds["n"].copy(data=np.log(3.0 * ds["n"]))
        T_i_sheath = ds.sheath.T_i_lower()
        assert len(T_i_sheath) == (
            len(lower_bndry_lower_regions) + len(lower_bndry_upper_regions)
        )
        for region in lower_bndry_upper_regions:
            npt.assert_allclose(T_i_sheath[region].values, 3.0)
        for region in lower_bndry_lower_regions:
            # stencil of extrapolation is
            # logT_i_sheath = 3*logT_i[y=0] - 3*logT_i[y=1] + logT_i[y=2]
            boundary_cell = 6.3**3 / 9.03**3 * 12.003
            sheath_value = boundary_cell**0.375 * 6.3**0.75 / 9.03**0.125
            npt.assert_allclose(T_i_sheath[region].values, sheath_value, rtol=2.0e-15)

        if keep_yboundaries and myg > 0:
            T_i_sheath = ds.sheath.T_i_upper()
            assert len(T_i_sheath) == (
                len(upper_bndry_lower_regions) + len(upper_bndry_upper_regions)
            )
            for region in upper_bndry_upper_regions:
                npt.assert_allclose(T_i_sheath[region].values, 3.0)
            for region in upper_bndry_lower_regions:
                # stencil of extrapolation is
                # logT_i_sheath = 3*logT_i[y=-1] - 3*logT_i[y=-2] + logT_i[y=-3]
                boundary_cell = 12.003**3 / 9.03**3 * 6.3
                sheath_value = boundary_cell**0.375 * 12.003**0.75 / 9.03**0.125
                npt.assert_allclose(
                    T_i_sheath[region].values, sheath_value, rtol=1.0e-15
                )
        else:
            with pytest.raises(ValueError):
                T_i_sheath = ds.sheath.T_i_upper()

        # reset logTi to simplify test for ambipolar_phi_sheath
        ds["logTi"] = ds["n"].copy(data=2.0 * np.ones(ds["n"].shape))
        phi_ambi = ds.sheath.ambipolar_potential_lower()
        expected = (
            -np.exp(1.0) * 0.5 * np.log(2.0 * np.pi * me_over_mi * (1.0 + np.exp(1.0)))
        )
        for phi in phi_ambi.values():
            npt.assert_allclose(phi, expected)
        if keep_yboundaries and myg > 0:
            phi_ambi = ds.sheath.ambipolar_potential_upper()
            for phi in phi_ambi.values():
                npt.assert_allclose(phi, expected)
        else:
            with pytest.raises(ValueError):
                phi_ambi = ds.sheath.ambipolar_potential_upper()

        # Test extrapolated |U| less than sound speed first.
        ds["U"] = ds["U"].copy(data=np.zeros(ds["U"].shape))
        U_sheath = ds.sheath.U_lower()
        expected = -np.sqrt(echarge * (np.exp(1.0) + np.exp(2.0)) / me_plus_mi)
        for U in U_sheath.values():
            npt.assert_allclose(U, expected)
        if keep_yboundaries and myg > 0:
            U_sheath = ds.sheath.U_upper()
            expected = np.sqrt(echarge * (np.exp(1.0) + np.exp(2.0)) / me_plus_mi)
            for U in U_sheath.values():
                npt.assert_allclose(U, expected)
        else:
            with pytest.raises(ValueError):
                U_sheath = ds.sheath.U_upper()

        # |U| greater than sound speed shouldn't really be affected, so don't test that
        # again
