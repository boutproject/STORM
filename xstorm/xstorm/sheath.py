from xarray import register_dataset_accessor, register_dataarray_accessor

import numpy as np

from xstorm import StormDatasetAccessor, StormDataArrayAccessor

from .plotfuncs import plot_sheath_currents
from .sheathutils import (
    apply_to_all_regions_lower,
    apply_to_all_regions_upper,
    ambipolar_potential_lowerupper,
    lower_from_field_aligned,
    upper_from_field_aligned,
)


@register_dataset_accessor("sheath")
class SheathDatasetAccessor(StormDatasetAccessor):
    """
    Calculate values and fluxes at sheath boundaries, replicating algorithms used in
    STORM.

    These methods are on a Dataset because they implement analysis which require multiple
    variables.
    """

    def __init__(self, ds):
        super().__init__(ds)

    def ambipolar_potential_lower(self, region=None, return_aligned=False):
        """
        Calculate ambipolar sheath potential at lower sheaths

        Parameters
        ----------
        region : str, optional
            Get the lower sheath value for a specific region. By default looks for all
            regions.
        return_aligned : bool, default False
            Return result on field-aligned grid (note, if input is field-aligned already
            then the result is always on the field-aligned grid)
        """
        return ambipolar_potential_lowerupper(
            self, region, "lower", return_aligned=return_aligned
        )

    def ambipolar_potential_upper(self, region=None, return_aligned=False):
        """
        Calculate ambipolar sheath potential at upper sheaths

        Parameters
        ----------
        region : str, optional
            Get the upper sheath value for a specific region. By default looks for all
            regions.
        return_aligned : bool, default False
            Return result on field-aligned grid (note, if input is field-aligned already
            then the result is always on the field-aligned grid)
        """
        return ambipolar_potential_lowerupper(
            self, region, "upper", return_aligned=return_aligned
        )

    def n_lower(self, region=None, *, return_aligned=False):
        """
        Use extrapolation replicated from STORM to calculate the density at lower (y
        increasing away from the boundary) sheath or sheaths.

        Parameters
        ----------
        region : str, optional
            Get the lower sheath value for a specific region. By default looks for all
            regions.
        return_aligned : bool, default False
            Return result on field-aligned grid (note, if input is field-aligned already
            then the result is always on the field-aligned grid)

        Returns
        -------
        DataArray or dict of DataArray
            If called with a 'region' argument or if only one sheath exists, then
            returns a DataArray with the result for that sheath. Otherwise returns a
            dict whose keys are region names and values are result for that region
        """
        if self.data.metadata["storm_normalised"]:
            raise ValueError(
                "n_lower() is only implemented for unnormalised Datasets so far"
            )

        if region is None:
            return apply_to_all_regions_lower(
                self.data, self.n_lower, return_aligned=return_aligned
            )

        ds = self.data

        return (
            np.exp(
                ds["logn"].sheath.extrapolate_to_lower(
                    region, return_aligned=return_aligned
                )
            )
            * ds.params["n_0"]
        )

    def n_upper(self, region=None, *, return_aligned=False):
        """
        Use extrapolation replicated from STORM to calculate the density at upper (y
        increasing towards the boundary) sheath or sheaths.

        Parameters
        ----------
        region : str, optional
            Get the upper sheath value for a specific region. By default looks for all
            regions.

        Returns
        -------
        DataArray or dict of DataArray
            If called with a 'region' argument or if only one sheath exists, then
            returns a DataArray with the result for that sheath. Otherwise returns a
            dict whose keys are region names and values are result for that region
        """
        if self.data.metadata["storm_normalised"]:
            raise ValueError(
                "n_upper() is only implemented for unnormalised Datasets so far"
            )

        if region is None:
            return apply_to_all_regions_upper(
                self.data, self.n_upper, return_aligned=return_aligned
            )

        ds = self.data

        return (
            np.exp(
                ds["logn"].sheath.extrapolate_to_upper(
                    region, return_aligned=return_aligned
                )
            )
            * ds.params["n_0"]
        )

    def T_e_lower(self, region=None, *, return_aligned=False):
        """
        Use extrapolation replicated from STORM to calculate the electron temperature at
        lower (y increasing away from the boundary) sheath or sheaths.

        Parameters
        ----------
        region : str, optional
            Get the lower sheath value for a specific region. By default looks for all
            regions.
        return_aligned : bool, default False
            Return result on field-aligned grid (note, if input is field-aligned already
            then the result is always on the field-aligned grid)

        Returns
        -------
        DataArray or dict of DataArray
            If called with a 'region' argument or if only one sheath exists, then
            returns a DataArray with the result for that sheath. Otherwise returns a
            dict whose keys are region names and values are result for that region
        """
        if self.data.metadata["storm_normalised"]:
            raise ValueError(
                "T_e_lower() is only implemented for unnormalised Datasets so far"
            )

        if region is None:
            return apply_to_all_regions_lower(
                self.data, self.T_e_lower, return_aligned=return_aligned
            )

        ds = self.data

        return (
            np.exp(
                ds["logT"].sheath.extrapolate_to_lower(
                    region, return_aligned=return_aligned
                )
            )
            * ds.params["T_e0"]
        )

    def T_e_upper(self, region=None, *, return_aligned=False):
        """
        Use extrapolation replicated from STORM to calculate the electron temperature at
        upper (y increasing towards the boundary) sheath or sheaths.

        Parameters
        ----------
        region : str, optional
            Get the upper sheath value for a specific region. By default looks for all
            regions.
        return_aligned : bool, default False
            Return result on field-aligned grid (note, if input is field-aligned already
            then the result is always on the field-aligned grid)

        Returns
        -------
        DataArray or dict of DataArray
            If called with a 'region' argument or if only one sheath exists, then
            returns a DataArray with the result for that sheath. Otherwise returns a
            dict whose keys are region names and values are result for that region
        """
        if self.data.metadata["storm_normalised"]:
            raise ValueError(
                "T_e_upper() is only implemented for unnormalised Datasets so far"
            )

        if region is None:
            return apply_to_all_regions_upper(
                self.data, self.T_e_upper, return_aligned=return_aligned
            )

        ds = self.data

        return (
            np.exp(
                ds["logT"].sheath.extrapolate_to_upper(
                    region, return_aligned=return_aligned
                )
            )
            * ds.params["T_e0"]
        )

    def T_i_lower(self, region=None, *, return_aligned=False):
        """
        Use extrapolation replicated from STORM to calculate the ion temperature at lower
        (y increasing away from the boundary) sheath or sheaths.

        Parameters
        ----------
        region : str, optional
            Get the lower sheath value for a specific region. By default looks for all
            regions.
        return_aligned : bool, default False
            Return result on field-aligned grid (note, if input is field-aligned already
            then the result is always on the field-aligned grid)

        Returns
        -------
        DataArray or dict of DataArray
            If called with a 'region' argument or if only one sheath exists, then
            returns a DataArray with the result for that sheath. Otherwise returns a
            dict whose keys are region names and values are result for that region
        """
        if self.data.metadata["storm_normalised"]:
            raise ValueError(
                "T_i_lower() is only implemented for unnormalised Datasets so far"
            )

        if region is None:
            return apply_to_all_regions_lower(
                self.data, self.T_i_lower, return_aligned=return_aligned
            )

        ds = self.data

        if "T_i" in ds:
            return (
                np.exp(
                    ds["logTi"].sheath.extrapolate_to_lower(
                        region, return_aligned=return_aligned
                    )
                )
                * ds.params["T_e0"]
            )
        elif (
            "finite_Ti" not in ds.options["filaments"]
            or ds.options["filaments"]["finite_Ti"] == "none"
        ):
            return 0.0
        else:
            # Isothermal, get constant value from options
            if "T_i" in ds.options and "function" in ds.options["T_i"]:
                return ds.options["T_i"].evaluate_scalar("function") * ds.params["T_e0"]
            else:
                return ds.params["T_e0"]

    def T_i_upper(self, region=None, *, return_aligned=False):
        """
        Use extrapolation replicated from STORM to calculate the ion temperature at upper
        (y increasing towards the boundary) sheath or sheaths.

        Parameters
        ----------
        region : str, optional
            Get the upper sheath value for a specific region. By default looks for all
            regions.
        return_aligned : bool, default False
            Return result on field-aligned grid (note, if input is field-aligned already
            then the result is always on the field-aligned grid)

        Returns
        -------
        DataArray or dict of DataArray
            If called with a 'region' argument or if only one sheath exists, then
            returns a DataArray with the result for that sheath. Otherwise returns a
            dict whose keys are region names and values are result for that region
        """
        if self.data.metadata["storm_normalised"]:
            raise ValueError(
                "T_i_upper() is only implemented for unnormalised Datasets so far"
            )

        if region is None:
            return apply_to_all_regions_upper(
                self.data, self.T_i_upper, return_aligned=return_aligned
            )

        ds = self.data

        if "T_i" in ds:
            return (
                np.exp(
                    ds["logTi"].sheath.extrapolate_to_upper(
                        region, return_aligned=return_aligned
                    )
                )
                * ds.params["T_e0"]
            )
        elif (
            "finite_Ti" not in ds.options["filaments"]
            or ds.options["filaments"]["finite_Ti"] == "none"
        ):
            return 0.0
        else:
            # Isothermal, get constant value from options
            if "T_i" in ds.options and "function" in ds.options["T_i"]:
                return ds.options["T_i"].evaluate_scalar("function") * ds.params["T_e0"]
            else:
                return ds.params["T_e0"]

    def U_lower(self, region=None, *, return_aligned=False):
        """
        Use boundary condition replicated from STORM to calculate the ion velocity at
        lower (y increasing away from the boundary) sheath or sheaths. Needed because
        parallel boundary condition is only applied to U_aligned, which is not always
        saved.

        Parameters
        ----------
        region : str, optional
            Get the lower sheath value for a specific region. By default looks for all
            regions.
        return_aligned : bool, default False
            Return result on field-aligned grid (note, if input is field-aligned already
            then the result is always on the field-aligned grid)

        Returns
        -------
        DataArray or dict of DataArray
            If called with a 'region' argument or if only one sheath exists, then
            returns a DataArray with the result for that sheath. Otherwise returns a
            dict whose keys are region names and values are result for that region
        """
        if self.data.metadata["storm_normalised"]:
            raise ValueError(
                "U_lower() is only implemented for unnormalised Datasets so far"
            )

        if region is None:
            return apply_to_all_regions_lower(
                self.data, self.U_lower, return_aligned=return_aligned
            )

        ds = self.data
        ycoord = ds.metadata["bout_ydim"]

        T_e = self.T_e_lower(region, return_aligned=True)
        T_i = self.T_i_lower(region, return_aligned=True)

        U_extrap = ds["U"].sheath.extrapolate_to_lower(
            region, ylow_force_extrapolate=True, return_aligned=True
        )

        e = ds.params["e"]
        m_i = ds.params["m_i"]
        m_e = ds.params["m_e"]

        sound_speed = -np.sqrt(e * (T_e + T_i) / (m_i + m_e))

        # sheath value is sound_speed when magnitude of sound_speed is greater than
        # magnitude is greater than magnitude of U_extrap, otherwise is U_extrap
        U_sheath = sound_speed.where(sound_speed < U_extrap, U_extrap)

        U_sheath.attrs["standard_name"] = "U_sheath lower"
        U_sheath.attrs["long_name"] = "Ion parallel velocity at lower sheath entrance"
        U_sheath.attrs["units"] = "ms^-1"
        U_sheath.attrs["direction_y"] = "Aligned"
        U_sheath.attrs["regions"] = {region: ds.regions[region]}
        U_sheath.attrs["metadata"] = ds.metadata

        if not return_aligned and "zShift" in self.data.coords:
            U_sheath = lower_from_field_aligned(ds["U"].sheath, U_sheath)

        return U_sheath

    def U_upper(self, region=None, *, return_aligned=False):
        """
        Use boundary condition replicated from STORM to calculate the ion velocity at
        upper (y increasing towards the boundary) sheath or sheaths. Needed because
        parallel boundary condition is only applied to U_aligned, which is not always
        saved.

        Parameters
        ----------
        region : str, optional
            Get the upper sheath value for a specific region. By default looks for all
            regions.
        return_aligned : bool, default False
            Return result on field-aligned grid (note, if input is field-aligned already
            then the result is always on the field-aligned grid)

        Returns
        -------
        DataArray or dict of DataArray
            If called with a 'region' argument or if only one sheath exists, then
            returns a DataArray with the result for that sheath. Otherwise returns a
            dict whose keys are region names and values are result for that region
        """
        if self.data.metadata["storm_normalised"]:
            raise ValueError(
                "U_upper() is only implemented for unnormalised Datasets so far"
            )

        if region is None:
            return apply_to_all_regions_upper(
                self.data, self.U_upper, return_aligned=return_aligned
            )

        ds = self.data
        ycoord = ds.metadata["bout_ydim"]

        T_e = self.T_e_upper(region, return_aligned=True)
        T_i = self.T_i_upper(region, return_aligned=True)

        U_extrap = ds["U"].sheath.extrapolate_to_upper(
            region, ylow_force_extrapolate=True, return_aligned=True
        )

        e = ds.params["e"]
        m_i = ds.params["m_i"]
        m_e = ds.params["m_e"]

        sound_speed = np.sqrt(e * (T_e + T_i) / (m_i + m_e))

        # sheath value is sound_speed when magnitude of sound_speed is greater than
        # magnitude is greater than magnitude of U_extrap, otherwise is U_extrap
        U_sheath = sound_speed.where(sound_speed > U_extrap, U_extrap)

        U_sheath.attrs["standard_name"] = "U_sheath upper"
        U_sheath.attrs["long_name"] = "Ion parallel velocity at upper sheath entrance"
        U_sheath.attrs["units"] = "ms^-1"
        U_sheath.attrs["direction_y"] = "Aligned"
        U_sheath.attrs["regions"] = {region: ds.regions[region]}
        U_sheath.attrs["metadata"] = ds.metadata

        if not return_aligned and "zShift" in self.data.coords:
            U_sheath = upper_from_field_aligned(ds["U"].sheath, U_sheath)

        return U_sheath

    def V_lower(self, region=None, *, return_aligned=False):
        """
        Use boundary condition replicated from STORM to calculate the electron velocity
        at lower (y increasing away from the boundary) sheath or sheaths. Needed because
        parallel boundary condition is only applied to V_aligned, which is not always
        saved.

        Assumes "phi_wall" is zero.

        Parameters
        ----------
        region : str, optional
            Get the lower sheath value for a specific region. By default looks for all
            regions.
        return_aligned : bool, default False
            Return result on field-aligned grid (note, if input is field-aligned already
            then the result is always on the field-aligned grid)

        Returns
        -------
        DataArray or dict of DataArray
            If called with a 'region' argument or if only one sheath exists, then
            returns a DataArray with the result for that sheath. Otherwise returns a
            dict whose keys are region names and values are result for that region
        """
        if self.data.metadata["storm_normalised"]:
            raise ValueError(
                "V_lower() is only implemented for unnormalised Datasets so far"
            )

        if "phi_wall" in self.data.settings:
            if self.data.settings["phi_wall"] != 0.0:
                raise ValueError(
                    f"Non-zero phi_wall is not supported. Got "
                    f"phi_wall={self.data.settings['phi_wall']}."
                )

        if region is None:
            return apply_to_all_regions_lower(
                self.data, self.V_lower, return_aligned=return_aligned
            )

        ds = self.data
        ycoord = ds.metadata["bout_ydim"]

        phi = ds["phi"].sheath.extrapolate_to_lower(region, return_aligned=True)
        T_e = self.T_e_lower(region, return_aligned=True)

        e = ds.params["e"]
        m_i = ds.params["m_i"]
        m_e = ds.params["m_e"]

        V_BC_prefactor = np.sqrt(m_i / (2.0 * np.pi * m_e) / (1.0 + m_e / m_i))
        exponential_factor = np.exp(-phi / T_e)
        # Limit the maximum electron flux into the sheath.
        # When phisheath < phi_wall the sheath attracts electrons, and the expression we
        # are using, derived for a strongly electron-repelling sheath, is not valid.
        exponential_factor = exponential_factor.where(exponential_factor < 1.0, 1.0)

        V_sheath = -V_BC_prefactor * np.sqrt(e * T_e / m_i) * exponential_factor

        V_sheath.attrs["standard_name"] = "V_sheath lower"
        V_sheath.attrs[
            "long_name"
        ] = "Electron parallel velocity at lower sheath entrance"
        V_sheath.attrs["units"] = "ms^-1"
        V_sheath.attrs["direction_y"] = "Aligned"
        V_sheath.attrs["regions"] = {region: ds.regions[region]}
        V_sheath.attrs["metadata"] = ds.metadata

        if not return_aligned and "zShift" in self.data.coords:
            V_sheath = lower_from_field_aligned(ds["V"].sheath, V_sheath)

        return V_sheath

    def V_upper(self, region=None, *, return_aligned=False):
        """
        Use boundary condition replicated from STORM to calculate the electron velocity
        at upper (y increasing towards the boundary) sheath or sheaths. Needed because
        parallel boundary condition is only applied to V_aligned, which is not always
        saved.

        Assumes "phi_wall" is zero.

        Parameters
        ----------
        region : str, optional
            Get the upper sheath value for a specific region. By default looks for all
            regions.
        return_aligned : bool, default False
            Return result on field-aligned grid (note, if input is field-aligned already
            then the result is always on the field-aligned grid)

        Returns
        -------
        DataArray or dict of DataArray
            If called with a 'region' argument or if only one sheath exists, then
            returns a DataArray with the result for that sheath. Otherwise returns a
            dict whose keys are region names and values are result for that region
        """
        if self.data.metadata["storm_normalised"]:
            raise ValueError(
                "V_upper() is only implemented for unnormalised Datasets so far"
            )

        if "phi_wall" in self.data.settings:
            if self.data.settings["phi_wall"] != 0.0:
                raise ValueError(
                    f"Non-zero phi_wall is not supported. Got "
                    f"phi_wall={self.data.settings['phi_wall']}"
                )

        if region is None:
            return apply_to_all_regions_upper(
                self.data, self.V_upper, return_aligned=return_aligned
            )

        ds = self.data
        ycoord = ds.metadata["bout_ydim"]

        phi = ds["phi"].sheath.extrapolate_to_upper(region, return_aligned=True)
        T_e = self.T_e_upper(region, return_aligned=True)

        e = ds.params["e"]
        m_i = ds.params["m_i"]
        m_e = ds.params["m_e"]

        V_BC_prefactor = np.sqrt(m_i / (2.0 * np.pi * m_e) / (1.0 + m_e / m_i))
        exponential_factor = np.exp(-phi / T_e)
        # Limit the maximum electron flux into the sheath.
        # When phisheath < phi_wall the sheath attracts electrons, and the expression we
        # are using, derived for a strongly electron-repelling sheath, is not valid.
        exponential_factor = exponential_factor.where(exponential_factor < 1.0, 1.0)

        V_sheath = V_BC_prefactor * np.sqrt(e * T_e / m_i) * exponential_factor

        V_sheath.attrs["standard_name"] = "V_sheath upper"
        V_sheath.attrs[
            "long_name"
        ] = "Electron parallel velocity at upper sheath entrance"
        V_sheath.attrs["units"] = "ms^-1"
        V_sheath.attrs["direction_y"] = "Aligned"
        V_sheath.attrs["regions"] = {region: ds.regions[region]}
        V_sheath.attrs["metadata"] = ds.metadata

        if not return_aligned and "zShift" in self.data.coords:
            V_sheath = upper_from_field_aligned(ds["V"].sheath, V_sheath)

        return V_sheath

    def current_lower(self, region=None, *, return_aligned=False):
        """
        Parallel current at the lower sheaths.

        Parameters
        ----------
        region : str, optional
            Get the lower sheath value for a specific region. By default looks for all
            regions.
        return_aligned : bool, default False
            Return result on field-aligned grid (note, if input is field-aligned already
            then the result is always on the field-aligned grid)

        Returns
        -------
        DataArray or dict of DataArray
            If called with a 'region' argument or if only one sheath exists, then
            returns a DataArray with the result for that sheath. Otherwise returns a
            dict whose keys are region names and values are result for that region
        """
        if self.data.metadata["storm_normalised"]:
            raise ValueError(
                "current_lower() is only implemented for unnormalised Datasets so far"
            )

        if region is None:
            return apply_to_all_regions_lower(
                self.data, self.current_lower, return_aligned=return_aligned
            )

        ds = self.data

        e = ds.params["e"]
        n = self.n_lower(region, return_aligned=return_aligned)
        U = self.U_lower(region, return_aligned=return_aligned)
        V = self.V_lower(region, return_aligned=return_aligned)

        J = e * n * (U - V)

        J.attrs["standard_name"] = "sheath current lower"
        J.attrs["long_name"] = "Parallel current at lower sheath"
        J.attrs["units"] = "Am^-2"

        return J

    def current_upper(self, region=None, *, return_aligned=False):
        """
        Parallel current at the upper sheaths.

        Parameters
        ----------
        region : str, optional
            Get the upper sheath value for a specific region. By default looks for all
            regions.
        return_aligned : bool, default False
            Return result on field-aligned grid (note, if input is field-aligned already
            then the result is always on the field-aligned grid)

        Returns
        -------
        DataArray or dict of DataArray
            If called with a 'region' argument or if only one sheath exists, then
            returns a DataArray with the result for that sheath. Otherwise returns a
            dict whose keys are region names and values are result for that region
        """
        if self.data.metadata["storm_normalised"]:
            raise ValueError(
                "current_upper() is only implemented for unnormalised Datasets so far"
            )

        if region is None:
            return apply_to_all_regions_upper(
                self.data, self.current_upper, return_aligned=return_aligned
            )

        ds = self.data

        e = ds.params["e"]
        n = self.n_upper(region, return_aligned=return_aligned)
        U = self.U_upper(region, return_aligned=return_aligned)
        V = self.V_upper(region, return_aligned=return_aligned)

        J = e * n * (U - V)

        J.attrs["standard_name"] = "sheath current upper"
        J.attrs["long_name"] = "Parallel current at upper sheath"
        J.attrs["units"] = "Am^-2"

        return J

    def plot_currents(self, *, output="total", lowerupper="both", ax=None, **kwargs):
        """
        Calculate parallel currents into the sheaths.

        Positive currents in the output are out of the domain, i.e. currents at lower
        sheaths are multiplied by -1 for output.

        Parameters
        ----------
        output : {["total"], "positive_negative", "animate2d"}
            - "total" : plot total current through each sheath (and sum over all
              sheaths) vs. time
            - "positive_negative": integrate the positive and negative parts of the
              current separately, and plot both
            - "animate2D": make a 2D animation of the current through each sheath
        lowerupper : str {["both"], "lower", "upper"}
            Include both sheaths in the analysis, or only lower (y increasing away from
            the boundary) or upper (y increasing towards the boundary) sheaths
        ax : Matplotlib.axes.Axes
            Axes object to plot on. If not passed in, new figure is created. Ignored if
            output=="animate2D"
        **kwargs
            Passed to animate_list if output=="animate2D"
        """
        return plot_sheath_currents(
            self, output=output, lowerupper=lowerupper, ax=ax, **kwargs
        )

    def qpar_lower(self, region=None, *, return_aligned=False):
        """
        Conductive electron heat flux at the lower sheaths.

        Parameters
        ----------
        region : str, optional
            Get the lower sheath value for a specific region. By default looks for all
            regions.
        return_aligned : bool, default False
            Return result on field-aligned grid (note, if input is field-aligned already
            then the result is always on the field-aligned grid)

        Returns
        -------
        DataArray or dict of DataArray
            If called with a 'region' argument or if only one sheath exists, then
            returns a DataArray with the result for that sheath. Otherwise returns a
            dict whose keys are region names and values are result for that region
        """
        if self.data.metadata["storm_normalised"]:
            raise ValueError(
                "qpar_lower() is only implemented for unnormalised Datasets so far"
            )

        if region is None:
            return apply_to_all_regions_lower(
                self.data, self.qpar_lower, return_aligned=return_aligned
            )

        ds = self.data

        echarge = ds.params["e"]
        me = ds.params["m_e"]
        eV_float_over_T = ds.params["dimless"]["eV_float_over_T"]

        T_e = self.T_e_lower(region, return_aligned=True)
        V = self.V_lower(region, return_aligned=True)
        n = self.n_lower(region, return_aligned=True)

        qpar_sheath = (
            ((np.abs(eV_float_over_T) - 0.5) * echarge * T_e - 0.5 * me * V**2)
            * n
            * V
        )

        qpar_sheath.attrs["standard_name"] = "qpar_sheath lower"
        qpar_sheath.attrs[
            "long_name"
        ] = "Conductive electron parallel heat flux at lower sheath entrance"
        qpar_sheath.attrs["units"] = "Wm^-2"
        qpar_sheath.attrs["direction_y"] = "Aligned"
        qpar_sheath.attrs["regions"] = {region: ds.regions[region]}
        qpar_sheath.attrs["metadata"] = ds.metadata

        if not return_aligned and "zShift" in self.data.coords:
            qpar_sheath = lower_from_field_aligned(ds["V"].sheath, qpar_sheath)

        return qpar_sheath

    def qpar_upper(self, region=None, *, return_aligned=False):
        """
        Conductive electron heat flux at the upper sheaths.

        Parameters
        ----------
        region : str, optional
            Get the upper sheath value for a specific region. By default looks for all
            regions.
        return_aligned : bool, default False
            Return result on field-aligned grid (note, if input is field-aligned already
            then the result is always on the field-aligned grid)

        Returns
        -------
        DataArray or dict of DataArray
            If called with a 'region' argument or if only one sheath exists, then
            returns a DataArray with the result for that sheath. Otherwise returns a
            dict whose keys are region names and values are result for that region
        """
        if self.data.metadata["storm_normalised"]:
            raise ValueError(
                "qpar_upper() is only implemented for unnormalised Datasets so far"
            )

        if region is None:
            return apply_to_all_regions_upper(
                self.data, self.qpar_upper, return_aligned=return_aligned
            )

        ds = self.data

        echarge = ds.params["e"]
        me = ds.params["m_e"]
        eV_float_over_T = ds.params["dimless"]["eV_float_over_T"]

        T_e = self.T_e_upper(region, return_aligned=True)
        V = self.V_upper(region, return_aligned=True)
        n = self.n_upper(region, return_aligned=True)

        qpar_sheath = (
            ((np.abs(eV_float_over_T) - 0.5) * echarge * T_e - 0.5 * me * V**2)
            * n
            * V
        )

        qpar_sheath.attrs["standard_name"] = "qpar_sheath upper"
        qpar_sheath.attrs[
            "long_name"
        ] = "Conductive electron parallel heat flux at upper sheath entrance"
        qpar_sheath.attrs["units"] = "Wm^-2"
        qpar_sheath.attrs["direction_y"] = "Aligned"
        qpar_sheath.attrs["regions"] = {region: ds.regions[region]}
        qpar_sheath.attrs["metadata"] = ds.metadata

        if not return_aligned and "zShift" in self.data.coords:
            qpar_sheath = upper_from_field_aligned(ds["V"].sheath, qpar_sheath)

        return qpar_sheath

    def electron_total_energy_flux_lower(self, region=None, *, return_aligned=False):
        """
        Total electron parallel energy flux at the lower sheaths. Uses equation (2.14) of
        Helander's book to convert the conductive heat flux, temperature and parallel
        velocity into a total energy flux, neglecting the electron viscosity as it is
        small in mass ratio, and neglecting the perpendicular kinetic energy as small in
        rho_e/L.

        This is the parallel energy flux, which is not the energy flux to the wall unless
        the wall is perpendicular to B.

        Parameters
        ----------
        region : str, optional
            Get the lower sheath value for a specific region. By default looks for all
            regions.
        return_aligned : bool, default False
            Return result on field-aligned grid (note, if input is field-aligned already
            then the result is always on the field-aligned grid)

        Returns
        -------
        DataArray or dict of DataArray
            If called with a 'region' argument or if only one sheath exists, then
            returns a DataArray with the result for that sheath. Otherwise returns a
            dict whose keys are region names and values are result for that region
        """
        if self.data.metadata["storm_normalised"]:
            raise ValueError(
                "electron_total_energy_flux_lower() is only implemented for "
                "unnormalised Datasets so far"
            )

        if region is None:
            return apply_to_all_regions_lower(
                self.data,
                self.electron_total_energy_flux_lower,
                return_aligned=return_aligned,
            )

        ds = self.data

        m_e = ds.params["m_e"]
        echarge = ds.params["e"]

        T_e = self.T_e_lower(region, return_aligned=return_aligned)
        V = self.V_lower(region, return_aligned=return_aligned)
        n = self.n_lower(region, return_aligned=return_aligned)
        qpar = self.qpar_lower(region, return_aligned=return_aligned)

        total_energy_flux = (
            qpar + 2.5 * n * (echarge * T_e) * V + 0.5 * m_e * n * V**3
        )

        total_energy_flux.attrs["standard_name"] = "Qpar_e sheath lower"
        total_energy_flux.attrs[
            "long_name"
        ] = "Total electron parallel energy flux at lower sheath entrance"
        total_energy_flux.attrs["units"] = "Wm^-2"
        total_energy_flux.attrs["direction_y"] = "Aligned"
        total_energy_flux.attrs["regions"] = {region: ds.regions[region]}
        total_energy_flux.attrs["metadata"] = ds.metadata

        return total_energy_flux

    def electron_total_energy_flux_upper(self, region=None, *, return_aligned=False):
        """
        Total electron parallel energy flux at the upper sheaths. Uses equation (2.14) of
        Helander's book to convert the conductive heat flux, temperature and parallel
        velocity into a total energy flux, neglecting the electron viscosity as it is
        small in mass ratio, and neglecting the perpendicular kinetic energy as small in
        rho_e/L.

        This is the parallel energy flux, which is not the energy flux to the wall unless
        the wall is perpendicular to B.

        Parameters
        ----------
        region : str, optional
            Get the upper sheath value for a specific region. By default looks for all
            regions.
        return_aligned : bool, default False
            Return result on field-aligned grid (note, if input is field-aligned already
            then the result is always on the field-aligned grid)

        Returns
        -------
        DataArray or dict of DataArray
            If called with a 'region' argument or if only one sheath exists, then
            returns a DataArray with the result for that sheath. Otherwise returns a
            dict whose keys are region names and values are result for that region
        """
        if self.data.metadata["storm_normalised"]:
            raise ValueError(
                "electron_total_energy_flux_upper() is only implemented for "
                "unnormalised Datasets so far"
            )

        if region is None:
            return apply_to_all_regions_upper(
                self.data,
                self.electron_total_energy_flux_upper,
                return_aligned=return_aligned,
            )

        ds = self.data

        m_e = ds.params["m_e"]
        echarge = ds.params["e"]

        T_e = self.T_e_upper(region, return_aligned=return_aligned)
        V = self.V_upper(region, return_aligned=return_aligned)
        n = self.n_upper(region, return_aligned=return_aligned)
        qpar = self.qpar_upper(region, return_aligned=return_aligned)

        total_energy_flux = (
            qpar + 2.5 * n * (echarge * T_e) * V + 0.5 * m_e * n * V**3
        )

        total_energy_flux.attrs["standard_name"] = "Qpar_e sheath upper"
        total_energy_flux.attrs[
            "long_name"
        ] = "Total electron parallel energy flux at upper sheath entrance"
        total_energy_flux.attrs["units"] = "Wm^-2"
        total_energy_flux.attrs["direction_y"] = "Aligned"
        total_energy_flux.attrs["regions"] = {region: ds.regions[region]}
        total_energy_flux.attrs["metadata"] = ds.metadata

        return total_energy_flux

    def metric_lower(self, component, region=None):
        """
        Get a metric component at the location of the lower sheath entrance

        Parameters
        ----------
        component : str
            Which metric component to get.
        region : str, optional
            Get the upper sheath value for a specific region. By default looks for all
            regions.
        """

        if region is None:
            return apply_to_all_regions_lower(
                self.data,
                self.metric_lower,
            )

        if self.data.regions[region].connection_lower_y is not None:
            raise ValueError(
                f"Requested lower sheath value for '{region}', but '{region}' has no "
                f"lower parallel boundary."
            )

        da = self.data[f"{component}_CELL_YLOW"].bout.from_region(region, with_guards=0)

        if self.metadata["keep_yboundaries"]:
            yind = self.metadata["MYG"]
        else:
            yind = 0

        return da.isel({da.metadata["bout_ydim"]: yind})

    def metric_upper(self, component, region=None):
        """
        Get a metric component at the location of the upper sheath entrance

        Parameters
        ----------
        component : str
            Which metric component to get.
        region : str, optional
            Get the upper sheath value for a specific region. By default looks for all
            regions.
        """
        if not self.data.metadata["keep_yboundaries"]:
            raise ValueError(
                f"keep_yboundaries=False but y-boundary cells are required to get the "
                f"value of {component} at the upper sheath"
            )
        if self.data.metadata["MYG"] == 0:
            raise ValueError(
                f"MYG=0 but y-boundary cells are required to get the value of "
                f"{component} at the upper sheath"
            )

        if region is None:
            return apply_to_all_regions_upper(
                self.data,
                self.metric_upper,
            )

        if self.data.regions[region].connection_upper_y is not None:
            raise ValueError(
                f"Requested upper sheath value for '{region}', but '{region}' has no "
                f"upper parallel boundary."
            )

        da = self.data[f"{component}_CELL_YLOW"].bout.from_region(region, with_guards=0)

        return da.isel({da.metadata["bout_ydim"]: -self.metadata["MYG"]})

    def normal_electron_particle_flux_lower(
        self, region=None, *, include_y_derivs=None
    ):
        """
        Electron particle flux projected normal to the wall, i.e. electron flux
        per unit area actually arriving on the wall.

        This method calculates the flux in the y-direction (not the outward
        flux) in order to be consistent with the sign of V_e,parallel.

        Parameters
        ----------
        region : str, optional
            Get the upper sheath value for a specific region. By default looks for all
            regions.
        include_y_derivs : bool, optional
            Include parallel derivatives in calculation of normal flux. In this
            case means including poloidal components of ExB flux and diffusive
            flux. These are often neglected as small. Default is False to match
            STORM.
        """
        if include_y_derivs:
            raise ValueError("Not implemented yet")

        if region is None:
            return apply_to_all_regions_lower(
                self.data,
                self.normal_electron_particle_flux_lower,
                include_y_derivs=include_y_derivs,
            )

        ds = self.data

        # Use same projection factor as in StormDataset.v_par_e() method
        g11 = self.metric_lower("g11", region)
        g_22 = self.metric_lower("g_22", region)
        g_23 = self.metric_lower("g_23", region)
        g_33 = self.metric_lower("g_33", region)
        J = self.metric_lower("J", region)
        result = (
            self.n_lower(region)
            * self.V_lower(region)
            * (-(g_23**2) + g_33 * g_22)
            / (J * np.sqrt(g11 * g_22 * g_33))
        )

        if include_y_derivs:
            # Todo: Add ExB and density diffusion contributions here??
            pass

        result.name = f"lower sheath electron particle flux"
        result.attrs["standard_name"] = result.name
        result.attrs["long_name"] = result.name
        if result.metadata["storm_normalised"]:
            result.attrs["units"] = "Omega_i0 rho_s0^-2"
        else:
            result.attrs["units"] = "s^-1 m^-2"

        return result

    def normal_electron_particle_flux_upper(
        self, region=None, *, include_y_derivs=None
    ):
        """
        Electron particle flux projected normal to the wall, i.e. electron flux
        per unit area actually arriving on the wall.

        This method calculates the flux in the y-direction (not the outward
        flux) in order to be consistent with the sign of V_e,parallel.

        Parameters
        ----------
        region : str, optional
            Get the upper sheath value for a specific region. By default looks for all
            regions.
        include_y_derivs : bool, optional
            Include parallel derivatives in calculation of normal flux. In this
            case means including poloidal components of ExB flux and diffusive
            flux. These are often neglected as small. Default is False to match
            STORM.
        """
        if include_y_derivs:
            raise ValueError("Not implemented yet")

        if region is None:
            return apply_to_all_regions_upper(
                self.data,
                self.normal_electron_particle_flux_upper,
                include_y_derivs=include_y_derivs,
            )

        ds = self.data

        # Use same projection factor as in StormDataset.v_par_e() method
        g11 = self.metric_upper("g11", region)
        g_22 = self.metric_upper("g_22", region)
        g_23 = self.metric_upper("g_23", region)
        g_33 = self.metric_upper("g_33", region)
        J = self.metric_upper("J", region)
        result = (
            self.n_upper(region)
            * self.V_upper(region)
            * (-(g_23**2) + g_33 * g_22)
            / (J * np.sqrt(g11 * g_22 * g_33))
        )

        if include_y_derivs:
            # Todo: Add ExB and density diffusion contributions here??
            pass

        result.name = f"upper sheath electron particle flux"
        result.attrs["standard_name"] = result.name
        result.attrs["long_name"] = result.name
        if result.metadata["storm_normalised"]:
            result.attrs["units"] = "Omega_i0 rho_s0^-2"
        else:
            result.attrs["units"] = "s^-1 m^-2"

        return result


@register_dataarray_accessor("sheath")
class SheathDataArrayAccessor(StormDataArrayAccessor):
    """
    Calculate values at sheath boundaries, replicating algorithms used in STORM.

    Implements methods for analysing single variables.
    """

    def __init__(self, da):
        super().__init__(da)

    def extrapolate_to_lower(
        self, region=None, *, ylow_force_extrapolate=False, return_aligned=False
    ):
        """
        Extrapolates the BoutDataArray in the parallel direction to the position of the
        sheath entrance at the lower boundary.

        Extrapolation like STORM's extrap_sheath_upper(var), with boundary cell value
        set by free_o3 boundary condition (this method does not require boundary cell
        values to be set).

        Parameters
        ----------
        region : str, optional
            Get the lower sheath value for a specific region. By default looks for all
            regions.
        ylow_force_extrapolate : bool, default False
            When getting sheath values for a CELL_YLOW field, extrapolate from points
            inside the domain instead of using the values on the boundary. Used for
            allowing supersonic flow in ion velocity boundary condition.
        return_aligned : bool, default False
            Return result on field-aligned grid (note, if input is field-aligned already
            then the result is always on the field-aligned grid)

        Returns
        -------
        DataArray or dict of DataArray
            If called with a 'region' argument or if only one sheath exists, then
            returns a DataArray with the result for that sheath. Otherwise returns a
            dict whose keys are region names and values are result for that region
        """
        if ylow_force_extrapolate and self.data.cell_location != "CELL_YLOW":
            raise ValueError(
                f"ylow_force_extrapolate option can only be used with CELL_YLOW "
                f"fields, but {self.data.name} is at {self.data.cell_location}."
            )

        if region is None:
            return apply_to_all_regions_lower(
                self.data,
                self.extrapolate_to_lower,
                ylow_force_extrapolate=ylow_force_extrapolate,
                return_aligned=return_aligned,
            )

        if self.data.regions[region].connection_lower_y is not None:
            raise ValueError(
                f"Requested lower sheath value for '{region}', but '{region}' has no "
                f"lower parallel boundary."
            )

        ycoord = self.data.metadata["bout_ydim"]
        zcoord = self.data.metadata["bout_zdim"]

        da = self.data.bout.from_region(region, with_guards=0)

        if da.cell_location == "CELL_YLOW" and not ylow_force_extrapolate:
            # Already have data on boundary location, no need to extrapolate
            if self.metadata["keep_yboundaries"]:
                yind = self.metadata["MYG"]
            else:
                yind = 0
            return da.isel({ycoord: yind})

        if (
            "zShift" in self.data.coords
            and zcoord in da.dims
            and da.direction_y != "Aligned"
        ):
            aligned_input = False
            da = da.bout.to_field_aligned()
        else:
            aligned_input = True

        if self.metadata["keep_yboundaries"]:
            # Do not use boundary cells
            da = da.isel({ycoord: slice(self.metadata["MYG"], None)})

        def extrap_to_sheath(f, *, ylow_force_extrapolate=False):
            if ylow_force_extrapolate:
                # values near boundary
                y0 = f.isel({ycoord: 1})
                y1 = f.isel({ycoord: 2})
                y2 = f.isel({ycoord: 3})

                return 3.0 * y0 - 3.0 * y1 + y2
            else:
                # values near boundary
                y0 = f.isel({ycoord: 0})
                y1 = f.isel({ycoord: 1})
                y2 = f.isel({ycoord: 2})

                # Extrapolate like free_o3 boundary condition
                boundary_cell_value = 3.0 * y0 - 3.0 * y1 + y2

                # Interpolate like extrap_sheath_lower(var)
                return 0.375 * boundary_cell_value + 0.75 * y0 - 0.125 * y1

        sheath_slice = extrap_to_sheath(
            da, ylow_force_extrapolate=ylow_force_extrapolate
        )

        if not aligned_input and not return_aligned:
            # Transform back to non-aligned coordinates
            sheath_slice = lower_from_field_aligned(self, sheath_slice)

        return sheath_slice

    def extrapolate_to_upper(
        self, region=None, *, ylow_force_extrapolate=False, return_aligned=False
    ):
        """
        Extrapolates the BoutDataArray in the parallel direction to the position of the
        sheath entrance at the upper boundary.

        Extrapolation like STORM's extrap_sheath_upper(var), with boundary cell value
        set by free_o3 boundary condition (this method does not require boundary cell
        values to be set).

        Parameters
        ----------
        region : str, optional
            Get the upper sheath value for a specific region. By default looks for all
            regions.
        ylow_force_extrapolate : bool, default False
            When getting sheath values for a CELL_YLOW field, extrapolate from points
            inside the domain instead of using the values on the boundary. Used for
            allowing supersonic flow in ion velocity boundary condition.
        return_aligned : bool, default False
            Return result on field-aligned grid (note, if input is field-aligned already
            then the result is always on the field-aligned grid)

        Returns
        -------
        DataArray or dict of DataArray
            If called with a 'region' argument or if only one sheath exists, then
            returns a DataArray with the result for that sheath. Otherwise returns a
            dict whose keys are region names and values are result for that region
        """
        if ylow_force_extrapolate and self.data.cell_location != "CELL_YLOW":
            raise ValueError(
                f"ylow_force_extrapolate option can only be used with CELL_YLOW "
                f"fields, but {self.data.name} is at {self.data.cell_location}."
            )

        if region is None:
            return apply_to_all_regions_upper(
                self.data,
                self.extrapolate_to_upper,
                ylow_force_extrapolate=ylow_force_extrapolate,
                return_aligned=return_aligned,
            )

        if self.data.regions[region].connection_upper_y is not None:
            raise ValueError(
                f"Requested upper sheath value for '{region}', but '{region}' has no "
                f"upper parallel boundary."
            )

        ycoord = self.data.metadata["bout_ydim"]
        zcoord = self.data.metadata["bout_zdim"]

        da = self.data.bout.from_region(region, with_guards=0)

        if da.cell_location == "CELL_YLOW" and not ylow_force_extrapolate:
            # Already have data on boundary location, no need to extrapolate
            if self.metadata["keep_yboundaries"] and self.metadata["MYG"] > 0:
                yind = -self.metadata["MYG"]
            else:
                raise ValueError(
                    "Cannot get sheath value from CELL_YLOW field without boundary "
                    "cells"
                )
            return da.isel({ycoord: yind})

        if (
            "zShift" in self.data.coords
            and zcoord in da.dims
            and da.direction_y != "Aligned"
        ):
            aligned_input = False
            da = da.bout.to_field_aligned()
        else:
            aligned_input = True

        myg = self.metadata["MYG"]
        if self.metadata["keep_yboundaries"] and myg > 0:
            # Do not use boundary cells
            da = da.isel({ycoord: slice(None, -myg)})

        def extrap_to_sheath(f, *, ylow_force_extrapolate=False):
            if ylow_force_extrapolate:
                # values near boundary
                y0 = f.isel({ycoord: -1})
                y1 = f.isel({ycoord: -2})
                y2 = f.isel({ycoord: -3})

                return 3.0 * y0 - 3.0 * y1 + y2
            else:
                # values near boundary
                y0 = f.isel({ycoord: -1})
                y1 = f.isel({ycoord: -2})
                y2 = f.isel({ycoord: -3})

                # Extrapolate like free_o3 boundary condition
                boundary_cell_value = 3.0 * y0 - 3.0 * y1 + y2

                # Interpolate like extrap_sheath_upper(var)
                return 0.375 * boundary_cell_value + 0.75 * y0 - 0.125 * y1

        sheath_slice = extrap_to_sheath(
            da, ylow_force_extrapolate=ylow_force_extrapolate
        )

        if not aligned_input and not return_aligned:
            # Transform back to non-aligned coordinates
            sheath_slice = upper_from_field_aligned(self, sheath_slice)

        return sheath_slice
