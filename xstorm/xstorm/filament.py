import numpy as np
import xarray as xr

from xbout import BoutDatasetAccessor

from .sheathutils import (
    apply_to_all_regions_lower,
    apply_to_all_regions_upper,
)


@xr.register_dataset_accessor("filament")
class FilamentDatasetAccessor(BoutDatasetAccessor):
    """
    Filament-specific analysis methods intended to be used with data obtained
    from a simulation using the STORM module for BOUT++.

    This accessor assumes the simulation was in xstorm's 'slab' coordinates in
    various places.
    """

    def __init__(self, ds):
        super().__init__(ds)

    @property
    def _spatial_dims(self):
        return list(set(self.data.dims) - set(["time"]))

    def _get_background_subtracted_fluctuation(
        self, name, return_mass=False, dims=None
    ):
        ds = self.data
        da = ds[name]

        if dims is None:
            dims = self._spatial_dims

        # Subtract background from parallel equilibrium, if one exists for this variable
        eq_name = name + "_eq"
        if eq_name in ds:
            da = da - ds[eq_name]
            da.name = name

        # Only include positive part of fluctuations, usually gives better measure of
        # filament position
        # The condition indicates where to keep values from da, the second argument the
        # value to put where the condition is false
        # da = da.where(da > 0.0, 0.0)

        if not return_mass:
            return da
        else:
            total_mass = da.sum(dim=dims)
            return da, total_mass

    def CoM(self, name, *, dims=None):
        """
        Calculate the 'centre-of-mass' as a function of time for the variable called
        'name'. Acts on as many spatial dimensions as the Dataset has. E.g. can be used
        on full 3d dataset, or on a slice at the midplane, etc.

        Parameters
        ----------
        name : str
            The variable to calculate the CoM for.
        dims : str or sequence of str, optional
            Dimensions to integrate over to calculate the CoM. By default uses all
            spatial dimensions in the Dataset.

        Returns
        -------
        dict - entry for each spatial dimension in the Dataset contains 'centre-of-mass'
               in that direction
        """
        ds = self.data
        da, total_mass = self._get_background_subtracted_fluctuation(
            name, return_mass=True, dims=dims
        )

        if dims is None:
            dims = self._spatial_dims

        result = {}
        for d in self._spatial_dims:
            coord = ds[d]
            CoM = (coord * da).sum(dim=dims) / total_mass
            CoM.attrs["units"] = coord.attrs["units"]
            CoM.attrs["standard_name"] = f"{d} centre of {name}"
            CoM.attrs["long_name"] = f"Centre of {name} in {d} direction"
            result[d] = CoM
            ds[f"{name}_CoM_{d}"] = CoM

        return result

    def CoM_velocity(self, name, direction, *, method=None, dims=None):
        """
        Centre-of-mass velocity for variable 'name' calculated from centre-of-mass
        position using finite-difference in time

        Parameters
        ----------
        name : str
            Name of the variable to calculate CoM velocity for
        direction : str
            Component of the CoM velocity to calculate. {"radial", "parallel", "binormal"}
        method : str, default "E-field"
            Method to use to calculate the velocity:
              - "finite-difference" uses finite difference in time of the 'direction'
                coordinate of the centre-of-mass
              - "E-field" uses the ExB velocity for perpendicular velocities and electron
                parallel velocity V for parallel velocity to calculate the instantaneous
                'name'-weighted velocity
        dims : str or sequence of str, optional
            Dimensions to integrate over to calculate the CoM velocity. By default uses
            all spatial dimensions in the Dataset.
        """
        ds = self.data

        if dims is None:
            dims = self._spatial_dims

        if method is None:
            # Set default value
            method = "E-field"

        if method == "finite-difference":
            CoM_name = f"{name}_CoM_{direction}"
            if not CoM_name in ds:
                ds.filament.CoM(name, dims=dims)

            CoM = ds[CoM_name]

            v = CoM.differentiate("time", edge_order=2)

            v.attrs["units"] = CoM.attrs["units"] + ds["time"].attrs["units"] + "-1"
        elif method == "E-field":
            da, total_mass = self._get_background_subtracted_fluctuation(
                name, return_mass=True, dims=dims
            )
            if direction == "radial":
                v_3d = ds.storm.v_ExB_slab_radial()
            elif direction == "binormal":
                v_3d = ds.storm.v_ExB_slab_binormal()
            elif direction == "parallel":
                v_3d = ds["V"]
            else:
                raise ValueError(
                    f"Unrecognised option direction='{direction}' for CoM_velocity()"
                )

            v = (da * v_3d).sum(dim=dims) / total_mass
            v.attrs["units"] = v_3d.attrs["units"]
        else:
            raise ValueError(
                f"Unrecognised option method='{method}' for CoM_velocity()"
            )

        v.attrs["standard_name"] = f"v {direction}"
        v.attrs["long_name"] = f"{direction} centre of mass velocity for {name}"

        return v

    def with_origin_at_initial_filament_position(self):
        """
        Set the origin of the "radial", "parallel" and "binormal" coordinates
        to be the centre of the filament's initial density perturbation.
        """

        if self.data.metadata["storm_normalised"]:
            raise ValueError(
                "with_origin_at_initial_filament_position() is only implemented for "
                "unnormalised Datasets so far"
            )

        ds = self.data.copy()

        def get_with_default(name, default):
            if name in ds.options["blob"]:
                return ds.options["blob"][name]
            else:
                return default

        xoffset = get_with_default("xoffset", 0.25)
        yoffset = get_with_default("yoffset", 0.5)
        zoffset = get_with_default("zoffset", 0.5)

        def get_L(name):
            try:
                return ds.options["filaments"].evaluate_scalar(name)
            except KeyError:
                try:
                    return ds.options["mesh"].evaluate_scalar(name)
                except KeyError:
                    raise ValueError(
                        f"with_origin_at_initial_filament_position() requires '{name}' "
                        f"to be in the options, but was not found in 'filaments' or "
                        f"'mesh' sections"
                    )

        Lx = get_L("Lx")
        Ly = get_L("Ly")
        Lz = get_L("Lz")

        radial_attrs = ds["radial"].attrs
        ds["radial"] = ds["radial"] - xoffset * Lx * ds.params["rho_s0"]
        ds["radial"].attrs = radial_attrs

        parallel_attrs = ds["parallel"].attrs
        ds["parallel"] = ds["parallel"] - (yoffset - 0.5) * Ly * ds.params["rho_s0"]
        ds["parallel"].attrs = parallel_attrs

        binormal_attrs = ds["binormal"].attrs
        ds["binormal"] = ds["binormal"] - zoffset * Lz * ds.params["rho_s0"]
        ds["binormal"].attrs = binormal_attrs

        return ds

    def time_index_of_max_CoM_velocity(
        self, variable, *, direction="radial", method=None
    ):
        """
        Find the index in the time dimension where filament reaches maximum
        centre-of-mass velocity.

        Parameters
        ----------
        See CoM_velocity
        """
        v = self.data.filament.CoM_velocity(
            variable, direction=direction, method=method
        )
        result = v.argmax(dim=...)

        # workaround for current bug in xarray
        # (https://github.com/pydata/xarray/issues/4276) if ds is using dask
        result["time"] = result["time"].load()

        return result

    def select_at_max_CoM_velocity(self, variable, *, direction="radial", method=None):
        """
        Find the index in the time dimension where filament reaches maximum
        centre-of-mass velocity.

        Parameters
        ----------
        See CoM_velocity
        """
        timeind = self.time_index_of_max_CoM_velocity(
            variable, direction=direction, method=method
        )

        return self.data.isel(timeind)

    def interpolate_relative_to_CoM(
        self,
        variable,
        *,
        radial_width=0.04,
        parallel_width=None,
        binormal_width=0.04,
        n_radial=32,
        n_parallel=32,
        n_binormal=32,
    ):
        """
        Interpolate a Dataset to a grid of coordinates relative to the filament
        centre-of-mass

        Parameters
        ----------
        variable : str
            Name of the variable to calculate the centre-of-mass from.
        radial_width : float or None
            Width of the radial grid relative to the centre-of-mass. If None,
            do not interpolate in the radial direction.
        parallel_width : float or None
            Width of the parallel grid relative to the centre-of-mass. If None,
            do not interpolate in the parallel direction.
        binormal_width : float or None
            Width of the binormal grid relative to the centre-of-mass. If None,
            do not interpolate in the binormal direction.
        n_radial : int
            Number of points in the radial interpolated grid, if one exists.
        n_parallel : int
            Number of points in the parallel interpolated grid, if one exists.
        n_binormal : int
            Number of points in the binormal interpolated grid, if one exists.

        Returns
        -------
        Dateset interpolated onto coordinates relative to the filament
        centre-of-mass
        """
        ds = self.data
        ds.filament.CoM(variable)

        interp_coords = {}
        if radial_width is not None:
            grid_r = xr.DataArray(
                np.linspace(
                    -radial_width / 2.0 + radial_width / n_radial / 2.0,
                    radial_width / 2.0 - radial_width / n_radial / 2.0,
                    n_radial,
                ),
                dims="relative_radial",
            )
            interp_coords["radial"] = ds[variable + "_CoM_radial"] + grid_r
        if parallel_width is not None:
            grid_p = xr.DataArray(
                np.linspace(
                    -parallel_width / 2.0 + parallel_width / n_parallel / 2.0,
                    parallel_width / 2.0 - parallel_width / n_parallel / 2.0,
                    n_parallel,
                ),
                dims="relative_parallel",
            )
            interp_coords["parallel"] = ds[variable + "_CoM_parallel"] + grid_p
        if binormal_width is not None:
            grid_b = xr.DataArray(
                np.linspace(
                    -binormal_width / 2.0 + binormal_width / n_binormal / 2.0,
                    binormal_width / 2.0 - binormal_width / n_binormal / 2.0,
                    n_binormal,
                ),
                dims="relative_binormal",
            )
            interp_coords["binormal"] = ds[variable + "_CoM_binormal"] + grid_b

        return ds.interp(interp_coords)

    def integrate_over_subregion(
        self, variable, *, radial_sel=0.05, parallel_sel=None, binormal_sel=None
    ):
        """
        Volume-integrate the amount of 'variable' with background subtracted,
        in a sub-region specified by the keyword arguments.

        Parameters
        ----------
        variable : str
            Name of the variable to integrate
        radial_sel : float, slice or None
            Radial coordinate selection, in units of m.
                - float - select the region at larger radius than this value.
                - slice - select the region specified by the slice. For example to
                  select the region inside a given value, pass slice(value).
                - None - keep the whole radial region.
        parallel_sel : float, slice or None
            Parallel coordinate selection, in units of m.
                 - float - select the region at larger parallel coordinate than this value.
                 - slice - select the region specified by the slice.
                 - None - keep the whole parallel region.
        binormal_sel : float, slice or None
            Binormal coordinate selection, in units of m.
                 - float - select the region at larger binormal coordinate than this value.
                 - slice - select the region specified by the slice.
                 - None - keep the whole binormal region.
        """

        f = self._get_background_subtracted_fluctuation(variable)

        if isinstance(radial_sel, float):
            radial_sel = slice(radial_sel, None)
        elif radial_sel is None:
            radial_sel = slice(None)
        elif not isinstance(radial_sel, slice):
            raise ValueError(
                f"radial_sel should be float, slice or None. Got {radial_sel}."
            )

        if isinstance(parallel_sel, float):
            parallel_sel = slice(parallel_sel, None)
        elif parallel_sel is None:
            parallel_sel = slice(None)
        elif not isinstance(parallel_sel, slice):
            raise ValueError(
                f"parallel_sel should be float, slice or None. Got {parallel_sel}."
            )

        if isinstance(binormal_sel, float):
            binormal_sel = slice(binormal_sel, None)
        elif binormal_sel is None:
            binormal_sel = slice(None)
        elif not isinstance(binormal_sel, slice):
            raise ValueError(
                f"binormal_sel should be float, slice or None. Got {binormal_sel}."
            )

        # select by coordinate value, not index value
        f = f.sel(radial=radial_sel, parallel=parallel_sel, binormal=binormal_sel)

        total = self.data.bout.integrate_midpoints(f)

        if variable == "n":
            total.attrs["units"] = "electrons"
        elif "units" in f.attrs:
            total.attrs["units"] = f.attrs["units"] + "m^3"

        return total

    def normal_electron_particle_flux_lower_background_subtracted(
        self, region=None, *, include_y_derivs=None
    ):
        """
        Electron particle flux only from the density fluctuation above the
        background, projected normal to the wall, i.e. electron flux per unit
        area actually arriving on the wall.

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
                self.normal_electron_particle_flux_lower_background_subtracted,
                include_y_derivs=include_y_derivs,
            )

        ds = self.data

        # Use same projection factor as in StormDataset.v_par_e() method
        g11 = ds.sheath.metric_lower("g11", region)
        g_22 = ds.sheath.metric_lower("g_22", region)
        g_23 = ds.sheath.metric_lower("g_23", region)
        g_33 = ds.sheath.metric_lower("g_33", region)
        J = ds.sheath.metric_lower("J", region)
        result = (
            self._get_background_subtracted_fluctuation(
                "n"
            ).sheath.extrapolate_to_lower(region)
            * ds.sheath.V_lower(region)
            * (-(g_23**2) + g_33 * g_22)
            / (J * np.sqrt(g11 * g_22 * g_33))
        )

        if include_y_derivs:
            # Todo: Add ExB and density diffusion contributions here??
            pass

        result.name = f"lower sheath electron density fluctuation particle flux"
        result.attrs["standard_name"] = result.name
        result.attrs["long_name"] = result.name
        if result.metadata["storm_normalised"]:
            result.attrs["units"] = "Omega_i0 rho_s0^-2"
        else:
            result.attrs["units"] = "s^-1 m^-2"

        return result

    def normal_electron_particle_flux_upper_background_subtracted(
        self, region=None, *, include_y_derivs=None
    ):
        """
        Electron particle flux only from the density fluctuation above the
        background, projected normal to the wall, i.e. electron flux per unit
        area actually arriving on the wall.

        This method calculates the flux in the y-direction (which is the
        outward flux) in order to be consistent with the sign of V_e,parallel.

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
                self.normal_electron_particle_flux_upper_background_subtracted,
                include_y_derivs=include_y_derivs,
            )

        ds = self.data

        # Use same projection factor as in StormDataset.v_par_e() method
        g11 = ds.sheath.metric_upper("g11", region)
        g_22 = ds.sheath.metric_upper("g_22", region)
        g_23 = ds.sheath.metric_upper("g_23", region)
        g_33 = ds.sheath.metric_upper("g_33", region)
        J = ds.sheath.metric_upper("J", region)
        result = (
            self._get_background_subtracted_fluctuation(
                "n"
            ).sheath.extrapolate_to_upper(region)
            * ds.sheath.V_upper(region)
            * (-(g_23**2) + g_33 * g_22)
            / (J * np.sqrt(g11 * g_22 * g_33))
        )

        if include_y_derivs:
            # Todo: Add ExB and density diffusion contributions here??
            pass

        result.name = f"upper sheath electron density fluctuation particle flux"
        result.attrs["standard_name"] = result.name
        result.attrs["long_name"] = result.name
        if result.metadata["storm_normalised"]:
            result.attrs["units"] = "Omega_i0 rho_s0^-2"
        else:
            result.attrs["units"] = "s^-1 m^-2"

        return result
