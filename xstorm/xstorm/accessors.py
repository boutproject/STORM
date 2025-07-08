from xarray import register_dataset_accessor, register_dataarray_accessor, apply_ufunc
import xarray as xr

from itertools import chain
from matplotlib import pyplot as plt
import numpy as np
from warnings import warn
import dask.array.stats as daskstats

from xbout import BoutDatasetAccessor, BoutDataArrayAccessor

from .load import calc_norms
from .timeseries import conditional_average, waiting_times
from .utils import _update_name, _update_units


@register_dataset_accessor("storm")
class StormDatasetAccessor(BoutDatasetAccessor):
    """
    Class specifically for accessing analysis methods which are to be
    performed on data obtained from a simulation using the STORM module
    for BOUT++.

    These methods are on a Dataset because they implement analysis which
    require multiple variables.
    """

    def __init__(self, ds):
        super().__init__(ds)

    def norms(self):
        return calc_norms(self.data)

    def save(self, *args, **kwargs):
        """
        Wrapper for bout.save() that handles STORM-specific attributes
        """
        # Shallow copy to ensure we do not change the attrs of the Dataset
        to_save = self.data.copy()

        # Delete 'settings' and 'params' as we cannot save a dict to netCDF and they will
        # be re-created when re-loading with open_stormdataset anyway.
        # When xBOUT implements a better way of handling and saving 'options', might be
        # worth copying that here for 'settings'
        del to_save.attrs["settings"]
        del to_save.attrs["params"]
        # Delete settings on each variable and coordinate
        for var in chain(to_save.data_vars, to_save.coords):
            try:
                del to_save[var].attrs["settings"]
            except KeyError:
                pass
            try:
                del to_save[var].attrs["params"]
            except KeyError:
                pass

        return to_save.bout.save(*args, **kwargs)

    def v_ExB_slab_radial(self):
        """Calculates local radial ExB velocity"""

        if "v_radial" not in self.data:
            E_z = self.data["phi"].differentiate(coord="binormal")
            v_radial = E_z / self.data["B"]
            v_radial.attrs["units"] = "ms-1"
            v_radial.attrs["standard_name"] = "radial velocity"
            v_radial.attrs["long_name"] = "Radial ExB Velocity"
            self.data["v_radial"] = v_radial
        return self.data["v_radial"]

    def v_ExB_slab_binormal(self):
        """Calculates local binormal ExB velocity"""

        if "v_binormal" not in self.data:
            E_x = self.data["phi"].differentiate(coord="radial")
            v_binormal = -E_x / self.data["B"]
            v_binormal.attrs["units"] = "ms-1"
            v_binormal.attrs["standard_name"] = "binormal velocity"
            v_binormal.attrs["long_name"] = "Binormal ExB Velocity"
            self.data["v_binormal"] = v_binormal
        return self.data["v_binormal"]

    def flux_slab(self, moment="particle", direction="radial"):
        """Calculates local flux of particles or energy"""

        if moment == "particle":
            return self.data.storm.ion_particle_flux_slab(direction)
        elif moment == "energy":
            return self.data.storm.energy_flux_slab(direction)
        else:
            raise ValueError

    def ion_particle_flux_slab(self, direction="radial"):
        """Calculates local particle flux"""

        warn(
            "ion_particle_flux_slab only includes ExB drift in perpendicular fluxes "
            "(not magnetic drift or polarisation drift)"
        )

        if direction == "radial":
            flux = self.data["n"] * self.data.storm.v_ExB_slab_radial()
            flux.attrs["units"] = "m-2s-1"
            flux.attrs["standard_name"] = "radial ion particle flux"
            flux.attrs["long_name"] = "Radial Ion Particle Flux"
            self.data["flux_radial"] = flux
        elif direction == "binormal":
            flux = self.data["n"] * self.data.storm.v_ExB_slab_binormal()
            flux.attrs["units"] = "m-2s-1"
            flux.attrs["standard_name"] = "binormal ion particle flux"
            flux.attrs["long_name"] = "Binormal Ion Particle Flux"
            self.data["flux_binormal"] = flux
        elif direction == "parallel":
            flux = self.data["n"] * self.data["U"]
            flux.attrs["units"] = "m-2s-1"
            flux.attrs["standard_name"] = "parallel ion particle flux"
            flux.attrs["long_name"] = "Parallel Ion Particle Flux"
            self.data["flux_parallel"] = flux
        else:
            raise ValueError
        return flux

    def energy_flux(self, dir="radial"):
        # TODO should this include thermal energy?
        raise NotImplementedError

    def v_ExB(self, direction=None, *, include_y_derivs=None):
        """
        Calculates ExB velocity

        The radial component is `v_ExB.Grad(psi)/|Grad(psi)|`.
        The poloidal component is `v_ExB.Grad(x) x e_z / |Grad(x)| |e_z|`.
        The toroidal component is `v_ExB.e_z/|e_z|`.

        Parameters
        ----------
        direction : str, optional
            By default, calculates radial, poloidal and toroidal components and returns
            them as a dict. Can pass "radial", "poloidal" or "toroidal" to return a
            specific component.
        include_y_derivs : bool, optional
            Include parallel derivatives in calculation of v_ExB. These are often
            neglected as small. Default is False to match STORM.
        """
        if direction is None:
            return {
                d: self.v_ExB(d, include_y_derivs=include_y_derivs)
                for d in ["radial", "poloidal", "toroidal"]
            }

        if include_y_derivs is None:
            # Default to STORM default behaviour
            include_y_derivs = False

        ds = self.data
        zcoord = ds.metadata["bout_zdim"]

        # From 'Field-aligned coordinates' section of BOUT++ manual
        # (https://bout-dev.readthedocs.io/en/latest/user_docs/coordinates.html)
        # v_ExB = 1/g_yy (- g_yy e_z d/dx + g_yz e_y d/dx + g_xy e_z d/dy
        #                 - g_yz e_x d/dy - g_xy e_y d/dz + g_yy e_x d/dz) phi
        # where g_ij are the covariant components of the metric tensor and e_i are the
        # tangent-basis vectors.
        # When include_yderivs is False, we need to neglect the y-component of v_ExB to
        # be consistent with the way advection is implemented in STORM. In that case, we
        # use the following 'approximate' form to calculate the components.
        # v_ExB = 1/g_yy (- g_yy e_z d/dx + g_xy e_z d/dy
        #                 - g_yz e_x d/dy + g_yy e_x d/dz) phi
        # Note that this does give a non-zero *poloidal* component.

        # The representation above relies on the coordinates being 'Clebsch
        # coordinates' where Grad(z) x Grad(x) = B. For slab simulations STORM does not
        # use Clebsch coordinates. Instead the metric is the identity so the above
        # expression reduces to (-e_z d/dx + e_x d/dz) phi where e_z and e_x are unit
        # vectors. To get the ExB velocity b x Grad(phi) / B for this slab geometry, we
        # just need to divide by B.

        phi = ds["phi"]
        g_22 = ds["g_22"]
        if direction == "radial":
            # x := psi
            # |Grad(x)| = sqrt(gxx)
            # Grad(x).e_i = delta^x_i
            # v_ExB_radial = 1/sqrt(gxx) Grad(x).v_ExB
            #              = 1/(g_yy sqrt(gxx)) (0 + 0 + 0
            #                                    - g_yz d/dy - 0 + g_yy d/dz) phi
            if zcoord in phi.dims:
                result = g_22 * phi.bout.ddz()
            else:
                result = 0.0
            if include_y_derivs:
                result = result - ds["g_23"] * phi.bout.ddy()
            result = result / (g_22 * np.sqrt(ds["g11"]))

            if ds.settings["realistic_geometry"] in ["none", "slab"]:
                # divide by B - see comment above
                result = result / ds["B"]

            result.attrs["standard_name"] = "v_ExB_radial"
            result.attrs["long_name"] = "Radial component of ExB velocity"
            result.attrs["units"] = ds["U"].units

            return result

        elif direction == "poloidal":
            if not include_y_derivs:
                # for non-orthogonal grids, drops all e_y and ddy()
                #  v_ExB_poloidal = 1 / (J g_yy sqrt(gxx g_zz))
                #                   * (
                #                       (- g_yy g_yz g_xz + g_yy g_zz g_xy) d/dz
                #                     ) phi
                result = (
                    phi.bout.ddz()
                    * (ds["g_12"] * ds["g_33"] - ds["g_23"] * ds["g_13"])
                    / (ds["J"] * np.sqrt(ds["g11"] * ds["g_33"]))
                )

                return result

            # Need a a unit vector hat{e}_pol in the poloidal direction. e_z is
            # in the toroidal direction and Grad(x) is orthogonal to flux
            # surfaces, so their cross product is in the poloidal direction
            # (within flux surfaces). e_z and Grad(x) are also always
            # orthogonal, so the magnitude of their cross product is the
            # product of their magnitudes. Therefore
            #   hat{e}_pol = (e_z x Grad(x))/|Grad(x)||e_z|
            #   hat{e}_pol = (e_z x Grad(x))/sqrt(gxx g_zz)
            #
            #   hat{e}_pol . e_i
            #    = (e_z x Grad(x)).e_i/sqrt(gxx g_zz)
            #    = (g_xz Grad(x) + g_yz Grad(y) + g_zz Grad(z))
            #       x Grad(x) . e_i / sqrt(gxx g_zz)
            # Using Grad(x^i) x Grad(x^j) = e_k/J  [i,j,k cyclic combination of x,y,z]
            #   hat{e}_pol . e_i
            #    = (- g_yz e_z + g_zz e_y) . e_i / (J sqrt(gxx g_zz))
            #    = (- g_yz g_iz + g_zz g_iy) / (J sqrt(gxx g_zz))
            #
            # v_ExB_poloidal = hat{e}_pol . v_ExB
            #                = 1 / (J g_yy sqrt(gxx g_zz))
            #                  * (  (  g_yy g_yz g_zz - g_yy g_zz g_yz) d/dx
            #                     + (- g_yz g_yz g_yz + g_yz g_zz g_yy) d/dx
            #                     + (- g_xy g_yz g_zz + g_xy g_zz g_yz) d/dy
            #                     + (  g_yz g_yz g_xz - g_yz g_zz g_xy) d/dy
            #                     + (  g_xy g_yz g_yz - g_xy g_zz g_yy) d/dz
            #                     + (- g_yy g_yz g_xz + g_yy g_zz g_xy) d/dz ) phi
            #                = 1 / (J g_yy sqrt(gxx g_zz))
            #                  * (+ (- g_yz g_yz g_yz + g_yz g_zz g_yy) d/dx
            #                     + (  g_yz g_yz g_xz - g_yz g_zz g_xy) d/dy
            #                     + (  g_xy g_yz g_yz - g_yy g_yz g_xz) d/dz ) phi
            #                = g_yz / (J g_yy sqrt(gxx g_zz))
            #                  * (+ (- g_yz g_yz + g_zz g_yy) d/dx
            #                     + (  g_yz g_xz - g_zz g_xy) d/dy
            #                     + (  g_xy g_yz - g_yy g_xz) d/dz ) phi
            #
            # In the case where the projection of the radial direction onto the
            # poloidal plane is orthogonal to flux surfaces (usual BOUT++ case)
            # this form should be equivalent to projecting onto
            # Grad(y)/sqrt(gyy), which is a unit vector poloidally along a flux
            # surface in that case.
            # |Grad(y)| = sqrt(gyy)
            # Grad(y).e_i = delta^y_i
            # v_ExB . Grad(y)/sqrt(gyy)
            #  = 1/sqrt(gyy) Grad(y).v_ExB
            #  = 1/(g_yy sqrt(gyy)) (0 + g_yz d/dx + 0 + 0 - g_xy d/dz + 0) phi
            #  = 1/(g_yy sqrt(gyy)) (g_yz d/dx - g_xy d/dz) phi
            #
            # In this special case (g_yz g_xz - g_zz g_xy)=0 (see
            # https://bout-dev.readthedocs.io/en/latest/user_docs/coordinates.html#equation-eq-covariantmetric
            # )
            # g_yz / (J g_yy sqrt(gxx g_zz)) (-g_yz g_yz + g_zz g_yy)
            #  = Bt hthe R/Bp Bp/hthe Bp**2/(B**2 hthe**2) 1/sqrt(R**2 Bp**2 R**2)
            #    * (-Bt**2 hthe**2 R**2/Bp**2 + R**2 B**2 hthe**2/Bp**2)
            #  = Bt Bp/(B**2 hthe**2 R)
            #    * R**2 hthe**2 (-Bt**2 + B**2) / Bp**2
            #  = Bt R Bp / (B**2)
            # and
            # g_yz / (g_yy sqrt(gyy))
            #  = Bt hthe R/Bp Bp**2/(B**2 hthe**2) hthe
            #  = Bt R Bp/(B**2)
            # which agree in this limit, and
            # g_yz / (J g_yy sqrt(gxx g_zz)) (g_xy g_yz - g_yy g_xz)
            #  = Bt hthe R/Bp Bp/hthe Bp**2/(B**2 hthe**2) 1/(R Bp) 1/R
            #     * (Bt hthe I R/Bp Bt hthe R/Bp - B**2 hthe**2/Bp**2 I R**2)
            #  = Bt Bp/(B**2 hthe**2 R)
            #     * I R**2 hthe**2/Bp**2 (Bt**2 - B**2)
            #  = - Bt Bp I R / B**2
            # - g_xy / (g_yy sqrt(gyy))
            #  = - Bt hthe I R/Bp Bp**2/(B**2 hthe**2) hthe
            #  = - Bt I R Bp / B**2
            # which also agree
            if zcoord in phi.dims:
                ddz = phi.bout.ddz()
            else:
                ddz = 0.0
            result = (
                ds["g_23"]
                / (ds["J"] * ds["g_22"] * np.sqrt(ds["g11"] * ds["g_33"]))
                * (
                    (-ds["g_23"] ** 2 + ds["g_33"] * ds["g_22"]) * phi.bout.ddx()
                    + (ds["g_23"] * ds["g_13"] - ds["g_33"] * ds["g_12"])
                    * phi.bout.ddy()
                    + (ds["g_12"] * ds["g_23"] - ds["g_22"] * ds["g_13"]) * ddz
                )
            )

            if ds.settings["realistic_geometry"] in ["none", "slab"]:
                # divide by B - see comment above
                result = result / ds["B"]

            result.attrs["standard_name"] = "v_ExB_poloidal"
            result.attrs["long_name"] = "Poloidal component of ExB velocity"
            result.attrs["units"] = ds["U"].units

            return result

        elif direction == "toroidal":
            # |e_z| = sqrt(g_zz)
            # e_z.e_i = g_iz
            # v_ExB_toroidal = 1/sqrt(g_zz) e_z.v_ExB
            #                = 1/(g_yy sqrt(g_zz))
            #                  * (- g_yy g_zz d/dx + g_yz g_yz d/dx + g_xy g_zz d/dy
            #                     - g_yz g_xz d/dy - g_xy g_yz d/dz + g_yy g_xz d/dz)
            #                    phi
            if zcoord in phi.dims:
                ddz = phi.bout.ddz()
            else:
                ddz = 0.0
            result = (-g_22 * ds["g_33"] + ds["g_23"] ** 2) * phi.bout.ddx() + (
                -ds["g_12"] * ds["g_23"] + g_22 * ds["g_13"]
            ) * ddz
            if include_y_derivs:
                result = (
                    result
                    + (ds["g_23"] ** 2) * phi.bout.ddx()
                    + (-ds["g_12"] * ds["g_23"]) * phi.bout.ddz()
                    + (ds["g_12"] * ds["g_33"] - ds["g_23"] * ds["g_13"])
                    * phi.bout.ddy()
                )
            result = result / (g_22 * np.sqrt(ds["g_33"]))

            if ds.settings["realistic_geometry"] in ["none", "slab"]:
                # divide by B - see comment above
                result = result / ds["B"]

            result.attrs["standard_name"] = "v_ExB_toroidal"
            result.attrs["long_name"] = "Toroidal component of ExB velocity"
            result.attrs["units"] = ds["U"].units

            return result

        else:
            raise ValueError(f'Unrecognised value "{direction}" for direction')

    def Curl_b_B(self, direction=None):
        """
        Calculates Curl(b/B)

        The radial component is `Curl(b/B).Grad(psi)/|Grad(psi)|`.
        The poloidal component is `Curl(b/B).Grad(x) x e_z / |Grad(x)| |e_z|`.
        The toroidal component is `Curl(b/B).e_z/|e_z|`.

        Parameters
        ----------
        direction : str, optional
            By default, calculates radial, poloidal and toroidal components and returns
            them as a dict. Can pass "radial", "poloidal" or "toroidal" to return a
            specific component.
        """
        if direction is None:
            return {
                d: self.Curl_b_B(d, include_y_derivs=include_y_derivs)
                for d in ["radial", "poloidal", "toroidal"]
            }

        ds = self.data

        B = ds["B"]

        if ds.settings["realistic_geometry"] in ["none", "slab"]:
            # In slab, only have a binormal (labelled "toroidal" here) component of Curl(b/B)
            if direction in ["radial", "poloidal"]:
                return 0.0
            elif direction == "toroidal":
                # g_0 := 2 rho_s0 / R_c
                # Curl(b/B) = 2/(B R_c) e_z = g_0/(B rho_s0) e_z
                # and as coordinates are orthonormal
                # Curl_b_B_toroidal = e_z.Curl(b/B) = g_0 / (B rho_s0)
                return ds.params["dimless"]["g0"] / (ds["B"] * ds.params["rho_s0"])
            raise ValueError(f'Unrecognised value "{direction}" for direction')
        else:
            # bxcvx, bxcvy and bxcvz are the contravariant components of B/2*Curl(b/B)
            # Curl(b/B) = 2/B * (bxcvx e_x + bxcvy e_y + bxcvz e_z)
            if direction == "radial":
                # x := psi
                # |Grad(x)| = sqrt(g11)
                # Grad(x).e_i = delta^x_i
                # Curl_b_B_radial = 1/sqrt(g11) Grad(x).Curl(b/B)
                #                 = 2/(B sqrt(g11)) bxcvx
                return 2.0 / (B * np.sqrt(ds["g11"])) * ds["bxcvx"]
            elif direction == "poloidal":
                # Need a a unit vector hat{e}_pol in the poloidal direction. e_z is
                # in the toroidal direction and Grad(x) is orthogonal to flux
                # surfaces, so their cross product is in the poloidal direction
                # (within flux surfaces). e_z and Grad(x) are also always
                # orthogonal, so the magnitude of their cross product is the
                # product of their magnitudes. Therefore
                #   hat{e}_pol = (e_z x Grad(x))/|Grad(x)||e_z|
                #   hat{e}_pol = (e_z x Grad(x))/sqrt(gxx g_zz)
                #
                #   hat{e}_pol . e_i
                #    = (e_z x Grad(x)).e_i/sqrt(g11 g_zz)
                #    = (g_xz Grad(x) + g_yz Grad(y) + g_zz Grad(z))
                #       x Grad(x) . e_i / sqrt(gxx g_zz)
                # Using Grad(x^i) x Grad(x^j) = e_k/J  [i,j,k cyclic combination of x,y,z]
                #   hat{e}_pol . e_i
                #    = (- g_yz e_z + g_zz e_y) . e_i / (J sqrt(gxx g_zz))
                #    = (- g_yz g_iz + g_zz e_iy) / (J sqrt(gxx g_zz))
                #    = (- g_yz g_iz + g_zz g_iy) / (J sqrt(gxx g_zz))
                # Curl_b_B_poloidal = hat{e}_pol.Curl(b/B)
                #                   = 2 / (B J sqrt(gxx g_zz))
                #                     * (   (-g_yz g_xz + g_zz g_xy) bxcvx
                #                         + (-g_yz g_yz + g_zz g_yy) bxcvy
                #                         + (-g_yz g_zz + g_zz g_yz) bxcvz)
                return (
                    2.0
                    / (B * ds["J"] * np.sqrt(ds["g11"] * ds["g_33"]))
                    * (
                        (-ds["g_23"] * ds["g_13"] + ds["g_33"] * ds["g_12"])
                        * ds["bxcvx"]
                        + (-ds["g_23"] ** 2 + ds["g_33"] * ds["g_22"]) * ds["bxcvy"]
                        + (-ds["g_23"] * ds["g_33"] + ds["g_33"] * ds["g_23"])
                        * ds["bxcvz"]
                    )
                )
            elif direction == "toroidal":
                # |e_z| = sqrt(g_33)
                # e_z.e_i = g_iz
                # Curl_b_B_toroidal = 1/sqrt(g_33) e_z.Curl(b/B)
                #              = 2/(B sqrt(g_33) (g_13 bxcvx + g_23 bxcvy + g_33 bxcvyz)
                return (
                    2.0
                    / (B * np.sqrt(ds["g_33"]))
                    * (
                        ds["g_13"] * ds["bxcvx"]
                        + ds["g_23"] * ds["bxcvy"]
                        + ds["g_33"] * ds["bxcvz"]
                    )
                )
            else:
                raise ValueError(f'Unrecognised value "{direction}" for direction')

    def v_mag_e(self, direction=None, *, include_y_derivs=None):
        """
        Calculates electron magnetic drift velocity v_mag_e=T_e/e*Curl(b/B)

        Want to calculate fluxes due to 'curvature' terms in STORM equations. For
        example term in density equation is
        +1/e*Curl(b/B).Grad(p_e)
        which comes from the divergence of the diamagnetic particle flux Div(n*v_dia,e).
        However, we can follow HERMES [Dudson&Leddy (2017)] and use a magnetic drift
        velocity instead, because the diamagnetic particle flux is

        .. code::

            n*v_dia_e = Grad(p_e) x b / e B
                      = Curl(p_e b / e B) - p_e / e * Curl(b/B)

        The first term is large but divergence-free so does not contribute to particle
        transport. The second, n*v_mag_e, is small but its divergence gives the
        'curvature' term in the density equation. When trying to integrate fluxes
        (especially as STORM's numerical method is not conservative), it is probably
        better to avoid integrating a large quantity that should exactly cancel, so
        better to calculate a flux from v_mag_e than v_dia_e.

        The radial component is `v_mag_e.Grad(psi)/|Grad(psi)|`.
        The poloidal component is `v_mag_e.Grad(x) x e_z / |Grad(x)| |e_z|`.
        The toroidal component is `v_mag_e.e_z/|e_z|`.

        Parameters
        ----------
        direction : str, optional
            By default, calculates radial, poloidal and toroidal components and returns
            them as a dict. Can pass "radial", "poloidal" or "toroidal" to return a
            specific component.
        include_y_derivs : bool, optional
            Parallel derivatives are often neglected as small. In this case poloidal
            advection may be neglected - here include_y_derivs=False sets the poloidal
            component to zero. Default is True to match STORM.
        """
        if direction is None:
            return {
                d: self.v_mag_e(d, include_y_derivs=include_y_derivs)
                for d in ["radial", "poloidal", "toroidal"]
            }

        if include_y_derivs is None:
            # Default to STORM default behaviour
            include_y_derivs = True

        ds = self.data

        if direction == "poloidal" and not include_y_derivs:
            return 0.0

        # v_mag_e = -T_e[J]/e Curl(b/B)
        #         = -T_e[eV]*e/e Curl(b/B) = -T_e[eV] Curl(b/B)
        return -ds["T"] * ds.storm.Curl_b_B(direction)

    def v_par_e(self, direction):
        """
        Calculates the components of the parallel velocity on orthonormal basis vectors

        The radial component is `v_par_e.Grad(psi)/|Grad(psi)|`.
        The poloidal component is `v_par_e.Grad(x) x e_z / |Grad(x)| |e_z|`.
        The toroidal component is `v_par_e.e_z/|e_z|`.

        Parameters
        ----------
        direction : str, optional
            By default, calculates radial, poloidal and toroidal components and returns
            them as a dict. Can pass "radial", "poloidal" or "toroidal" to return a
            specific component.
        """
        if direction is None:
            return {d: self.v_par_e(d) for d in ["radial", "poloidal", "toroidal"]}

        ds = self.data

        if direction == "radial":
            # b.Grad(psi) = 0
            return 0.0
        elif direction == "poloidal":
            # Need a a unit vector hat{e}_pol in the poloidal direction. e_z is
            # in the toroidal direction and Grad(x) is orthogonal to flux
            # surfaces, so their cross product is in the poloidal direction
            # (within flux surfaces). e_z and Grad(x) are also always
            # orthogonal, so the magnitude of their cross product is the
            # product of their magnitudes. Therefore
            #   hat{e}_pol = (e_z x Grad(x))/|Grad(x)||e_z|
            #   hat{e}_pol = (e_z x Grad(x))/sqrt(gxx g_zz)
            #
            #   hat{e}_pol . e_i
            #    = (e_z x Grad(x)).e_i/sqrt(g11 g_zz)
            #    = (g_xz Grad(x) + g_yz Grad(y) + g_zz Grad(z))
            #       x Grad(x) . e_i / sqrt(gxx g_zz)
            # Using Grad(x^i) x Grad(x^j) = e_k/J  [i,j,k cyclic combination of x,y,z]
            #   hat{e}_pol . e_i
            #    = (- g_yz e_z + g_zz e_y) . e_i / (J sqrt(gxx g_zz))
            #    = (- g_yz g_iz + g_zz e_iy) / (J sqrt(gxx g_zz))
            #    = (- g_yz g_iz + g_zz g_iy) / (J sqrt(gxx g_zz))
            # hat{e}_pol . v_par_e = hat{e}_pol . b |v_par_e|
            #                      = hat{e}_pol . e_y/sqrt(g_yy) |v_par_e|
            #                      = 1/(J sqrt(gxx g_yy g_zz))
            #                        * (-g_yz g_yz + g_zz g_yy) |v_par_e|
            return (
                ds["V"]
                * (-ds["g_23"] ** 2 + ds["g_33"] * ds["g_22"])
                / (ds["J"] * np.sqrt(ds["g11"] * ds["g_22"] * ds["g_33"]))
            )
        elif direction == "toroidal":
            # v_par_e.e_z/|e_z| = |v_par_e|b.e_z/sqrt(g_zz)
            #                   = |v_par_e| e_y/sqrt(g_yy) . e_z/sqrt(g_zz)
            #                   = |v_par_e| g_yz / sqrt(g_yy g_zz)
            return ds["V"] * ds["g_23"] / np.sqrt(ds["g_22"] * ds["g_33"])
        else:
            raise ValueError(f'Unrecognised value "{direction}" for direction')

    def particle_flux_diss_e(self, direction, *, include_y_derivs=None):
        """
        Calculate the components of the particle flux that gives the density diffusion
        term, -mu_n Grad_perp(n).

        mu_n is calculated either including density and temperature variation if
        uniform_diss_paras=True or just as mu_n0 if uniform_diss_paras=False.

        The radial component is `particle_flux_diss_e.Grad(psi)/|Grad(psi)|`.
        The poloidal component is `particle_flux_diss_e.Grad(x) x e_z / |Grad(x)| |e_z|`.
        The toroidal component is `particle_flux_diss_e.e_z/|e_z|`.

        Parameters
        ----------
        direction : str, optional
            By default, calculates radial, poloidal and toroidal components and returns
            them as a dict. Can pass "radial", "poloidal" or "toroidal" to return a
            specific component.
        include_y_derivs : bool, optional
            Include parallel derivatives in calculation of particle_flux_diss_e. These
            are often neglected as small. Default is False to match STORM.
        """
        if direction is None:
            return {
                d: self.particle_flux_diss_e(d)
                for d in ["radial", "poloidal", "toroidal"]
            }

        if include_y_derivs is None:
            # Default to STORM default behaviour
            include_y_derivs = False

        ds = self.data

        mu_n0 = ds.params["dimless"]["mu_n0"]
        if ds.metadata["storm_normalised"] == False:
            mu_n0 *= ds.params["rho_s0"] ** 2 * ds.params["omega_i0"]

        n = ds["n"]

        if ds.settings["uniform_diss_paras"]:
            mu_n = mu_n0
        else:
            mu_n = mu_n0 * n / ds.params["n_0"] * np.sqrt(ds.params["T_e0"] / ds["T"])

        # Grad_perp(n) = Grad(n) - b b.Grad(n)
        if direction == "radial":
            # particle_flux_diss_e.Grad(psi)/|Grad(psi)|
            #  = -mu_n Grad_perp(n) . Grad(psi) / |Grad(psi)|
            #  = -mu_n (Grad(n) - b b.Grad(n)) . Grad(psi) / |Grad(psi)|
            #    note b.Grad(psi) is always 0 because e_i.Grad(x^j) = delta_i^j
            #  = -mu_n Grad(n) . Grad(psi) / |Grad(psi)|
            #  = -mu_n Grad_i(n) g^ij Grad_j(psi) / |Grad(psi)|
            #  = -mu_n (gxx dn/dx + gxy dn/dy + gxz dn/dz) |Grad(psi)| / |Grad(psi)|
            #  = -mu_n (gxx dn/dx + gxy dn/dy + gxz dn/dz)
            result = -(ds["g11"] * n.bout.ddx() + ds["g13"] * n.bout.ddz())
            if include_y_derivs:
                result += -ds["g12"] * n.bout.ddy()
            # result *= mu_n
            result = result / np.sqrt(ds["g11"]) * mu_n
            return result

        elif direction == "poloidal":
            if not include_y_derivs:
                # Poloidal component neglected by STORM
                return 0.0
            else:
                # Need a a unit vector hat{e}_pol in the poloidal direction. e_z is
                # in the toroidal direction and Grad(x) is orthogonal to flux
                # surfaces, so their cross product is in the poloidal direction
                # (within flux surfaces). e_z and Grad(x) are also always
                # orthogonal, so the magnitude of their cross product is the
                # product of their magnitudes. Therefore
                #   hat{e}_pol = (e_z x Grad(x))/|Grad(x)||e_z|
                #   hat{e}_pol = (e_z x Grad(x))/sqrt(gxx g_zz)
                #
                #   hat{e}_pol . e_i
                #    = (e_z x Grad(x)).e_i/sqrt(g11 g_zz)
                #    = (g_xz Grad(x) + g_yz Grad(y) + g_zz Grad(z))
                #       x Grad(x) . e_i / sqrt(gxx g_zz)
                # Using Grad(x^i) x Grad(x^j) = e_k/J  [i,j,k cyclic combination of x,y,z]
                #   hat{e}_pol . e_i
                #    = (- g_yz e_z + g_zz e_y) . e_i / (J sqrt(gxx g_zz))
                #    = (- g_yz g_iz + g_zz g_iy) / (J sqrt(gxx g_zz))
                #
                #   hat{e}_pol . Grad(x^i)
                #    = (- g_yz e_z + g_zz e_y) . Grad(x^i) / (J sqrt(gxx g_zz))
                #    = (- g_yz delta^i_z + g_zz delta^i_y) / (J sqrt(gxx g_zz))
                # hat{e}_pol . particle_flux_diss_e
                #  = -hat{e}_pol . mu_n Grad_perp(n)
                #  = -hat{e}_pol . mu_n (Grad(n) - b b.Grad(n))
                #  = -hat{e}_pol . mu_n (dn/dx Grad(x) + dn/dy Grad(y) + dn/dz Grad(z)
                #                      - b b.Grad(n))
                #  = -hat{e}_pol . mu_n (dn/dx Grad(x) + dn/dy Grad(y) + dn/dz Grad(z)
                #                      - e_y e_y.Grad(n)/g_yy)
                #  = -hat{e}_pol . mu_n (dn/dx Grad(x) + dn/dy Grad(y) + dn/dz Grad(z)
                #                      - e_y/g_yy dn/dy)
                #  = -mu_n (0 + g_zz dn/dy - g_yz dn/dz
                #          - (-g_yz g_yz + g_zz g_yy)/g_yy dn/dy) / (J sqrt(gxx g_zz))
                #  = -mu_n (g_zz dn/dy - g_yz dn/dz
                #          + g_yz g_yz/g_yy dn/dy - g_zz dn/dy) / (J sqrt(gxx g_zz))
                #  = -mu_n (- g_yz dn/dz + g_yz g_yz/g_yy dn/dy) / (J sqrt(gxx g_zz))
                #  = -mu_n g_yz (- dn/dz + g_yz/g_yy dn/dy) / (J sqrt(gxx g_zz))
                result = -n.bout.ddz()
                if include_y_derivs:
                    result = result + ds["g_23"] / ds["g_22"] * n.bout.ddy()
                result = (
                    -mu_n
                    * ds["g_23"]
                    * result
                    / (ds["J"] * np.sqrt(ds["g11"] * ds["g_33"]))
                )
                return result
        elif direction == "toroidal":
            # particle_flux_diss_e.e_z/|e_z| = mu_n Grad_perp(n) . e_z / |e_z|
            #  = -mu_n (Grad(n) - b b.Grad(n)) . e_z / |e_z|
            #  = -mu_n (Grad(n) - b b.Grad(n)) . e_z / sqrt(g_zz)
            #  = -mu_n (Grad(n) - e_y e_y.Grad(n)/g_yy) . e_z/sqrt(g_zz)
            #  = -mu_n (dn/dz - g_yz dn/dy/g_yy) / sqrt(g_zz)
            result = -mu_n * n.bout.ddz()
            if include_y_derivs:
                result = result + mu_n * ds["g_23"] * n.bout.ddy() / ds["g_22"]
            result = result / np.sqrt(ds["g_33"])
            return result
        else:
            raise ValueError(f'Unrecognised value "{direction}" for direction')

    def electron_particle_flux(self, direction, *, include_y_derivs=None):
        """
        Calculates particle flux of electrons

        `electron particle flux = n * (v_ExB + v_mag_e + v_par_e + v_diss_e)`

        The radial component is `particle_flux_e.Grad(psi)/|Grad(psi)|`.
        The poloidal component is `particle_flux_e.Grad(x) x e_z / |Grad(x)| |e_z|`.
        The toroidal component is `particle_flux_e.e_z/|e_z|`.

        Parameters
        ----------
        direction : str, optional
            By default, calculates radial, poloidal and toroidal components and returns
            them as a dict. Can pass "radial", "poloidal" or "toroidal" to return a
            specific component.
        include_y_derivs : None or bool, optional
            Include parallel derivatives in calculation of v_diss_e. These are often
            neglected as small. Default (None) matches STORM defaults for v_ExB,
            v_mag_e and v_diss_e.
        """
        if direction is None:
            return {
                d: self.electron_particle_flux(d, include_y_derivs=include_y_derivs)
                for d in ["radial", "poloidal", "toroidal"]
            }

        result = self.data["n"] * (
            self.v_ExB(direction, include_y_derivs=include_y_derivs)
            + self.v_mag_e(direction, include_y_derivs=include_y_derivs)
            + self.v_par_e(direction)
        ) + self.particle_flux_diss_e(direction, include_y_derivs=include_y_derivs)

        result.name = f"{direction} electron particle flux"
        result.attrs["standard_name"] = result.name
        result.attrs["long_name"] = result.name
        if result.metadata["storm_normalised"]:
            result.attrs["units"] = "Omega_i0 rho_s0^-2"
        else:
            result.attrs["units"] = "s^-1 m^-2"

        return result

    def electron_mfp(self, n=None, T_e=None):
        """
        Calculates the electron collisional mean-free-path.

        As STORM assumes a pure, Z=1 plasma the mean-free-paths for electron-electron or
        electron-ion collisions are the same, so this method is valid for either.

        Parameters
        ----------

        n : float or array, optional
            Density to use. If not supplied will use the density array in this
            StormDataset.
        T_e : float or array, optional
            Electron temperature to use (in eV). If not supplied will use the temperature
            array in this StormDataset.
        """

        ds = self.data

        if ds.metadata["storm_normalised"]:
            raise ValueError(
                "electron_mfp() is only implemented for unnormalised Datasets so far"
            )

        if T_e is None:
            T_e = ds["T"]
        if n is None:
            n = ds["n"]
        m_e = ds.params["m_e"]
        epsilon0 = ds.params["epsilon_0"]
        echarge = ds.params["e"]
        v_te = np.sqrt(T_e * echarge / m_e)
        loglambda = ds.params["loglambda"]
        nu_e = (
            n
            * echarge**4
            * loglambda
            / (
                np.sqrt(m_e)
                * epsilon0**2
                * 3.0
                * (2.0 * np.pi * echarge * T_e) ** 1.5
            )
        )
        mfp_e = v_te / nu_e
        if isinstance(mfp_e, xr.DataArray):
            mfp_e.attrs["units"] = "m"
            mfp_e.attrs["standard_name"] = "electron mfp"
            mfp_e.attrs["long_name"] = "electron mean-free-path"

        return mfp_e

    def ion_mfp(self, n=None, T_i=None):
        """
        Calculates the ion collisional mean-free-path.

        Parameters
        ----------

        n : float or array, optional
            Density to use. If not supplied will use the density array in this
            StormDataset.
        T_i : float or array, optional
            Electron temperature to use (in eV). If not supplied will use the temperature
            array in this StormDataset.
        """

        ds = self.data

        if ds.metadata["storm_normalised"]:
            raise ValueError(
                "ion_mfp() is only implemented for unnormalised Datasets so far"
            )

        if T_i is None:
            T_i = ds["T_i"]
        if n is None:
            n = ds["n"]
        m_i = ds.params["m_i"]
        epsilon0 = ds.params["epsilon_0"]
        echarge = ds.params["e"]
        v_ti = np.sqrt(T_i * echarge / m_i)
        loglambda = ds.params["loglambda"]
        nu_i = (
            n
            * echarge**4
            * loglambda
            / (
                np.sqrt(m_i)
                * epsilon0**2
                * 3.0
                * (2.0 * np.pi * echarge * T_i) ** 1.5
            )
        )
        mfp_i = v_ti / nu_i
        if isinstance(mfp_i, xr.DataArray):
            mfp_i.attrs["units"] = "m"
            mfp_i.attrs["standard_name"] = "ion mfp"
            mfp_i.attrs["long_name"] = "ion mean-free-path"

        return mfp_i

    def quiver_plot(
        self,
        v=None,
        *,
        radial=None,
        poloidal=None,
        toroidal=None,
        isel=None,
        step=None,
        background=None,
        ax=None,
        poloidal_plot=True,
        quiver_kwargs=None,
        **kwargs,
    ):
        """
        Make a quiver plot of a vector (for example the ExB flow velocity). Should be
        used on a Dataset with 2 spatial dimensions and no time dimension.

        When plots are made against grid coordinates (e.g "radial" and "poloidal"), the
        vector components are still drawn according to the physical sizes in the
        corresponding directions, so may not follow the streamlines in the plot
        coordinates.

        Parameters
        ----------
        v : dict of DataArray, optional
            Dict of DataArray as returned, for example, by `StormDataset.v_ExB()`. If
            not given, the velocity components needed for the plot must be passed in
            `radial`, `poloidal`, and `toroidal` arguments.
        radial : DataArray, optional
            If `v` is not passed, and the radial velocity component is required, it
            should be passed to this argument. Ignored if `v` is not `None`.
        poloidal : DataArray, optional
            If `v` is not passed, and the poloidal velocity component is required, it
            should be passed to this argument. Ignored if `v` is not `None`.
        toroidal : DataArray, optional
            If `v` is not passed, and the toroidal velocity component is required, it
            should be passed to this argument. Ignored if `v` is not `None`.
        isel : dict, optional
            `dict` to be passed to `isel()` to select a slice of the Dataset and
            velocity components for plotting. Not needed if both Dataset and velocity
            components already have the correct dimensions (2 spatial, no time).
        step : int or (int, int), optional
            A `step`, so the vector is drawn only every `step` points. Can help to make
            arrows more visible. If `step` is an integer it is used for both dimensions,
            if it is `(int, int)` the values are applied to the first and second
            dimensions.
        background : str, optional
            Variable to use for color-map plot in the background.
        ax : Matplotlib.Axes.Axes
            Pass an existing Axes object to add this plot to it. By default creates a
            new `Figure`.
        poloidal_plot : bool, default True
            If plotting x- and y-dimensions, make a 'poloidal plot' with R-Z axes.
        quiver_kwargs : dict
            Optional keyword arguments passed to xarray.Dataset.plot.quiver()
        **kwargs : dict
            Keyword arguments passed to the background plotting function:
            `BoutDataArray.pcolormesh()` (for poloidal plots) or
            `xarray.DataArray.plot.pcolormesh()` (otherwise).
        """

        # Make a copy so we don't add new variables to original ds
        ds = self.data.copy()

        if v is None:
            v = {}

        if step is None:
            xstep, ystep = None, None
        elif isinstance(step, int):
            xstep, ystep = step, step
        else:
            xstep, ystep = step

        if ax is None:
            _, ax = plt.subplots()

        if quiver_kwargs is None:
            quiver_kwargs = {}

        ds = ds.isel(isel)

        tcoord = ds.metadata["bout_tdim"]
        xcoord = ds.metadata["bout_xdim"]
        ycoord = ds.metadata["bout_ydim"]
        zcoord = ds.metadata["bout_zdim"]

        dims = list(ds["phi"].dims)

        if len(dims) > 2:
            raise ValueError(f"Too many dimensions for quiver_plot(). Got {list(dims)}")
        elif len(dims) < 2:
            raise ValueError(
                f"Not enough dimensions for quiver_plot(). Got {list(dims)}"
            )
        elif tcoord in dims:
            raise ValueError("quiver_plot() does not work with a time dimension")
        elif xcoord in dims and ycoord in dims and poloidal_plot:
            # R-Z plot

            # velocity components relative to flux surface
            vx = v.get("radial", radial).isel(isel)
            vy = v.get("poloidal", poloidal).isel(isel)

            # Poloidal magnetic field components
            br = ds["Brxy"] / ds["Bpxy"]
            bz = ds["Bzxy"] / ds["Bpxy"]

            # Unit vector in 'radial' direction is (-bz, br)
            # Unit vector in 'poloidal' direction is (br, bz)

            vr = -vx * bz + vy * br
            vz = vx * br + vy * bz

            ds["quiver_plot_u"] = vr.fillna(0.0)
            ds["quiver_plot_v"] = vz.fillna(0.0)

            if background is not None:
                ds[background].bout.pcolormesh(ax=ax, **kwargs)

            return ds.isel(
                {xcoord: slice(None, None, xstep), ycoord: slice(None, None, ystep)}
            ).plot.quiver(
                x="R",
                y="Z",
                ax=ax,
                u="quiver_plot_u",
                v="quiver_plot_v",
                angles="xy",
                **quiver_kwargs,
            )
        elif xcoord in dims and ycoord in dims:
            # x-y plot

            # velocity components relative to flux surface
            vx = v.get("radial", radial).isel(isel)
            vy = v.get("poloidal", poloidal).isel(isel)

            ds["quiver_plot_u"] = vx.fillna(0.0)
            ds["quiver_plot_v"] = vy.fillna(0.0)

            if background is not None:
                ds[background].plot.pcolormesh(x=xcoord, y=ycoord, ax=ax, **kwargs)

            return ds.isel(
                {xcoord: slice(None, None, xstep), ycoord: slice(None, None, ystep)}
            ).plot.quiver(
                x=xcoord,
                y=ycoord,
                ax=ax,
                u="quiver_plot_u",
                v="quiver_plot_v",
                angles="xy",
                **quiver_kwargs,
            )
        elif xcoord in dims and zcoord in dims:
            # x-z plot

            # velocity components
            vx = v.get("radial", radial).isel(isel)
            vz = v.get("toroidal", toroidal).isel(isel)

            ds["quiver_plot_u"] = vx.fillna(0.0)
            ds["quiver_plot_v"] = vz.fillna(0.0)

            if background is not None:
                ds[background].plot.pcolormesh(x=xcoord, y=zcoord, ax=ax, **kwargs)

            return ds.isel(
                {xcoord: slice(None, None, xstep), zcoord: slice(None, None, ystep)}
            ).plot.quiver(
                x=xcoord,
                y=zcoord,
                ax=ax,
                u="quiver_plot_u",
                v="quiver_plot_v",
                angles="xy",
                **quiver_kwargs,
            )
        elif ycoord in dims and zcoord in dims:
            # y-z plot

            # velocity components
            vy = v.get("poloidal", poloidal).isel(isel)
            vz = v.get("toroidal", toroidal).isel(isel)

            ds["quiver_plot_u"] = vy.fillna(0.0)
            ds["quiver_plot_v"] = vz.fillna(0.0)

            if background is not None:
                ds[background].plot.pcolormesh(x=ycoord, y=zcoord, ax=ax, **kwargs)

            return ds.isel(
                {ycoord: slice(None, None, xstep), zcoord: slice(None, None, ystep)}
            ).plot.quiver(
                x=ycoord,
                y=zcoord,
                ax=ax,
                u="quiver_plot_u",
                v="quiver_plot_v",
                angles="xy",
                **quiver_kwargs,
            )
        else:
            raise ValueError(f"Dimensions {d} not supported")


@register_dataarray_accessor("storm")
class StormDataArrayAccessor(BoutDataArrayAccessor):
    """
    Class specifically for accessing analysis methods which are to be
    performed on data obtained from a simulation using the STORM module
    for BOUT++.

    Implements methods for analysing single variables.
    """

    def __init__(self, da):
        super().__init__(da)

    def norms(self):
        return calc_norms(self.data)

    def std_dev(self, dim=None):
        """
        Compute standard deviation along the given dimensions.

        Parameters
        ----------
        dim : str or list of str, optional

        Returns
        -------
        xarray.DataArray
        """

        if isinstance(dim, str):
            dim = [dim]

        result = apply_ufunc(
            np.std,
            self.data,
            input_core_dims=[dim],
            kwargs={"axis": tuple(np.arange(-len(dim), 0))},
            dask="parallelized",
            output_dtypes=[np.float64],
            keep_attrs=True,
        )

        result = _update_name(result, suffix="_std_dev")
        return result

    def std_score(self, dim=None):
        """
        Compute standard score (normalised standard deviation) along the given
        dimensions.

        i.e. `z = (x-mu)/sigma`, where z is standard score (or z-score),
        mu is the mean, and sigma is the standard deviation.

        Parameters
        ----------
        dim : str or list of str, optional

        Returns
        -------
        xarray.DataArray
        """

        data = self.data
        sigma = data.storm.std_dev(dim=dim)
        mean = data.mean(dim=dim)
        result = (data - mean) / sigma

        result = _update_name(result, suffix="_std_score")
        return _update_units(result, replace="std score")

    def skewness(self, dim):
        """
        Compute skewness along the given dimension.

        Parameters
        ----------
        dim : str

        Returns
        -------
        xarray.DataArray
        """
        result = apply_ufunc(
            daskstats.skew,
            self.data,
            input_core_dims=[[dim]],
            kwargs={"axis": -1},
            dask="allowed",
            output_dtypes=[np.float64],
            keep_attrs=True,
        )

        result = _update_name(result, suffix="_skewness")
        return _update_units(result, replace="")

    def kurtosis(self, dim):
        """
        Compute kurtosis along the given dimension.

        Parameters
        ----------
        dim : str

        Returns
        -------
        xarray.DataArray
        """

        result = apply_ufunc(
            daskstats.kurtosis,
            self.data,
            input_core_dims=[[dim]],
            kwargs={"axis": -1},
            dask="allowed",
            output_dtypes=[np.float64],
            keep_attrs=True,
        )

        result = _update_name(result, suffix="_kurtosis")
        return _update_units(result, replace="")

    def pdf(self, dim=None, bins=None):
        """
        Calculate the normalised probability distribution function over one or more dimensions.

        Wraps xhistogram.xarray.histogram.
        Requires PR #17 to xhistogram to compute density correctly.

        Parameters
        ----------
        dim : str or list of str, optional
            Dimensions over which which the histogram is computed. The default is to
            compute the histogram of the flattened array.
        bins : int or array_like, optional
            If int, the number of bins for all arguments in args.
            If array_like, the bin edges for all arguments in args.

        Returns
        -------
        xarray.DataArray
        """

        from xhistogram.xarray import histogram

        da = self.data
        mean = da.mean(dim)
        std = da.std(dim)
        fluc = (da - mean) / std

        hist = histogram(fluc, dim=[dim], bins=[bins], density=True)

        var = da.name
        bin_coord = f"{var}_bin"
        hist[bin_coord].attrs["long_name"] = f"({var} - <{var}>_{dim}) / sigma_{var}"
        hist[bin_coord].attrs["units"] = ""
        hist.name = "PDF"
        return hist

    def power_spectra(self, dim=None, **kwargs):
        """
        Compute power spectrum along the given dimensions.

        Parameters
        ----------
        dim : str or list of str, optional
        kwargs : Passed on to xrft.power_spectrum

        Returns
        -------
        xarray.DataArray
        """

        from xrft import power_spectrum

        ps = power_spectrum(self.data, dim=dim, **kwargs)

        if isinstance(dim, str):
            dim = [dim]
        freq_dim = ["freq_" + d for d in dim]

        abs_ps = np.abs(ps)
        abs_ps.name = f"|{self.data.name} power spectrum|"

        # We know the normalisation for the perpendicular directions
        # TODO the same but for 'time'
        for d, f_d in zip(dim, freq_dim):
            if d in ["radial", "binormal"]:
                if self.data.metadata["storm_normalised"]:
                    raise ValueError(
                        "power_spectra() is only implemented for unnormalised Datasets "
                        "so far"
                    )

                k_perp_rho_s = abs_ps.coords[f_d] * self.data.attrs["params"]["rho_s0"]
                new_perp_freq_coord = f"k_{d}_rho_s"
                abs_ps = abs_ps.assign_coords(**{new_perp_freq_coord: k_perp_rho_s})

        return abs_ps

    def conditional_average(self, coord=None, level=2.5, delta=None):
        """
        Conditional average along a single dimension.

        Finds positions where signal exceeds more than value level,
        and are separated by more than delta along coord.

        Assumes input DataArray has already had the mean subtracted
        and been divided by its standard deviation.

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
        crossings : xarray.DataArray
            Indices of points where level was exceeded, which are separated
            by more than delta.
        """
        return conditional_average(self.data, coord=coord, level=level, delta=delta)

    def waiting_times(self, coord=None, level=2.5, delta=None):
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
        return waiting_times(self.data, coord=coord, level=level, delta=delta)
