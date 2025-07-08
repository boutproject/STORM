from warnings import warn

from boutdata.data import BoutOptions
import numpy as np

from xbout import open_boutdataset
from xbout.load import _open_grid
from xbout.geometries import apply_geometry
from xbout.region import _create_regions_toroidal
from xbout.utils import _set_attrs_on_all_vars

from .coords import add_time, add_slab_coords


# TODO register other possible geometries (storm_salpha, storm_doublenull)


# Particular options to load into settings
STORM_SETTINGS = {"boussinesq": 1}
# storm3d-specific settings
STORM3D_SETTINGS = {
    "normalise_lengths": "false",
    "normalise_all": "false",
    "uniform_diss_paras": "false",
}


def open_stormdataset(
    datapath="./BOUT.dmp.*.nc",
    inputfilepath=None,
    gridfilepath=None,
    chunks={},
    run_name=None,
    info=True,
    unnormalise=True,
    **kwargs,
):
    """
    Load a dataset from a set of STORM output files, including the input
    options file.

    Parameters
    ----------
    datapath : str, optional
    chunks : dict, optional
    inputfilepath : str, optional
    gridfilepath : str, optional
    run_name : str, optional
    info : bool, optional
    unnormalise : bool, default True
        Convert variables to SI units - set to False to leave them normalised.
    kwargs : optional
        Passed on to xbout.open_boutdataset

    Returns
    -------
    ds : xarray.Dataset
    """

    # Geometry is None because we don't know it yet
    ds = open_boutdataset(
        datapath=datapath,
        inputfilepath=inputfilepath,
        geometry=None,
        chunks=chunks,
        run_name=run_name,
        info=False,
        **kwargs,
    )

    if not ds.options:
        raise ValueError(
            "Can't add units or calculate Storm normalisations without the input file "
            "options"
        )

    if "metadata" not in ds.attrs:
        raise ValueError("Library error: metadata missing from BoutDataset")

    # This will be set to False if apply_unnormalise() is called.
    # Need to use int instead of bool so xarray can save to netcdf.
    # Keep current value if one exists, because we might be reloading a Dataset that
    # was saved by xSTORM already
    ds.metadata["storm_normalised"] = ds.metadata.get("storm_normalised", 1)

    # Read various STORM options from the input file or output files
    settings = read_storm_settings(ds)
    ds.attrs["settings"] = settings

    # Now apply relevant type of BOUT++ geometry
    # TODO option for s-alpha
    if gridfilepath is not None:
        grid = _open_grid(
            gridfilepath,
            chunks=chunks,
            keep_xboundaries=ds.metadata["keep_xboundaries"],
            keep_yboundaries=ds.metadata["keep_yboundaries"],
            mxg=ds.metadata["MXG"],
        )
    else:
        grid = None
    if ds.settings["realistic_geometry"] != "slab":
        ds = apply_geometry(ds, geometry_name="toroidal", grid=grid)

    # Calculate and store STORM-specific plasma parameters
    params = calc_plasma_params(ds.options)
    _set_attrs_on_all_vars(ds, "params", params)

    if "time" in ds.dims:
        # Reloading data that was already processed by xSTORM: do not need to add
        # coordinates or unnormalise
        if ds.settings["realistic_geometry"] == "slab":
            # Create regions even with slab geometry as these are used by some xBOUT/xSTORM
            # methods
            ds = _create_regions_toroidal(ds)
        if info:
            print("Read in:\n")
            print(ds)
        return ds

    if ds.settings["realistic_geometry"] == "slab":
        # Add slab coordinates
        ds = add_slab_coords(ds)
    else:
        # theta and zeta are angle or angle-like coordinates so do not need any units or
        # conversion.
        # psi_poloidal is read in physical units from the grid file, so also does not need
        # any conversion here.
        pass

    # Un-normalise data variables so they are in physical units, following
    # STORM conventions
    # Have to special-case B because of how storm2d treats it
    if "B" not in ds:
        if "Bxy" in ds:
            # storm3d
            ds = ds.rename(Bxy="B")
        else:
            raise ValueError("No magnetic field data found")
    if unnormalise:
        ds = apply_unnormalise(ds)
    else:
        ds = add_normalised_units(ds)

    if ds.settings["realistic_geometry"] == "slab":
        # Create regions even with slab geometry as these are used by some xBOUT/xSTORM
        # methods
        ds = _create_regions_toroidal(ds)

    # Drop y dimension if necessary
    ds = ds.squeeze(drop=True)

    if info:
        print("Read in:\n")
        print(ds)

    return ds


def read_storm_settings(ds):
    settings = BoutOptions()

    if "filaments" in ds.options:
        # It's 3D
        if "realistic_geometry" in ds.options["filaments"]:
            settings["realistic_geometry"] = ds.options["filaments"][
                "realistic_geometry"
            ]
        else:
            settings["realistic_geometry"] = "slab"
    else:
        settings["realistic_geometry"] = "slab"

    storm_section = "filaments" if "filaments" in ds.options else "storm"

    tdim = ds.metadata["bout_tdim"]

    def _get_setting(ds, name, default):
        if name in ds:
            # Setting may be saved as a time-dependent variable
            value = ds[name]
            if tdim in value.dims:
                # Check the setting was not changed during the simulation
                # TODO: decide how to handle settings that did change
                value = value.isel({tdim: 0})
                if not (ds[name] == value).all():
                    raise ValueError(f"Setting {name} changed during the run")
            return value.load().item()
        elif name in ds.metadata:
            return ds.metadata[name]
        elif name in ds.options[storm_section]:
            return ds.options[storm_section][name]
        else:
            return default

    settings_dict = STORM_SETTINGS
    if storm_section == "filaments":
        settings_dict.update(STORM3D_SETTINGS)
    for name, default in settings_dict.items():
        settings[name] = _get_setting(ds, name, default)

    return settings


def calc_plasma_params(options):
    """
    Calculates a dictionary of important plasma parameters calculated from
    STORM input file options. Saved in SI units, except for temperatures
    which are in eV, and a set of dimensionless parameters.
    """

    # Set physical constants
    e = 1.602176565e-19  # Coulombs
    u = 1.66053892e-27  # kg
    m_e = 9.10938291e-31  # kg
    epsilon_0 = 8.854187817e-12  # Fm^-1
    mu_0 = 4 * np.pi * 1e-7  # Hm^-1

    attrs = {}

    # Determine if data is from storm2d or 3d
    if "filaments" in options.keys():
        opts_key = "filaments"
        attrs["bout_module"] = "STORM3D"
    elif "storm" in options.keys():
        opts_key = "storm"
        attrs["bout_module"] = "STORM2D"
    else:
        raise ValueError(
            "Cannot find options in input file data corresponding to storm2d or storm3d"
        )

    # Read plasma parameters
    module_opts = options[opts_key]
    Z = module_opts.evaluate_scalar("z")  # Atomic number
    B_0 = module_opts.evaluate_scalar("b_0")  # Tesla
    m_i = module_opts.evaluate_scalar("m_i") * u  # kg
    T_i = module_opts.evaluate_scalar("t_i0")  # eV
    T_e = module_opts.evaluate_scalar("t_e0")  # eV
    n_0 = module_opts.evaluate_scalar("n_0") * 100 * 100 * 100  # m^-3
    try:
        R_c = module_opts.evaluate_scalar("R_c")  # m
    except KeyError:
        # R_c may not always be present
        R_c = None
    try:
        q = module_opts.evaluate_scalar("q")
    except KeyError:
        # copy default from storm3d and storm2d
        q = 7.0

    # Calculate characteristic times/speeds etc.
    c_s = np.sqrt(T_e * e / m_i)  # ms^1
    Omega_i = Z * e * B_0 / m_i  # s^-1
    Omega_e = e * B_0 / m_e  # s^-1
    V_thi = np.sqrt(T_i * e / m_i)  # ms^-1
    V_the = np.sqrt(T_e * e / m_e)  # ms^-1
    rho_i = V_thi / Omega_i  # m
    rho_s = c_s / Omega_i  # m
    rho_e = V_the / Omega_e
    delta_t = 1.0 / Omega_i  # s

    def get_scalar(key, default):
        try:
            value = module_opts.evaluate_scalar(key)
        except KeyError:
            value = -1.0
        if value >= 0.0:
            return value
        else:
            return default

    loglambda = 18.0 - np.log(
        np.sqrt(n_0 / 1.0e19) * np.power(T_e / (1000.0), -1.5)
    )  # Coulomb logarithm
    loglambda = get_scalar("loglambda", default=loglambda)

    nu_ee = (
        n_0
        * (e**4)
        * loglambda
        / (np.sqrt(m_e) * (epsilon_0**2) * 3.0 * np.power(2 * np.pi * T_e * e, 1.5))
    )  # Hz
    nu_ei = (
        n_0
        * (Z**2)
        * (e**4)
        * loglambda
        / (np.sqrt(m_e) * (epsilon_0**2) * 3.0 * np.power(2 * np.pi * T_e * e, 1.5))
    )
    nu_ii = (
        n_0
        * ((Z * e) ** 4)
        * loglambda
        / (
            np.sqrt(m_i)
            * (epsilon_0**2)
            * 3.0
            * np.power(2 * np.pi * T_i * e, 1.5)
            * np.sqrt(2)
        )
    )

    params = {
        "Z": Z,
        "B_0": B_0,
        "m_i": m_i,
        "T_i0": T_i,
        "T_e0": T_e,
        "n_0": n_0,
        "c_s0": c_s,
        "omega_i0": Omega_i,
        "omega_e0": Omega_e,
        "V_thi0": V_thi,
        "V_the0": V_the,
        "rho_s0": rho_s,
        "rho_i0": rho_i,
        "delta_t": delta_t,
        "e": e,
        "u": u,
        "m_e": m_e,
        "epsilon_0": epsilon_0,
        "mu_0": mu_0,
        "loglambda": loglambda,
        "nu_ee0": nu_ee,
        "nu_ei0": nu_ei,
        "nu_ii0": nu_ii,
    }

    mu_n_base = (
        (1.0 + 1.3 * q**2) * (1.0 + T_i / T_e) * (rho_e**2) * nu_ei
    )  # Base Neo-Classical Particle diffusion, m^2s-1
    mu_n = mu_n_base / (rho_s * rho_s * Omega_i)  # Normalised
    mu_n = get_scalar("mu_n0", default=mu_n)

    mu_vort_base = (
        (1.0 + 1.6 * (q**2)) * (6.0 / 8.0) * (rho_i**2) * nu_ii
    )  # Base Neo-Classical Ion viscosity, m^2s-1
    mu_vort = mu_vort_base / (rho_s * rho_s * Omega_i)  # Normalised
    mu_vort = get_scalar("mu_vort0", default=mu_vort)

    if R_c is not None:
        params["R_c"] = R_c
        g = get_scalar("g0", default=2.0 * rho_s / R_c)
    else:
        g = get_scalar("g0", default=None)

    kappa_perp = (
        (1.0 + 1.6 * (q**2))
        * (4.66)
        * (rho_e**2)
        * nu_ee
        / (rho_s * rho_s * Omega_i)
    )
    kappa_perp = get_scalar("kappa0_perp", default=kappa_perp)

    kappa = 3.16 * ((V_the**2) / nu_ei) / (rho_s * rho_s * Omega_i)
    kappa = get_scalar("kappa0", default=kappa)

    if attrs["bout_module"] == "STORM2D":
        l = options["storm"].evaluate_scalar("L")
        lamb = rho_s / l

    mu = m_i / m_e
    eV_float_over_T = 0.5 * np.log(2 * np.pi / mu)

    dimless_params = {
        "mu_n0": mu_n,
        "mu_vort0": mu_vort,
        "kappa_perp0": kappa_perp,
        "kappa0": kappa,
        "eV_float_over_T": eV_float_over_T,
    }
    if g is not None:
        dimless_params["g0"] = g

    if attrs["bout_module"] == "STORM2D":
        dimless_params["lamb"] = lamb

    params["dimless"] = dimless_params

    return params


def apply_unnormalise(ds):
    """Convert the output data to SI units"""

    # TODO geometry is currently unused, should add coordinate normalisation for non-slab

    # TODO as BOUT v4.2 can now write attributes to the output NetCDF files
    # would it be better to calculate all normalisations within the BOUT++ module itself?

    if not ds.metadata["storm_normalised"]:
        warn("Dataset has already been unnormalised, doing nothing")
        return ds

    # Time treated separately because it needs renaming as well as unnormalising
    ds = add_time(ds)

    ds.metadata["storm_normalised"] = 0

    norms = calc_norms(ds)

    for data_var in norms.keys():
        # Unnormalise variables and coordinates
        if data_var in ds.variables or data_var in ds.coords:
            ds[data_var] = ds[data_var] * norms[data_var]["conversion"]
            ds[data_var].attrs.update(norms[data_var])
        # also unnormalise variables in metadata
        if data_var in ds.metadata:
            ds.metadata[data_var] = (
                ds.metadata[data_var] * norms[data_var]["conversion"]
            )

    return ds


def calc_norms(ds):
    """Create dictionary of normalisations of all possible variables of interest"""

    # Add physical normalisations
    params = ds.attrs["params"]
    e_field_norm = params["T_e0"] / params["rho_s0"]

    if "boussinesq" in ds.settings and (
        ds.settings["boussinesq"] == 1 or ds.settings["boussinesq"] == "true"
    ):
        # Original STORM Boussinesq approximation, divides 1/e*Div(J)=0 through by n
        vort_norm = (params["omega_i0"], "s^-1")
    else:
        # If boussinesq approximation not employed then vort is really more like a vorticity density
        vort_norm = (params["omega_i0"] * params["n_0"], "m^-3 s^-1")

    norms = {
        "n": {
            "conversion": params["n_0"],
            "units": "m^-3",
            "standard_name": "density",
            "long_name": "Ion number density",
        },
        "n_eq": {
            "conversion": params["n_0"],
            "units": "m^-3",
            "standard_name": "equilibrium density",
            "long_name": "Equilibrium ion number density",
        },
        "phi": {
            "conversion": params["T_e0"],
            "units": "V",
            "standard_name": "potential",
            "long_name": "Plasma potential",
        },
        "phi_eq": {
            "conversion": params["T_e0"],
            "units": "V",
            "standard_name": "equilibrium potential",
            "long_name": "Equilibrium plasma potential",
        },
        "vort": {
            "conversion": vort_norm[0],
            "units": vort_norm[1],
            "standard_name": "vorticity",
            "long_name": "Vorticity",
        },
        "vort_eq": {
            "conversion": vort_norm[0],
            "units": vort_norm[1],
            "standard_name": "equilibrium vorticity",
            "long_name": "Equilibrium vorticity",
        },
        "T": {
            "conversion": params["T_e0"],
            "units": "eV",
            "standard_name": "temperature",
            "long_name": "Electron temperature",
        },
        "T_eq": {
            "conversion": params["T_e0"],
            "units": "eV",
            "standard_name": "equilibrium temperature",
            "long_name": "Equilibrium electron temperature",
        },
        "T_i": {
            "conversion": params["T_e0"],
            "units": "eV",
            "standard_name": "ion temperature",
            "long_name": "Ion temperature",
        },
        "Ti_eq": {
            "conversion": params["T_e0"],
            "units": "eV",
            "standard_name": "equilibrium ion temperature",
            "long_name": "Equilibrium ion temperature",
        },
        "U": {
            "conversion": params["c_s0"],
            "units": "m s^-1",
            "standard_name": "ion parallel velocity",
            "long_name": "Ion parallel velocity",
        },
        "U_eq": {
            "conversion": params["c_s0"],
            "units": "m s^-1",
            "standard_name": "equilibrium ion parallel velocity",
            "long_name": "Equilibrium ion parallel velocity",
        },
        "V": {
            "conversion": params["c_s0"],
            "units": "m s^-1",
            "standard_name": "electron parallel velocity",
            "long_name": "Electron parallel velocity",
        },
        "V_eq": {
            "conversion": params["c_s0"],
            "units": "m s^-1",
            "standard_name": "equilibrium electron parallel velocity",
            "long_name": "Equilibrium electron parallel velocity",
        },
        "qpar": {
            "conversion": params["n_0"] * params["T_e0"] * params["e"] * params["c_s0"],
            "units": "W m^-2",
            "standard_name": "heat flux",
            "long_name": "Heat Flux",
        },
        "S": {
            "conversion": params["n_0"] / params["delta_t"],
            "units": "m^-3 s^-1",
            "standard_name": "particle source",
            "long_name": "Particle source",
        },
        "S_eq": {
            "conversion": params["n_0"] / params["delta_t"],
            "units": "m^-3 s^-1",
            "standard_name": "equilibrium particle source",
            "long_name": "Equilibrium particle source",
        },
        "S_E": {
            "conversion": params["n_0"]
            * params["e"]
            * params["T_e0"]
            / params["delta_t"],
            "units": "W m^-3",
            "standard_name": "energy source",
            "long_name": "Energy source",
        },
        "S_E_eq": {
            "conversion": params["n_0"]
            * params["e"]
            * params["T_e0"]
            / params["delta_t"],
            "units": "W m^-3",
            "standard_name": "equilibrium energy source",
            "long_name": "Equilibrium energy source",
        },
        "B": {
            "conversion": params["B_0"],
            "units": "T",
            "standard_name": "B-field",
            "long_name": "Magnetic field strength",
        },
        "E_z": {
            "conversion": e_field_norm,
            "units": "V m^-1",
            "standard_name": "binormal E-field",
            "long_name": "Binormal Electric field",
        },
        "t_array": {
            "conversion": params["delta_t"],
            "units": "s",
            "standard_name": "time",
            "long_name": "Simulation time elapsed",
        },
        "sigma_n": {
            "conversion": params["omega_i0"],
            "units": "s^-1",
            "standard_name": "parallel density loss rate",
            "long_name": "Parallel density loss rate",
        },
        "sigma_vort": {
            "conversion": params["omega_i0"],
            "units": "s^-1",
            "standard_name": "parallel vorticity loss rate",
            "long_name": "Parallel vorticity loss rate",
        },
        "sigma_T": {
            "conversion": params["omega_i0"],
            "units": "s^-1",
            "standard_name": "parallel temperature loss rate",
            "long_name": "Parallel temperature loss rate",
        },
        "loss": {
            "conversion": 1.0,
            "units": "None",
            "standard_name": "parallel loss",
            "long_name": "Parallel loss fraction",
        },
    }
    norms["n_aligned"] = norms["n"]
    norms["T_aligned"] = norms["T"]
    norms["U_aligned"] = norms["U"]
    norms["V_aligned"] = norms["V"]
    norms["vort_aligned"] = norms["vort"]
    norms["phi_aligned"] = norms["phi"]

    if ds.settings.get_bool("normalise_lengths", False):
        for g in ["g_11", "g_22", "g_33", "g_12", "g_13", "g_23"]:
            norms[g] = {
                "conversion": params["rho_s0"] ** 2,
                "units": "m^2",
                "standard_name": g,
                "long_name": f"Covariant metric component {g}",
            }
            norms[g + "_CELL_YLOW"] = norms[g]
        for g in ["g11", "g22", "g33", "g12", "g13", "g23"]:
            norms[g] = {
                "conversion": 1.0 / params["rho_s0"] ** 2,
                "units": "m^-2",
                "standard_name": g,
                "long_name": f"Contravariant metric component {g}",
            }
            norms[g + "_CELL_YLOW"] = norms[g]
        for G in ["G1", "G2", "G3"]:
            norms[g] = {
                "conversion": 1.0 / params["rho_s0"] ** 2,
                "units": "m^-2",
                "standard_name": G,
                "long_name": G,
            }
            norms[g + "_CELL_YLOW"] = norms[g]

        if ds.settings["realistic_geometry"] not in ("none", "slab"):
            warn(
                "not handling bxcv as normalise_lengths=true does not normalise bxcv - "
                "this may be incorrect!"
            )

    elif ds.settings.get_bool("normalise_all", False):
        norms["g_11"] = {
            "conversion": 1.0 / (params["rho_s0"] ** 2 * params["B_0"] ** 2),
            "units": "m^-2 T^-2",
            "standard_name": "g_11",
            "long_name": f"Covariant metric component g_11",
        }
        norms["g_11_CELL_YLOW"] = norms["g_11"]
        norms["g_22"] = {
            "conversion": params["rho_s0"] ** 2,
            "units": "m^2",
            "standard_name": "g_22",
            "long_name": f"Covariant metric component g_22",
        }
        norms["g_22_CELL_YLOW"] = norms["g_22"]
        norms["g_33"] = {
            "conversion": params["rho_s0"] ** 2,
            "units": "m^2",
            "standard_name": "g_33",
            "long_name": f"Covariant metric component g_33",
        }
        norms["g_33_CELL_YLOW"] = norms["g_33"]
        norms["g_12"] = {
            "conversion": 1.0 / params["B_0"],
            "units": "T^-1",
            "standard_name": "g_12",
            "long_name": f"Covariant metric component g_12",
        }
        norms["g_12_CELL_YLOW"] = norms["g_12"]
        norms["g_13"] = {
            "conversion": 1.0 / params["B_0"],
            "units": "T^-1",
            "standard_name": "g_13",
            "long_name": f"Covariant metric component g_13",
        }
        norms["g_13_CELL_YLOW"] = norms["g_13"]
        norms["g_23"] = {
            "conversion": params["rho_s0"] ** 2,
            "units": "m^2",
            "standard_name": "g_23",
            "long_name": f"Covariant metric component g_23",
        }
        norms["g_23_CELL_YLOW"] = norms["g_23"]
        norms["g11"] = {
            "conversion": params["rho_s0"] ** 2 * params["B_0"] ** 2,
            "units": "m^2 T^2",
            "standard_name": "g11",
            "long_name": f"Contravariant metric component g11",
        }
        norms["g11_CELL_YLOW"] = norms["g11"]
        norms["g22"] = {
            "conversion": 1.0 / params["rho_s0"] ** 2,
            "units": "m^-2",
            "standard_name": "g22",
            "long_name": f"Covariant metric component g22",
        }
        norms["g22_CELL_YLOW"] = norms["g22"]
        norms["g33"] = {
            "conversion": 1.0 / params["rho_s0"] ** 2,
            "units": "m^-2",
            "standard_name": "g33",
            "long_name": f"Contravariant metric component g33",
        }
        norms["g33_CELL_YLOW"] = norms["g33"]
        norms["g12"] = {
            "conversion": params["B_0"],
            "units": "T",
            "standard_name": "g12",
            "long_name": f"Contravariant metric component g12",
        }
        norms["g12_CELL_YLOW"] = norms["g12"]
        norms["g13"] = {
            "conversion": params["B_0"],
            "units": "T",
            "standard_name": "g13",
            "long_name": f"Contravariant metric component g13",
        }
        norms["g13_CELL_YLOW"] = norms["g13"]
        norms["g23"] = {
            "conversion": 1.0 / params["rho_s0"] ** 2,
            "units": "m^-2",
            "standard_name": "g23",
            "long_name": f"Contravariant metric component g23",
        }
        norms["g23_CELL_YLOW"] = norms["g23"]
        norms["dx"] = {
            "conversion": params["rho_s0"] ** 2 * params["B_0"],
            "units": "m^2 T",
            "standard_name": "dx",
            "long_name": "Grid spacing in poloidal magnetic flux",
        }
        norms["dx_CELL_YLOW"] = norms["dx"]
        norms["dy"] = {
            "conversion": 1.0,
            "units": "",
            "standard_name": "dy",
            "long_name": "dy",
        }
        norms["dy_CELL_YLOW"] = norms["dy"]
        norms["dz"] = {
            "conversion": 1.0,
            "units": "",
            "standard_name": "dz",
            "long_name": "dz",
        }
        norms["dz_CELL_YLOW"] = norms["dz"]
        norms["J"] = {
            "conversion": params["rho_s0"] / params["B_0"],
            "units": "T^-1 m",
            "standard_name": "Jacobian",
            "long_name": "Jacobian of field-aligned coordinate system",
        }
        norms["J_CELL_YLOW"] = norms["J"]
        norms["G1"] = {
            "conversion": params["B_0"],
            "units": "T",
            "standard_name": "G1",
            "long_name": "G1",
        }
        norms["G1_CELL_YLOW"] = norms["G1"]
        norms["G2"] = {
            "conversion": 1.0 / params["rho_s0"] ** 2,
            "units": "m^-2",
            "standard_name": "G2",
            "long_name": "G2",
        }
        norms["G2_CELL_YLOW"] = norms["G2"]
        norms["G3"] = {
            "conversion": 1.0 / params["rho_s0"] ** 2,
            "units": "m^-2",
            "standard_name": "G3",
            "long_name": "G3",
        }
        norms["G3_CELL_YLOW"] = norms["G3"]

        # psi_poloidal is read from the grid file, so does not need conversion, but does
        # have units in the 'normalise_all=true' conventions
        norms["psi_poloidal"] = {
            "conversion": 1.0,
            "units": "T m^2",
            "standard_name": "psi",
            "long_name": "Poloidal magnetic flux function",
        }

        # B/2*Curl(b/B).Grad(x) has units T
        # as x is poloidal magnetic flux, has units [T m^2]
        norms["bxcvx"] = {
            "conversion": params["B_0"],
            "units": "T",
            "standard_name": "bxcvx",
            "long_name": "x-component of B/2*Curl(b/B)",
        }
        # B/2*Curl(b/B).Grad(y) has units m^-2
        # as y is dimensionless
        norms["bxcvy"] = {
            "conversion": 1.0 / params["rho_s0"] ** 2,
            "units": "m^-2",
            "standard_name": "bxcvy",
            "long_name": "y-component of B/2*Curl(b/B)",
        }
        # B/2*Curl(b/B).Grad(z) has units m^-2
        # as z is dimensionless
        norms["bxcvz"] = {
            "conversion": 1.0 / params["rho_s0"] ** 2,
            "units": "m^-2",
            "standard_name": "bxcvz",
            "long_name": "z-component of B/2*Curl(b/B)",
        }

    if ds.settings["realistic_geometry"] == "slab":
        norms["radial"] = {
            "conversion": params["rho_s0"],
            "units": "m",
            "standard_name": "radial",
            "long_name": "Radial",
        }
        norms["parallel"] = {
            "conversion": params["rho_s0"],
            "units": "m",
            "standard_name": "parallel",
            "long_name": "Parallel",
        }
        norms["binormal"] = {
            "conversion": params["rho_s0"],
            "units": "m",
            "standard_name": "binormal",
            "long_name": "Binormal",
        }

        norms["dx"] = {
            "conversion": params["rho_s0"],
            "units": "m",
            "standard_name": "dx",
            "long_name": "dx",
        }
        norms["dy"] = {
            "conversion": params["rho_s0"],
            "units": "m",
            "standard_name": "dy",
            "long_name": "dy",
        }
        norms["dz"] = {
            "conversion": params["rho_s0"],
            "units": "m",
            "standard_name": "dz",
            "long_name": "dz",
        }

        if ds.settings.get_bool("normalise_lengths", False):
            warn(
                "slab simulation run with normalise_lengths=true, which is not expected"
            )
        if ds.settings.get_bool("normalise_all", False):
            warn("slab simulation run with normalise_all=true, which is not expected")

    return norms


def add_normalised_units(ds):
    """Convert the output data to SI units"""

    # Use values from 't_array' as the time coordinate
    if "t_array" in ds:
        ds = ds.rename(t_array="t")
    ds = ds.rename(t="time")

    attrs = create_normalised_attrs(ds)

    for data_var, var_attrs in attrs.items():
        if data_var in ds.variables or data_var in ds.coords:
            ds[data_var].attrs.update(var_attrs)

    return ds


def create_normalised_attrs(ds):
    """Add units of normalised variables"""

    if "boussinesq" in ds.settings and (
        ds.settings["boussinesq"] == 1 or ds.settings["boussinesq"] == "true"
    ):
        # Original STORM Boussinesq approximation, divides 1/e*Div(J)=0 through by n
        vort_units = "omega_i0"
    else:
        # If boussinesq approximation not employed then vort is really more like a vorticity density
        vort_units = "omega_i0.n_0"

    attrs = {
        "n": {
            "units": "n_0",
            "standard_name": "density",
            "long_name": "Ion number density",
        },
        "n_eq": {
            "units": "n_0",
            "standard_name": "equilibrium density",
            "long_name": "Equilibrium ion number density",
        },
        "phi": {
            "units": "T_e0.e^-1",
            "standard_name": "potential",
            "long_name": "Plasma potential",
        },
        "phi_eq": {
            "units": "T_e0.e^-1",
            "standard_name": "equilibrium potential",
            "long_name": "Equilibrium plasma potential",
        },
        "vort": {
            "units": vort_units,
            "standard_name": "vorticity",
            "long_name": "Vorticity",
        },
        "vort_eq": {
            "units": vort_units,
            "standard_name": "equilibrium vorticity",
            "long_name": "Equilibrium vorticity",
        },
        "T": {
            "units": "T_e0",
            "standard_name": "temperature",
            "long_name": "Electron temperature",
        },
        "T_eq": {
            "units": "T_e0",
            "standard_name": "equilibrium temperature",
            "long_name": "Equilibrium electron temperature",
        },
        "T_i": {
            "units": "T_e0",
            "standard_name": "ion temperature",
            "long_name": "Ion temperature",
        },
        "Ti_eq": {
            "units": "T_e0",
            "standard_name": "equilibrium ion temperature",
            "long_name": "Equilibrium ion temperature",
        },
        "U": {
            "units": "c_s0",
            "standard_name": "ion parallel velocity",
            "long_name": "Ion parallel velocity",
        },
        "U_eq": {
            "units": "c_s0",
            "standard_name": "equilibrium ion parallel velocity",
            "long_name": "Equilibrium ion parallel velocity",
        },
        "V": {
            "units": "c_s0",
            "standard_name": "electron parallel velocity",
            "long_name": "Electron parallel velocity",
        },
        "V_eq": {
            "units": "c_s0",
            "standard_name": "equilibrium electron parallel velocity",
            "long_name": "Equilibrium electron parallel velocity",
        },
        "qpar": {
            "units": "n_0.T_e0.c_s0",
            "standard_name": "heat flux",
            "long_name": "Heat Flux",
        },
        "S": {
            "units": "n_0.Omega_i0",
            "standard_name": "particle source",
            "long_name": "Particle source",
        },
        "S_eq": {
            "units": "n_0.Omega_i0",
            "standard_name": "equilibrium particle source",
            "long_name": "Equilibrium particle source",
        },
        "S_E": {
            "units": "n_0.T_e0.Omega_i0",
            "standard_name": "energy source",
            "long_name": "Energy source",
        },
        "S_E_eq": {
            "units": "n_0.T_e0.Omega_i0",
            "standard_name": "equilibrium energy source",
            "long_name": "Equilibrium energy source",
        },
        "B": {
            "units": "B_0",
            "standard_name": "B-field",
            "long_name": "Magnetic field strength",
        },
        "E_z": {
            "units": "T_e0.e^-1.rho_s0^-1",
            "standard_name": "binormal E-field",
            "long_name": "Binormal Electric field",
        },
        "t_array": {
            "units": "Omega_i0^-1",
            "standard_name": "time",
            "long_name": "Simulation time elapsed",
        },
        "sigma_n": {
            "units": "Omega_i0",
            "standard_name": "parallel density loss rate",
            "long_name": "Parallel density loss rate",
        },
        "sigma_vort": {
            "units": "Omega_i0",
            "standard_name": "parallel vorticity loss rate",
            "long_name": "Parallel vorticity loss rate",
        },
        "sigma_T": {
            "units": "Omega_i0",
            "standard_name": "parallel temperature loss rate",
            "long_name": "Parallel temperature loss rate",
        },
        "loss": {
            "units": "None",
            "standard_name": "parallel loss",
            "long_name": "Parallel loss fraction",
        },
    }
    attrs["n_aligned"] = attrs["n"]
    attrs["T_aligned"] = attrs["T"]
    attrs["U_aligned"] = attrs["U"]
    attrs["V_aligned"] = attrs["V"]
    attrs["vort_aligned"] = attrs["vort"]
    attrs["phi_aligned"] = attrs["phi"]

    if ds.settings.get_bool("normalise_lengths", False):
        for g in ["g_11", "g_22", "g_33", "g_12", "g_13", "g_23"]:
            attrs[g] = {
                "units": "rho_s0^2",
                "standard_name": g,
                "long_name": f"Covariant metric component {g}",
            }
            attrs[g + "_CELL_YLOW"] = attrs[g]
        for g in ["g11", "g22", "g33", "g12", "g13", "g23"]:
            attrs[g] = {
                "units": "rho_s0^-2",
                "standard_name": g,
                "long_name": f"Contravariant metric component {g}",
            }
            attrs[g + "_CELL_YLOW"] = attrs[g]
        for G in ["G1", "G2", "G3"]:
            attrs[g] = {"units": "rho_s0^-2", "standard_name": G, "long_name": G}
            attrs[g + "_CELL_YLOW"] = attrs[g]

        if ds.settings["realistic_geometry"] not in ("none", "slab"):
            warn(
                "not handling bxcv as normalise_lengths=true does not normalise bxcv - "
                "this may be incorrect!"
            )

    elif ds.settings.get_bool("normalise_all", False):
        attrs["g_11"] = {
            "units": "rho_s0^-2.B_0^-2",
            "standard_name": "g_11",
            "long_name": f"Covariant metric component g_11",
        }
        attrs["g_11_CELL_YLOW"] = attrs["g_11"]
        attrs["g_22"] = {
            "units": "rho_s0^2",
            "standard_name": "g_22",
            "long_name": f"Covariant metric component g_22",
        }
        attrs["g_22_CELL_YLOW"] = attrs["g_22"]
        attrs["g_33"] = {
            "units": "rho_s0^2",
            "standard_name": "g_33",
            "long_name": f"Covariant metric component g_33",
        }
        attrs["g_33_CELL_YLOW"] = attrs["g_33"]
        attrs["g_12"] = {
            "units": "B_0^-1",
            "standard_name": "g_12",
            "long_name": f"Covariant metric component g_12",
        }
        attrs["g_12_CELL_YLOW"] = attrs["g_12"]
        attrs["g_13"] = {
            "units": "B_0^-1",
            "standard_name": "g_13",
            "long_name": f"Covariant metric component g_13",
        }
        attrs["g_13_CELL_YLOW"] = attrs["g_13"]
        attrs["g_23"] = {
            "units": "rho_s0^2",
            "standard_name": "g_23",
            "long_name": f"Covariant metric component g_23",
        }
        attrs["g_23_CELL_YLOW"] = attrs["g_23"]
        attrs["g11"] = {
            "units": "rho_s0^2.B_0^2",
            "standard_name": "g11",
            "long_name": f"Contravariant metric component g11",
        }
        attrs["g11_CELL_YLOW"] = attrs["g11"]
        attrs["g22"] = {
            "units": "rho_s0^-2",
            "standard_name": "g22",
            "long_name": f"Covariant metric component g22",
        }
        attrs["g22_CELL_YLOW"] = attrs["g22"]
        attrs["g33"] = {
            "units": "rho_s0^-2",
            "standard_name": "g33",
            "long_name": f"Contravariant metric component g33",
        }
        attrs["g33_CELL_YLOW"] = attrs["g33"]
        attrs["g12"] = {
            "units": "B_0",
            "standard_name": "g12",
            "long_name": f"Contravariant metric component g12",
        }
        attrs["g12_CELL_YLOW"] = attrs["g12"]
        attrs["g13"] = {
            "units": "B_0",
            "standard_name": "g13",
            "long_name": f"Contravariant metric component g13",
        }
        attrs["g13_CELL_YLOW"] = attrs["g13"]
        attrs["g23"] = {
            "units": "rho_s0^-2",
            "standard_name": "g23",
            "long_name": f"Contravariant metric component g23",
        }
        attrs["g23_CELL_YLOW"] = attrs["g23"]
        attrs["dx"] = {
            "units": "rho_s0^2.B_0",
            "standard_name": "dx",
            "long_name": "Grid spacing in poloidal magnetic flux",
        }
        attrs["dx_CELL_YLOW"] = attrs["dx"]
        attrs["J"] = {
            "units": "rho_s0.B_0^-1",
            "standard_name": "Jacobian",
            "long_name": "Jacobian of field-aligned coordinate system",
        }
        attrs["J_CELL_YLOW"] = attrs["J"]
        attrs["G1"] = {"units": "B_0", "standard_name": "G1", "long_name": "G1"}
        attrs["G1_CELL_YLOW"] = attrs["G1"]
        attrs["G2"] = {"units": "rho_s0^-2", "standard_name": "G2", "long_name": "G2"}
        attrs["G2_CELL_YLOW"] = attrs["G2"]
        attrs["G3"] = {"units": "rho_s0^-2", "standard_name": "G3", "long_name": "G3"}
        attrs["G3_CELL_YLOW"] = attrs["G3"]

        # psi_poloidal is read from the grid file, so does not need conversion, but does
        # have units in the 'normalise_all=true' conventions
        attrs["psi_poloidal"] = {
            "units": "T m^2",
            "standard_name": "psi",
            "long_name": "Poloidal magnetic flux function",
        }

        attrs["bxcvx"] = {
            "units": "B_0",
            "standard_name": "bxcvx",
            "long_name": "x-component of B/2*Curl(b/B)",
        }
        attrs["bxcvy"] = {
            "units": "rho_s0^-2",
            "standard_name": "bxcvy",
            "long_name": "y-component of B/2*Curl(b/B)",
        }
        attrs["bxcvz"] = {
            "units": "rho_s0^-2",
            "standard_name": "bxcvz",
            "long_name": "z-component of B/2*Curl(b/B)",
        }

    if ds.settings["realistic_geometry"] == "slab":
        attrs["radial"] = {
            "units": "rho_s0",
            "standard_name": "radial",
            "long_name": "Radial",
        }
        attrs["parallel"] = {
            "units": "rho_s0",
            "standard_name": "parallel",
            "long_name": "Parallel",
        }
        attrs["binormal"] = {
            "units": "rho_s0",
            "standard_name": "binormal",
            "long_name": "Binormal",
        }

        attrs["dx"] = {
            "units": "rho_s0",
            "standard_name": "dx",
            "long_name": "dx",
        }
        attrs["dy"] = {
            "units": "rho_s0",
            "standard_name": "dy",
            "long_name": "dy",
        }
        attrs["dz"] = {
            "units": "rho_s0",
            "standard_name": "dz",
            "long_name": "dz",
        }

        if ds.settings.get_bool("normalise_lengths", False):
            warn(
                "slab simulation run with normalise_lengths=true, which is not expected"
            )
        if ds.settings.get_bool("normalise_all", False):
            warn("slab simulation run with normalise_all=true, which is not expected")

    return attrs
