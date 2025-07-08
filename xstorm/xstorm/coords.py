import numpy as np
from xbout.utils import _add_attrs_to_var, _1d_coord_from_spacing


# TODO generalise to non-uniform and non-slab geometries


def add_time(ds):
    """Add time coordinate in physical units"""

    # Time dimension
    delta_t = ds.attrs["params"]["delta_t"]
    if "t_array" in ds:
        t_normalised = ds["t_array"]
    else:
        t_normalised = ds.coords["t"]
    time = t_normalised.values * delta_t
    ds = ds.rename({"t": "time"})
    ds = ds.assign_coords(time=("time", time))
    ds["time"].attrs["units"] = "s"
    ds["time"].attrs["standard_name"] = "time"
    ds["time"].attrs["long_name"] = "Time"
    ds.metadata["bout_tdim"] = "time"
    ds["time"].attrs["options"] = ds.options
    ds["time"].attrs["params"] = ds.params
    if "geometry" in ds.attrs:
        ds["time"].attrs["geometry"] = ds.geometry
    if "regions" in ds.attrs:
        ds["time"].attrs["regions"] = ds.regions

    return ds


def add_slab_coords(ds):
    """
    Add coordinates in normalised units.

    Coordinates will be unnormalised if apply_unnormalise is called.
    """

    # Get mesh lengths from input file
    opts = ds.attrs["options"]
    mesh = opts["mesh"]
    storm_opts = opts["storm"] if "storm" in opts else opts["filaments"]
    if "lx" in mesh:
        l_x = mesh.evaluate_scalar("lx")
    else:
        l_x = storm_opts.evaluate_scalar("lx")
    if "ly" in mesh:
        l_y = mesh.evaluate_scalar("ly")
    else:
        l_y = storm_opts.evaluate_scalar("ly")
    if "lz" in mesh:
        l_z = mesh.evaluate_scalar("lz")
    elif "lz" in storm_opts:
        l_z = storm_opts.evaluate_scalar("lz")
    else:
        z_max = opts.evaluate_scalar("zmax")
        z_min = opts.evaluate_scalar("zmin")
        l_z = (z_max - z_min) * 2 * np.pi

    # Radial dimension
    n_x = ds.dims["x"]
    ds = ds.rename({"x": "radial"})
    radial = _1d_coord_from_spacing(ds["dx"], "radial")
    ds["radial"] = radial
    ds["radial"].attrs["units"] = "m"
    ds["radial"].attrs["standard_name"] = "radial distance"
    ds["radial"].attrs["long_name"] = "Radial Distance"
    ds.metadata["bout_xdim"] = "radial"
    _add_attrs_to_var(ds, "radial")

    # Parallel dimension
    n_y = ds.dims["y"]
    ds = ds.rename({"y": "parallel"})
    parallel = _1d_coord_from_spacing(ds["dy"], "parallel", origin_at="centre")
    ds["parallel"] = parallel
    ds["parallel"].attrs["units"] = "m"
    ds["parallel"].attrs["standard_name"] = "parallel distance"
    ds["parallel"].attrs["long_name"] = "Parallel Distance"
    ds.metadata["bout_ydim"] = "parallel"
    _add_attrs_to_var(ds, "parallel")

    # Binormal dimension
    n_z = ds.dims["z"]
    ds = ds.rename({"z": "binormal"})
    binormal = np.linspace(start=0.0, stop=l_z, num=n_z, endpoint=False)
    ds["binormal"] = ("binormal", binormal)
    ds["binormal"].attrs["units"] = "m"
    ds["binormal"].attrs["standard_name"] = "binormal distance"
    ds["binormal"].attrs["long_name"] = "Binormal Distance"
    ds.metadata["bout_zdim"] = "binormal"
    _add_attrs_to_var(ds, "binormal")

    # Include params attr on all coords
    for coord in ds.coords:
        ds[coord].attrs["params"] = ds.params

    # Add zShift as a coordinate, so that it gets interpolated along with a variable
    try:
        ds = ds.set_coords("zShift")
    except ValueError:
        # zShift is 0.0 for slab geometry anyway
        ds = ds.assign_coords(zShift=ds["dx"].copy())
        ds["zShift"].values[:] = 0.0
    try:
        ds = ds.set_coords("zShift_CELL_XLOW")
    except ValueError:
        pass
    try:
        ds = ds.set_coords("zShift_CELL_YLOW")
    except ValueError:
        pass
    try:
        ds = ds.set_coords("zShift_CELL_ZLOW")
    except ValueError:
        pass

    return ds
