import re

import xarray as xr

from xbout.load import _expand_filepaths, _add_options

from .load import (
    read_storm_settings,
    calc_plasma_params,
    apply_unnormalise,
    add_normalised_units,
)
from .coords import add_time, add_slab_coords


def open_raw_fastoutput(datapath):
    """
    Opens fast output data and combines into a single dataset.

    Returns data still using simulation normalisation.
    """

    # Get list of all files
    filepaths, filetype = _expand_filepaths(datapath)

    # Iterate over all files, extracting DataArrays ready for combining
    fo_data = []
    for i, filepath in enumerate(filepaths):

        fo = xr.open_dataset(filepath)

        if i == 0:
            # Get time coordinate from first file
            time = fo["time"]

        # Time is global, and we already extracted it
        fo = fo.drop_vars("time", errors="ignore")

        # There might be no virtual probe in this region
        if len(fo.data_vars) > 0:

            for name, da in fo.items():

                # Save the physical position (in index units)
                da = da.expand_dims(x=1, y=1, z=1)
                da = da.assign_coords(
                    x=xr.DataArray([da.attrs["ix"]], dims=["x"]),
                    y=xr.DataArray([da.attrs["iy"]], dims=["y"]),
                    z=xr.DataArray([da.attrs["iz"]], dims=["z"]),
                )

                # Re-attach the time coordinate
                da = da.assign_coords(time=time)

                # We saved the position, so don't care what number the variable was
                # Only need it's name (i.e. n, T, etc.)
                regex = re.compile("(\D+)([0-9]+)")
                match = regex.match(name)
                if match is None:
                    raise ValueError(f"Regex could not parse the variable named {name}")
                var, num = match.groups()
                da.name = var

                # Must promote DataArrays to Datasets until xarray GH #3248 is fixed
                ds = xr.Dataset({var: da})
                fo_data.append(ds)

        fo.close()

    # This will merge different variables, and arrange by physical position
    full_fo = xr.combine_by_coords(fo_data)

    return full_fo


def open_fastoutput(datapath="BOUT.fast.*.nc", inputfilepath=None, unnormalise=True):
    """
    Load a time series dataset from a set of STORM fast output files,
    including the input options file.

    Parameters
    ----------
    datapath: str, default "BOUT.fast.*.nc"
        Paths to files to load FastOutput data from
    inputfilepath: str, optional
        Input file to read parameters needed for unnormalising from
    unnormalise: bool, default True
        By default convert the output to SI units, set to False to leave normalised
    """

    fo = open_raw_fastoutput(datapath)

    fo = _add_options(fo, inputfilepath)

    if unnormalise and not fo.options:
        raise ValueError(
            "Can't calculate Storm normalisations without the " "input file options"
        )

    if "metadata" not in fo.attrs:
        fo.attrs["metadata"] = {}

    # Read various STORM options from the input file or output files
    settings = read_storm_settings(fo)
    fo.attrs["settings"] = settings

    # Calculate and store STORM-specific plasma parameters
    params = calc_plasma_params(fo.options)
    fo.attrs["params"] = params
    for var in fo:
        fo[var].attrs["params"] = params

    # Add physical coordinates
    fo = fo.rename({"time": "t_array"})
    fo = add_time(fo)
    fo = fo.drop("t_array")
    if fo.settings["realistic_geometry"] == "slab":
        fo = add_slab_coords(fo)
    else:
        # theta and zeta are angle or angle-like coordinates so do not need any units or
        # conversion.
        # psi_poloidal is read in physical units from the grid file, so also does not need
        # any conversion here.
        pass

    if unnormalise:
        fo = apply_unnormalise(fo)
    else:
        fo = add_normalised_units(fo)

    return fo
