def _update_name(da, *, prefix="", suffix="", replace=""):
    old_name = da.attrs["standard_name"]
    if replace:
        da.name = replace
        da.attrs["standard_name"] = replace
        da.attrs["long_name"] = replace
    else:
        da.name = prefix + old_name + suffix
        da.attrs["standard_name"] = prefix + old_name + suffix
        da.attrs["long_name"] = prefix + old_name + suffix
    return da


def _update_units(da, *, prefix="", suffix="", replace=""):
    old_units = da.attrs["units"]
    if replace:
        da.attrs["units"] = replace
    else:
        da.attrs["units"] = prefix + old_units + suffix
    return da


def get_dim_from_coord(coord, da):
    if coord in da.dims:
        return coord
    else:
        dim, *extra = da.coords[coord].dims

        if extra is None:
            return dim
        else:
            raise ValueError(
                f"Coordinate {coord} has multiple dimensions: " f"{dim, *extra}"
            )
