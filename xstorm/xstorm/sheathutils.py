# Utility methods for SheathDataset and SheathDataArray accessors

import numpy as np


def apply_to_all_regions_lower(da, method, **kwargs):
    result = {}
    for region_name, region in da.regions.items():
        if region.connection_lower_y is None:
            # region has lower y boundary
            result[region_name] = method(region_name, **kwargs)
    if len(result) == 1:
        # Just return the single value instead of dict with one entry
        return next(iter(result.values()))
    else:
        return result


def apply_to_all_regions_upper(da, method, **kwargs):
    result = {}
    for region_name, region in da.regions.items():
        if region.connection_upper_y is None:
            # region has upper y boundary
            result[region_name] = method(region_name, **kwargs)
    if len(result) == 1:
        # Just return the single value instead of dict with one entry
        return next(iter(result.values()))
    else:
        return result


def ambipolar_potential_lowerupper(
    sheathaccessor, region=None, lowerupper=None, *, return_aligned=False
):
    """
    Calculate ambipolar sheath potential at sheaths

    Parameters
    ----------
    lowerupper : {"lower", "upper"}
        Set to calculate the ambipolar sheath potential at lower or upper sheath
    region : str, optional
        Get the sheath value for a specific region. By default looks for all
        regions.
    return_aligned : bool, default False
        Return result on field-aligned grid (note, if input is field-aligned already
        then the result is always on the field-aligned grid)
    """
    if sheathaccessor.data.metadata["storm_normalised"]:
        raise ValueError(
            "ambipolar_potential_lowerupper() is only implemented for unnormalised "
            "Datasets so far"
        )

    if region is None:
        result = {}
        if lowerupper == "lower":
            for region_name, region in sheathaccessor.data.regions.items():
                if region.connection_lower_y is None:
                    # region has lower y boundary
                    result[region_name] = ambipolar_potential_lowerupper(
                        sheathaccessor,
                        region_name,
                        lowerupper,
                        return_aligned=return_aligned,
                    )
        elif lowerupper == "upper":
            for region_name, region in sheathaccessor.data.regions.items():
                if region.connection_upper_y is None:
                    # region has upper y boundary
                    result[region_name] = ambipolar_potential_lowerupper(
                        sheathaccessor,
                        region_name,
                        lowerupper,
                        return_aligned=return_aligned,
                    )
        else:
            raise ValueError(
                f"lowerupper must be 'lower' or 'upper'. Got '{lowerupper}'"
            )

        if len(result) == 1:
            # Just return the single value instead of dict with one entry
            return next(iter(result.values()))
        else:
            return result

    ds = sheathaccessor.data

    if lowerupper == "lower":
        T_e = sheathaccessor.T_e_lower(region, return_aligned=return_aligned)
        T_i = sheathaccessor.T_i_lower(region, return_aligned=return_aligned)
    elif lowerupper == "upper":
        T_e = sheathaccessor.T_e_upper(region, return_aligned=return_aligned)
        T_i = sheathaccessor.T_i_upper(region, return_aligned=return_aligned)
    else:
        raise ValueError(f"lowerupper must be 'lower' or 'upper'. Got '{lowerupper}'")

    m_e = ds.params["m_e"]
    m_i = ds.params["m_i"]
    # From Stangeby's (2.60). We have T_e in eV, so T_e=k*T_e[Kelvin]/e.
    # Stangeby's expression is for the wall potential when the sheath edge potential
    # is taken as zero, so we need the negative (which is the sheath potential when
    # the wall potential is taken as zero).
    # STORM's sheath BCs include some m_e/m_i correction to sound speed, but these
    # cancel in setting U_sheath=V_sheath to get ambipolar sheath potential
    phi_ambipolar = -T_e * 0.5 * np.log((2.0 * np.pi * m_e / m_i) * (1.0 + T_i / T_e))
    phi_ambipolar.name = "phi_ambipolar"
    phi_ambipolar.attrs["name"] = "phi_ambipolar"
    phi_ambipolar.attrs["long_name"] = "Ambipolar phi_sheath"
    return phi_ambipolar


def lower_from_field_aligned(sheathaccessor, sheath_slice):
    ycoord = sheathaccessor.metadata["bout_ydim"]
    zcoord = sheathaccessor.metadata["bout_zdim"]

    if zcoord not in sheath_slice.dims:
        # no toroidal/binormal dimension so nothing to do
        return sheath_slice

    if not sheath_slice.direction_y == "Aligned":
        raise ValueError(
            "Trying to transform non-aligned field from field-aligned coordinates"
        )

    if (
        sheathaccessor.data.cell_location == "CELL_CENTRE"
        or sheathaccessor.data.cell_location == "CELL_YLOW"
        or sheathaccessor.data.cell_location == "CELL_ZLOW"
    ):
        region = next(iter(sheath_slice.regions.keys()))
        myg = (
            sheathaccessor.metadata["MYG"]
            if sheathaccessor.metadata["keep_yboundaries"]
            else 0
        )
        zShift = (
            sheathaccessor.data["zShift_CELL_YLOW"]
            .bout.from_region(region, with_guards=0)
            .isel({ycoord: myg})
        )
    elif da.cell_location == "CELL_XLOW":
        raise ValueError("CELL_XLOW fields not supported yet")
    else:
        raise ValueError(f"Unsupported location {da.cell_location}.")

    result = sheath_slice.bout._shift_z(-zShift)
    result.attrs["direction_y"] = "Standard"

    return result


def upper_from_field_aligned(sheathaccessor, sheath_slice):
    ycoord = sheathaccessor.metadata["bout_ydim"]
    zcoord = sheathaccessor.metadata["bout_zdim"]

    if zcoord not in sheath_slice.dims:
        # no toroidal/binormal dimension so nothing to do
        return sheath_slice

    if not sheath_slice.direction_y == "Aligned":
        raise ValueError(
            "Trying to transform non-aligned field from field-aligned coordinates"
        )

    if (
        sheathaccessor.data.cell_location == "CELL_CENTRE"
        or sheathaccessor.data.cell_location == "CELL_YLOW"
        or sheathaccessor.data.cell_location == "CELL_ZLOW"
    ):
        myg = sheathaccessor.metadata["MYG"]
        region = next(iter(sheath_slice.regions))
        if not (sheathaccessor.metadata["keep_yboundaries"] and myg > 0):
            raise ValueError(
                "Missing upper sheath-edge value of zShift_CELL_YLOW because "
                "keep_yboundaries==False"
            )
        zShift = (
            sheathaccessor.data["zShift_CELL_YLOW"]
            .bout.from_region(region, with_guards=0)
            .isel({ycoord: -myg})
        )
    elif sheathaccessor.data.cell_location == "CELL_XLOW":
        raise ValueError("CELL_XLOW fields not supported yet")
    else:
        raise ValueError(f"Unsupported location {sheathaccessor.data.cell_location}.")

    result = sheath_slice.bout._shift_z(-zShift)
    result.attrs["direction_y"] = "Standard"

    return result
