from matplotlib import pyplot as plt


def plot_sheath_currents(
    sheathaccessor, *, output="total", lowerupper="both", ax=None, **kwargs
):
    # documented in SheathDataset.plot_currents

    ds = sheathaccessor.data
    if ds.settings["realistic_geometry"] != "slab":
        # Would be nice to have an optional argument to current_lower() and
        # current_upper() like 'projection' which could be None for the current
        # behaviour and 'wall-normal' to get the component normal to the wall, which
        # we could then integrate below
        raise ValueError(
            "Non-slab geometry not supported yet. Need to handle projection of "
            "parallel current onto targets when B is not normal to the targets."
        )

    if lowerupper == "both":
        include_lower = True
        include_upper = True
    elif lowerupper == "upper":
        include_lower = False
        include_upper = True
    elif lowerupper == "lower":
        include_lower = True
        include_upper = False
    else:
        raise ValueError(f"Unrecognised value '{lowerupper}' for lowerupper argument")

    lower_currents = {}
    upper_currents = {}
    for region_name, region in ds.regions.items():
        if include_lower and region.connection_lower_y is None:
            # flip sign so J is positive outward
            current = sheathaccessor.current_lower(region_name)
            lower_currents[region_name] = -current
            lower_currents[region_name].attrs = current.attrs
        if include_upper and region.connection_upper_y is None:
            upper_currents[region_name] = sheathaccessor.current_upper(region_name)

    if output == "total":
        if ax is None:
            _, ax = plt.subplots()

        # Convert to total currents through each sheath
        for region_name in lower_currents:
            lower_currents[region_name] = lower_currents[region_name].integrate(
                dim=["radial", "binormal"]
            )
        for region_name in upper_currents:
            upper_currents[region_name] = upper_currents[region_name].integrate(
                dim=["radial", "binormal"]
            )
        sum_lower = sum(lower_currents.values())
        sum_upper = sum(upper_currents.values())
        sum_total = sum_upper - sum_lower

        ax.set_title("Integrated currents through sheaths")

        sum_total.plot(ls="-", c="k", ax=ax, label="Total current")
        if len(lower_currents) > 1:
            sum_lower.plot(ax=ax, label="Lower sheaths total")
        if len(upper_currents) > 1:
            sum_upper.plot(ax=ax, label="Upper sheaths total")

        for region_name, current in lower_currents.items():
            current.plot(ls=":", ax=ax, label=f"{region_name} lower")

        for region_name, current in upper_currents.items():
            current.plot(ls="--", ax=ax, label=f"{region_name} upper")

        ax.set_ylabel("Parallel current / A")

        l = ax.legend()
        l.set_draggable(True)

        return ax

    elif output == "positive_negative":
        if ax is None:
            _, ax = plt.subplots()

        # Integrate positive currents and negative currents through each sheath
        lower_positive_currents = {}
        lower_negative_currents = {}
        for region_name in lower_currents:
            J = lower_currents[region_name]
            positive = J.where(J > 0.0, 0.0)
            lower_positive_currents[region_name] = positive.integrate(
                dim=["radial", "binormal"]
            )
            negative = J.where(J < 0.0, 0.0)
            lower_negative_currents[region_name] = negative.integrate(
                dim=["radial", "binormal"]
            )

        upper_positive_currents = {}
        upper_negative_currents = {}
        for region_name in upper_currents:
            J = upper_currents[region_name]
            positive = J.where(J > 0.0, 0.0)
            upper_positive_currents[region_name] = positive.integrate(
                dim=["radial", "binormal"]
            )
            negative = J.where(J < 0.0, 0.0)
            upper_negative_currents[region_name] = negative.integrate(
                dim=["radial", "binormal"]
            )

        sum_lower_positive = sum(lower_positive_currents.values())
        sum_lower_negative = sum(lower_negative_currents.values())
        sum_upper_positive = sum(upper_positive_currents.values())
        sum_upper_negative = sum(upper_negative_currents.values())
        total = (
            sum_upper_positive
            - sum_upper_negative
            - sum_lower_positive
            + sum_lower_negative
        )

        ax.set_title("Integrated positive and negative currents through sheaths")

        prop_cycle = iter(plt.rcParams["axes.prop_cycle"])

        total.plot(ls="-", c="k", ax=ax, label="Total current")

        if len(lower_positive_currents) > 1:
            color = next(prop_cycle)["color"]
            sum_lower_positive.plot(
                ls="-", c=color, ax=ax, label="Lower sheaths positive total"
            )
            sum_lower_negative.plot(
                ls="--", c=color, ax=ax, label="Lower sheaths negative total"
            )

        if len(upper_positive_currents) > 1:
            color = next(prop_cycle)["color"]
            sum_upper_positive.plot(
                ls="-", c=color, ax=ax, label="Upper sheaths positive total"
            )
            sum_upper_negative.plot(
                ls="--", c=color, ax=ax, label="Upper sheaths negative total"
            )

        for region_name, positive_current, negative_current in zip(
            lower_positive_currents,
            lower_positive_currents.values(),
            lower_negative_currents.values(),
        ):
            color = next(prop_cycle)["color"]
            positive_current.plot(
                ls=":", c=color, ax=ax, label=f"{region_name} lower positive"
            )
            negative_current.plot(
                ls="-.", c=color, ax=ax, label=f"{region_name} lower negative"
            )

        for region_name, positive_current, negative_current in zip(
            upper_positive_currents,
            upper_positive_currents.values(),
            upper_negative_currents.values(),
        ):
            color = next(prop_cycle)["color"]
            positive_current.plot(
                ls=":", c=color, ax=ax, label=f"{region_name} upper positive"
            )
            negative_current.plot(
                ls="-.", c=color, ax=ax, label=f"{region_name} upper negative"
            )

        ax.set_ylabel("Parallel current / A")

        l = ax.legend()
        l.set_draggable(True)

        return ax

    elif output.lower() == "animate2d":
        plot_list = []
        for region_name, current in lower_currents.items():
            current.attrs[
                "long_name"
            ] = f"{current.attrs['long_name']} at {region_name} lower"
            plot_list.append(current)
        for region_name, current in upper_currents.items():
            current.attrs[
                "long_name"
            ] = f"{current.attrs['long_name']} at {region_name} upper"
            plot_list.append(current)
        return ds.bout.animate_list(plot_list, animate_over="time", **kwargs)
    else:
        raise ValueError(f"Unrecognised argument output='{output}'")
