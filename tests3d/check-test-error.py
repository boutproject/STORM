#!/usr/bin/env python3

from xbout import open_boutdataset

# x-boundary cells are not set or used for these variables, so don't check them
noXBoundaryVariables = ['qpar', 'chiU', 'chiV']
# y-boundary cells are not set or used for these variables, so don't check them
noYBoundaryVariables = ['n', 'T', 'uE2', 'chiU', 'chiV']

result = open_boutdataset("data/BOUT.dmp.*.nc", keep_xboundaries=True, keep_yboundaries=True, info=False).isel(t=-1)
expected = open_boutdataset("data/expectedResults/BOUT.dmp.nc", keep_xboundaries=True, keep_yboundaries=True, info=False).isel(t=-1, x=slice(1, -1))

to_plot = []

for v in noXBoundaryVariables:
    if v not in result:
        continue
    r = result[v].isel(x=slice(1, -1))
    e = expected[v].isel(x=slice(1, -1))
    to_plot.append(r)
    to_plot.append(e)
    to_plot.append(r - e)

for v in noYBoundaryVariables:
    if v not in result:
        continue
    r = result[v].isel(y=slice(2, -2))
    e = expected[v].isel(y=slice(2, -2))
    to_plot.append(r)
    to_plot.append(e)
    to_plot.append(r - e)

result.bout.animate_list(to_plot, animate_over="z", aspect="auto", ncols=3, show=True)
