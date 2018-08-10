#!/usr/bin/env python3

# Copyright L. Easy, F. Militello, T. Nicholas, J. Omotani, F. Riva, N.
# Walkden, UKAEA, 2017, 2018
# email: fulvio.militello@ukaea.uk
#
# This file is part of the STORM module of BOUT++.
#
# STORM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STORM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with STORM.  If not, see <https://www.gnu.org/licenses/>.

from boutdata.data import BoutOutputs
from boututils.showdata import showdata
import numpy
from sys import exit, argv

# get object that provides access to output from BOUT++
outputs = BoutOutputs(path='.', xguards=False, yguards=False)

# choose variables to plot
varlist = ['n', 'T', 'vort', 'phi', 'U', 'V']

# get y-index (giving position in direction parallel to B)
if len(argv)>1:
    # get from command line argument
    yindex = int(argv[1])
else:
    # choose midpoint of grid as default
    yindex = outputs['n'].shape[2]//2

# make a numpy index slice
# indices are {t,x,y,z} so we will make an animation of values on the x-z plane
# at y=yindex
plot_indices = numpy.index_exp[:,:,yindex,:]

# show animation
showdata([outputs[var][plot_indices] for var in varlist], titles=varlist)

exit(0)
