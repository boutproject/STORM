#!/usr/bin/env python3

from boututils.showdata import showdata
from boutdata import collect
from numpy import *
from matplotlib import pyplot as plt
import sys
import os

# T = collect('T', tind = [5500,6000], yguards = False)
phi = collect('phi', yguards = False)
n = collect('n',  yguards = False)
vort = collect('vort', yguards = False)

x = slice(2,-2)
y = 0
z = slice(None)
t = slice(None)

n =       n[t,x,y,z]
phi =   phi[t,x,y,z]
# T =       T[t,x,y,z]
vort = vort[t,x,y,z]

showdata([[n], [phi], [vort]], titles = ['n', 'phi', 'vort'])
