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

from boututils.run_wrapper import shell
import numpy
import os
import glob
import time
from datetime import datetime, timedelta

tolerance = 1.e-11
abs_tolerance = tolerance

testname = os.path.basename(os.getcwd())
equilibDir = "data/equilibrium/"
expectedEquilibDir = "data/expected_equilib"
testFiles = ["n_eq.dat","phi_eq.dat","U_eq.dat","V_eq.dat"]
ylowFiles = ["U_eq.dat", "V_eq.dat"]
optionalTestFiles = ["T_eq.dat"]
equilibFiles = [os.path.join(equilibDir,f) for f in testFiles+["background.mat","equilibrium.nc"] ]
outputDir = "data/"
outputFiles = [x for f in ["BOUT.dmp.*","BOUT.restart.*","BOUT.log.*"] for x in glob.glob(os.path.join(outputDir,f)) ]
createbgscript = "../../storm3d/create_bg_1d"

def test(retestOutput=False,nproc=1):

    start_time = time.monotonic()

    print("***************************************")
    print(testname)
    print("testing at "+str(datetime.now()))
    print("***************************************")

    if not retestOutput:
        run_test(nproc=nproc)

    numFailures,numTests = check_test()

    end_time = time.monotonic()
    print ("Test took "+str(timedelta(seconds=end_time-start_time)))

    return numFailures,numTests

def run_test(nproc=1):

    print("Running simulation")

    # Delete all output files before starting
    for filename in equilibFiles + outputFiles:
        if os.path.isfile(filename):
            os.remove(filename)
    s,out = shell(createbgscript+" --nproc "+str(nproc)+" > createbg.log")
    if s is not 0:
        raise ValueError("Failed to run "+createbgscript)

def check_test():

    print("Checking output")
    numFailures = 0
    numTests = 0

    # Check if optional test files exist, e.g. T_eq.dat for non-isothermal simulations
    for filename in optionalTestFiles:
        if os.path.isfile(os.path.join(expectedEquilibDir,filename)):
            testFiles.append(filename)

    # Check each equilibrium output against its expected value
    for filename in testFiles:
        data = numpy.fromfile(os.path.join(equilibDir,filename))
        expectedData = numpy.fromfile(os.path.join(expectedEquilibDir,filename))
        diff,norm = testfield(data, expectedData, filename in ylowFiles)
        numTests = numTests+1
        if diff/norm>tolerance and diff>abs_tolerance:
            numFailures = numFailures+1
            print("FAILURE: Test "+str(numTests)+" ("+filename+") failed, with diff/norm="+str(diff/norm)+" and diff="+str(diff))
        else:
            print("Test "+str(numTests)+" ("+filename+") passed, with diff/norm="+str(diff/norm)+" and diff="+str(diff))

    print(str(numTests-numFailures)+"/"+str(numTests)+" tests passed in "+testname)

    return numFailures,numTests

def testfield(data, expectedData, ylow):
    # don't check y-boundary cells as these are not set for non-aligned fields
    if not ylow:
        indices = numpy.index_exp[2:-2]
    else:
        indices = numpy.index_exp[3:-2]

    norm = numpy.sqrt(expectedData[indices]**2).mean()
    if norm==0.:
        norm==1.
    diff = numpy.sqrt(((data[indices] - expectedData[indices])**2).max())
    return diff,norm

if __name__=="__main__":

    import argparse
    from sys import exit

    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("np",nargs='?',type=int,default=1)
    parser.add_argument("--retest",action="store_true",default=False)
    args = parser.parse_args()

    numFailures,numTests = test(retestOutput=args.retest,nproc=args.np)
    exit(numFailures)
