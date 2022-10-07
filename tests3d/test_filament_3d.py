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
import os
import time
from datetime import timedelta
from collections import namedtuple

Test = namedtuple('test', ['name', 'runtime'])
def Duration(minutes, seconds):
    return timedelta(minutes=minutes, seconds=seconds)

tests = [Test("test-3d", Duration(3, 23))]

testcommand = "./runtest.py"
retestOption = " --retest"
time_tolerance = timedelta(seconds=4.)

def test_filament_3d(numProcs,retest=False):
    global tests,testcommand,retestOption

    start_time = time.monotonic()

    numTests = 0
    numFailures = 0
    failedTests = []
    currentDir = os.getcwd()
    testcommand = testcommand + " " + str(numProcs)
    if retest:
        print("Rechecking results - NOT running test examples")
        testcommand = testcommand+retestOption

    warnings = []
    for test in tests:
        testdir = test.name
        os.chdir(currentDir)
        os.chdir(testdir)

        test_start = time.monotonic()
        s,out = shell(testcommand,pipe=False)
        test_time = timedelta(seconds=time.monotonic() - test_start)

        numTests = numTests+1
        if s is not 0:
            numFailures = numFailures+1
            failedTests.append(testdir)

        this_warning = None
        if test_time - test.runtime > time_tolerance:
            this_warning = testdir + ' took '+str(test_time)+'. This is longer than the expected '+str(test.runtime)+'.'
        elif test.runtime - test_time  > time_tolerance:
            this_warning = testdir + ' took '+str(test_time)+'. This is faster than the expected '+str(test.runtime)+'.'
        if this_warning is not None:
            print(this_warning, flush=True)
            warnings.append(this_warning)

        print("", flush=True)

    if numFailures is 0:
        print("All "+str(numTests)+" tests passed")
    else:
        print(str(numFailures)+"/"+str(numTests)+" failed.")
        print("Failed tests:")
        for name in failedTests:
            print("    "+name)

    if warnings:
        for warning in warnings:
            print(warning)
        print('Expected times are from running on 1 node (48 cores) on the A3 (SKL) partition of Marconi, with optimized configuration of BOUT++. If a test is slower for the same case, check for performance regression. If it is faster, the expected time may need updating to account for improved performance.')

    end_time = time.monotonic()
    print ("Tests took "+str(timedelta(seconds=end_time-start_time)))

if __name__=="__main__":

    import argparse
    from sys import exit

    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("np",type=int)
    parser.add_argument("--retest",action="store_true",default=False)
    args = parser.parse_args()

    test_filament_3d(args.np,args.retest)

    exit(0)
