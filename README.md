STORM
-----

STORM plasma simulation model, built using BOUT++ framework.

This is the public release version of STORM. The model assumes cold ions, but
includes electron temperature evolution. The current development version, used
at UKAEA/CCFE may have significant features added that are not yet included in
this public release. If you are interested in using features that have not been
released yet, please contact us to set up a collaboration.

This version of STORM has been tested with BOUT++ v4.4.2.

Table of Contents
-----------------

* [Options](/doc/options.md)

License
-------

If you use STORM, please cite the relevant papers (see
[References](#references)) in any publication. We also request that you contact
us (email:fulvio.militello@ukaea.uk) if you want to use or modify STORM.

Full text of the license is in the file COPYING.

  Copyright L. Easy, F. Militello, T. Nicholas, J. Omotani, F. Riva, N.
  Walkden, UKAEA, 2017, 2018
  email: fulvio.militello@ukaea.uk

  This file is part of the STORM module of BOUT++.

  STORM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  STORM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with STORM.  If not, see <https://www.gnu.org/licenses/>.

Quickstart
----------

To compile STORM:
1. [download and compile
   BOUT++](https://bout-dev.readthedocs.io/en/stable/user_docs/quickstart.html#building-bout).
2. Configure STORM to link to your copy of BOUT++. In this directory, run
   ```bash
   cmake . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=../build_bout
   ```
   This command assumes that you followed the BOUT++ instructions to build
   BOUT++ in a directory called `build_bout`, which is located one level up
   from this STORM directory. Otherwise modify the path passed to
   `-DCMAKE_PREFIX_PATH`.
3. Compile STORM
   ```bash
   cmake --build build
   ```

Test
----

To test STORM is working correctly, you can use an example 3d filament simulation.
1. Install the `boutdata` Python package (used to load output from the simulation)
   ```bash
   pip install --user boutdata
   ```
   You can also install `boutdata` using `conda`.
2. Run the test. For example on a machine with 8 cores
   ```bash
   cd tests3d
   ./test_filament_3d.py 8
   ```
   The terminal output should report
   ```
   8/8 tests passed in test-3d
   ```
   The test should take ~1 hour in serial (on a 3Ghz Intel Core i7 CPU), but
   also should scale well at least up to 48 cores.

References
----------

[1] L. Easy, F. Militello, J. Omotani, B. Dudson, E. Havlíčková, P. Tamain, V.
Naulin, and A. H. Nielsen, Physics of Plasmas 21, 122515 (2014),
https://doi.org/10.1063/1.4904207

[2] L. Easy, F. Militello, J. Omotani, N.R. Walkden, and B. Dudson, Physics of
Plasmas 23, 012512 (2016), https://doi.org/10.1063/1.4940330

[3] N.R. Walkden, L. Easy, F. Militello and J.T. Omotani, Plasma Physics and
Controlled Fusion 58, 115010 (2016),
https://doi.org/10.1088/0741-3335/58/11/115010

[4] F. Militello, B. Dudson, L. Easy, A. Kirk and P. Naylor, Plasma Physics and
Controlled Fusion 59, 125013 (2017), https://doi.org/10.1088/1361-6587/aa9252

[5] D. Hoare, F. Militello, J.T. Omotani, F. Riva, S. Newton, T. Nicholas, D.
Ryan and N.R. Walkden, Plasma Physics and Controlled Fusion 61, 105013 (2019),
https://doi.org/10.1088/1361-6587/ab34f8
