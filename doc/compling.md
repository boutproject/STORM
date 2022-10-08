# Compiling

**Contents**
* [Storm](#storm)
* [Storm2D](#storm2d)

Storm
-----

The recommended method to compile STORM is using CMake. You will first need
BOUT++ compiled with CMake (github.com/boutproject/BOUT-configs might help).

Then to make an optimised build (in the subdirectory `build`, using up to 16
process to compile), run for example the following commands from the top level
`STORM` directory
```
cmake . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/BOUT++/build/directory
cmake --build build -j 16
```

or to make a debug build
```
cmake . -B build-debug -DCMAKE_BUILD_TYPE=Debug -DCMAKE_PREFIX_PATH=/path/to/BOUT++/debug/build/directory
cmake --build build-debug -j 16
```

An autotools build is also possible, see `storm3d/make.config.example`.

Storm2D
-------

Storm2D does not have a CMake setup yet, so to build:
* go to the `storm2d` subdirectory
* copy `make.config.example` to `make.config`
* edit `make.config` with the path to your compiled version of BOUT++
* build with `make -j 4` (to build using up to 4 processes)
