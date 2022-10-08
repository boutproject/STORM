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

# Set user-specified values
ifneq ("$(wildcard make.config)","")
  include make.config
else
  $(error "The 'make.config' file does not exist. You must create it for STORM to compile. See the example in 'make.config.example'")
endif

ifeq ($(wildcard ../shared/BoutEquation/equation.hxx),)
  $(error "BoutEquation header missing - try running `git submodule update --init --recursive`")
endif
ifeq ($(wildcard ../shared/BoutFastOutput/fast_output.hxx),)
  $(error "FastOutput header missing - try running `git submodule update --init --recursive`")
endif

SOURCEC = storm.cxx boundaries.cxx initialise.cxx monitors.cxx operators.cxx utilities.cxx ../shared/BoutFastOutput/fast_output.cxx ../shared/BoutEquation/equation.cxx neutral-rates.cxx D-vpar.cxx neutral-model.cxx
TARGET = storm

save_git_version_status := $(shell ../build_tools/save_git_version.py)
ifeq ($(save_git_version_status), 11)
  $(error "gitpython not installed. You can install with 'conda install gitpython' or 'pip3 install --user gitpython'")
endif
ifeq ($(save_git_version_status), 12)
  $(error "save_git_version.py failed to find a delimiter that is not in the git diff")
endif
ifneq ($(save_git_version_status), 0)
  $(error "save_git_version.py failed with an unrecognised error $(save_git_version_status)")
endif

$(shell [ -e storm.o ] && rm storm.o)

include $(BOUT_TOP)/make.config
