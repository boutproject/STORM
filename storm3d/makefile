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

SOURCEH = storm.hxx
SOURCEC = storm.cxx boundaries.cxx initialise.cxx operators.cxx utilities.cxx
TARGET = storm

GIT_VERSION := $(shell git describe --abbrev=0 --dirty --always --tags)
CXXFLAGS += -DSTORM_VERSION=\"$(GIT_VERSION)\"

include $(BOUT_TOP)/make.config
