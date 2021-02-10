/*
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
*/

#ifndef __FAST_OUTPUT_H__
#define __FAST_OUTPUT_H__

#include <bout/monitor.hxx>
#include <bout/physicsmodel.hxx>
#include <bout/region.hxx>
#include <datafile.hxx>
#include <globals.hxx>
#include <deque>
class Solver;

// Subclass of BOUT++ Monitor class so we can pass to Solver::
class FastOutput : public Monitor {
  public:
    FastOutput();
    FastOutput(FastOutput &f) = delete;
    FastOutput &operator=(FastOutput &f) = delete;

    /// Add a point to the output: the element of f with global indices
    /// {ix,iy,iz} will be written when monitor_method() is called
    void add(const std::string name, Field3D &f, const int ix, const int iy, const int iz);

    /// Writes added points to output_file
    int monitor_method(BoutReal simtime);

    bool enabled=false, enable_monitor=false, enable_timestep=false;

    /// provide Monitor::call method which is called by the Solver
    int call(Solver* UNUSED(solver), BoutReal time, int UNUSED(iter), int UNUSED(nout)) {
      return monitor_method(time);
    }
  private:
    Datafile output_file;
    std::deque<Field3D*> field3d_list;
    std::deque<Ind3D> field3d_inds;
    BoutReal current_time;
    std::deque<BoutReal> output_vals;
};

#endif //__FAST_OUTPUT_H__
