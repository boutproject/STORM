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

#include <boutcomm.hxx>
#include "fast_output.hxx"
#include <../src/fileio/formatfactory.hxx>
#include <bout/solver.hxx>
#include <string>

FastOutput::FastOutput() {
  // Get options
  Options* global_options = Options::getRoot();
  Options* options = global_options->getSection("fast_output");
  std::string type;
  OPTION(options, type, "none");
  if (type=="monitor") {
    enable_monitor = true;
    // Calculate output frequency
    BoutReal output_timestep;
    global_options->get("timestep", output_timestep, 1.);
    int frequency_multiplier;
    OPTION(options, frequency_multiplier, 100); // multiple of the output frequency to call fast_output at
    setTimestep(output_timestep / double(frequency_multiplier));
  } else if (type == "timestep") {
    enable_timestep = true;
    // Check that solver:monitor_timestep is true, otherwise type=timestep does nothing
    bool monitor_timestep;
    OPTION(global_options, monitor_timestep, false);
    if (!monitor_timestep) {
      throw BoutException("fast_output:type=timestep but solver:monitor_timestep=false, so FastOutput would do nothing. You probably want to set solver:monitor_timestep=true.");
    }
  } else if (type == "none") {
    // don't enable anything
  } else {
    throw BoutException("Unrecognized option for FastOutput: type=%s", type.c_str());
  }
  if (enable_monitor || enable_timestep) {
    enabled = true;
    // Read more options
    std::string prefix;
    OPTION(options, prefix, "BOUT.fast");
    std::string dump_ext;
    OPTION(options, dump_ext, "nc");
    std::string datadir;
    OPTION(global_options, datadir, "data");
    bool append;
    OPTION(global_options, append, false);

    std::string filename = "%s/" + prefix + ".%s";

    // Initialize output_file
    output_file = Datafile(options);
    if (append) {
      output_file.opena(filename.c_str(), datadir.c_str(), dump_ext.c_str());
    } else {
      output_file.openw(filename.c_str(), datadir.c_str(), dump_ext.c_str());
    }

    // Add the time to the output
    int MYPE;
    MPI_Comm_rank(BoutComm::get(), &MYPE);
    if (MYPE == 0)
      output_file.add(current_time, "time", true);
  }
}

void FastOutput::add(const std::string name, Field3D &f, const int ix_global, const int iy_global, const int iz_global) {
  int ix = ix_global - mesh->OffsetX;
  int iy = iy_global - mesh->OffsetY;
  int iz = iz_global - mesh->OffsetZ;

  if (ix>=mesh->xstart && ix<=mesh->xend && iy>=mesh->ystart && iy<=mesh->yend && iz>=0 && iz<mesh->LocalNz) {
    // Store a reference to the field
    field3d_list.push_back(&f);

    // Store the location of the element to output
    Ind3D i(0, mesh->LocalNy, mesh->LocalNz);
    i = i.offset(ix, iy, iz);
    field3d_inds.push_back(i);

    // Extend output_vals
    output_vals.push_back(0.);

    // Add to the output_file
    output_file.add(output_vals.back(), name.c_str(), true);
    // Store the indices as attributes of the variable
    output_file.setAttribute(name, "ix", ix_global);
    output_file.setAttribute(name, "iy", iy_global);
    output_file.setAttribute(name, "iz", iz_global);
  }
}

int FastOutput::monitor_method(BoutReal simtime) {
  // Set time
  current_time = simtime;

  // Set values to output
  for (size_t i=0; i<output_vals.size(); i++) {
    output_vals[i] = field3d_list[i]->operator[](field3d_inds[i]);
  }

  output_file.write();

  return 0;
}
