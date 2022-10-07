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

#include "storm.hxx"

int STORM::outputMonitor(BoutReal UNUSED(simtime), int UNUSED(iteration), int UNUSED(nout)) {
  if (hydrodynamic) {
    phisolver_1d();
    phi.applyBoundary();
  }
  if (monitor_minmaxmean) {
    output.write("\nmin(phi) = %e, max(phi) = %e, mean(phi) = %e\n",
      min(phi,true,"RGN_NOBNDRY"),max(phi,true,"RGN_NOBNDRY"),mean(phi,true,"RGN_NOBNDRY"));
    output.write("min(n) = %e, max(n) = %e, mean(n) = %e\n",
      min(n,true,"RGN_NOBNDRY"),max(n,true,"RGN_NOBNDRY"),mean(n,true,"RGN_NOBNDRY"));
    output.write("min(T) = %e, max(T) = %e, mean(T) = %e\n",
      min(T,true,"RGN_NOBNDRY"),max(T,true,"RGN_NOBNDRY"),mean(T,true,"RGN_NOBNDRY"));
    output.write("min(U) = %e, max(U) = %e, mean(U) = %e\n",
      min(U_aligned,true,"RGN_NOBNDRY"),max(U_aligned,true,"RGN_NOBNDRY"),mean(U_aligned,true,"RGN_NOBNDRY"));
    output.write("min(V) = %e, max(V) = %e, mean(V) = %e\n",
      min(V_aligned,true,"RGN_NOBNDRY"),max(V_aligned,true,"RGN_NOBNDRY"),mean(V_aligned,true,"RGN_NOBNDRY"));
    output.write("min(vort) = %e, max(vort) = %e, mean(vort) = %e\n",
      min(vort,true,"RGN_NOBNDRY"),max(vort,true,"RGN_NOBNDRY"),mean(vort,true,"RGN_NOBNDRY"));
  }
  return 0;
}
