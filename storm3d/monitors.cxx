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

void STORM::printMinMaxMean() {
  output.write("\nmin(phi) = %e, max(phi) = %e, mean(phi) = %e\n",
    min(phi,true,"RGN_NOBNDRY"),max(phi,true,"RGN_NOBNDRY"),mean(phi,true,"RGN_NOBNDRY"));

  output.write("min(n) = %e, max(n) = %e, mean(n) = %e, ",
    min(n,true,"RGN_NOBNDRY"),max(n,true,"RGN_NOBNDRY"),mean(n,true,"RGN_NOBNDRY"));
  output.write("min(ddt(logn)) = %e, max(ddt(logn)) = %e, mean(ddt(logn)) = %e\n",
    min(ddt(logn),true,"RGN_NOBNDRY"),max(ddt(logn),true,"RGN_NOBNDRY"),mean(ddt(logn),true,"RGN_NOBNDRY"));

  if (not isothermal) {
    output.write("min(T) = %e, max(T) = %e, mean(T) = %e",
      min(T,true,"RGN_NOBNDRY"),max(T,true,"RGN_NOBNDRY"),mean(T,true,"RGN_NOBNDRY"));
    output.write(", min(ddt(logp)) = %e, max(ddt(logp)) = %e, mean(ddt(logp)) = %e\n",
      min(ddt(logp),true,"RGN_NOBNDRY"),max(ddt(logp),true,"RGN_NOBNDRY"),mean(ddt(logp),true,"RGN_NOBNDRY"));
  }

  output.write("min(U) = %e, max(U) = %e, mean(U) = %e, ",
    min(U_aligned,true,"RGN_NOBNDRY"),max(U_aligned,true,"RGN_NOBNDRY"),mean(U_aligned,true,"RGN_NOBNDRY"));
  output.write("min(ddt(chiU)) = %e, max(ddt(chiU)) = %e, mean(ddt(chiU)) = %e\n",
    min(ddt(chiU),true,"RGN_NOBNDRY"),max(ddt(chiU),true,"RGN_NOBNDRY"),mean(ddt(chiU),true,"RGN_NOBNDRY"));

  output.write("min(V) = %e, max(V) = %e, mean(V) = %e",
    min(V_aligned,true,"RGN_NOBNDRY"),max(V_aligned,true,"RGN_NOBNDRY"),mean(V_aligned,true,"RGN_NOBNDRY"));
  output.write(", min(ddt(chiV)) = %e, max(ddt(chiV)) = %e, mean(ddt(chiV)) = %e\n",
    min(ddt(chiV),true,"RGN_NOBNDRY"),max(ddt(chiV),true,"RGN_NOBNDRY"),mean(ddt(chiV),true,"RGN_NOBNDRY"));

  output.write("min(vort) = %e, max(vort) = %e, mean(vort) = %e, ",
    min(vort,true,"RGN_NOBNDRY"),max(vort,true,"RGN_NOBNDRY"),mean(vort,true,"RGN_NOBNDRY"));
  output.write("min(ddt(vort)) = %e, max(ddt(vort)) = %e, mean(ddt(vort)) = %e\n",
    min(ddt(vort),true,"RGN_NOBNDRY"),max(ddt(vort),true,"RGN_NOBNDRY"),mean(ddt(vort),true,"RGN_NOBNDRY"));

  if (electromagnetic) {
    output.write("min(psi) = %e, max(psi) = %e, mean(psi) = %e\n",
      min(psi,true,"RGN_NOBNDRY"),max(psi,true,"RGN_NOBNDRY"),mean(psi,true,"RGN_NOBNDRY"));
  }

  output.write("\n");
}

int STORM::outputMonitor(BoutReal UNUSED(simtime), int UNUSED(iteration), int UNUSED(nout)) {
  if (hydrodynamic) {
    phisolver_1d();
    phi.applyBoundary();
  }

  return 0;
}
