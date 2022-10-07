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
#include <bout/constants.hxx>

int i, j, k ; 


// Set default options for BOUT++ library-applied boundary conditions.
// We only support one option, so this avoids mistakes in input files and also makes them
// shorter.
// Values set here can still be overridden in input files, in case this is useful.
// Read and assign values here so that the default is taken as an explicitly set value, so
// it won't be replaced by a different default later.
void STORM::setBoundaryConditionsOptions() {
  auto& opt = *globalOptions;

  // density
  if (!average_radial_boundaries_core_SOL) {
    opt["logn"]["bndry_xin"] = opt["logn"]["bndry_xin"].withDefault("neumann_o2");
    opt["logn"]["bndry_xout"] = opt["logn"]["bndry_xout"].withDefault("neumann_o2");
  } else {
    opt["logn"]["bndry_xin"] = opt["logn"]["bndry_xin"].withDefault("none");
    opt["logn"]["bndry_xout"] = opt["logn"]["bndry_xout"].withDefault("none");
  }
  // y-boundary guard cells only used from n_aligned
  opt["logn"]["bndry_ydown"] = opt["logn"]["bndry_ydown"].withDefault("none");
  opt["logn"]["bndry_yup"] = opt["logn"]["bndry_yup"].withDefault("none");

  // x-boundary guard cells of n_aligned not used
  opt["logn_aligned"]["bndry_xin"] = opt["logn_aligned"]["bndry_xin"].withDefault("none");
  opt["logn_aligned"]["bndry_xout"] = opt["logn_aligned"]["bndry_xout"].withDefault("none");
  if (symmetry_plane) {
    opt["logn_aligned"]["bndry_ydown"] = opt["logn_aligned"]["bndry_ydown"].withDefault("neumann_o2");
  } else {
    opt["logn_aligned"]["bndry_ydown"] = opt["logn_aligned"]["bndry_ydown"].withDefault("free_o3");
  }
  opt["logn_aligned"]["bndry_yup"] = opt["logn_aligned"]["bndry_yup"].withDefault("free_o3");

  // vorticity
  if (!average_radial_boundaries_core_SOL) {
    opt["vort"]["bndry_xin"] = opt["vort"]["bndry_xin"].withDefault("neumann_o2");
    opt["vort"]["bndry_xout"] = opt["vort"]["bndry_xout"].withDefault("neumann_o2");
  } else {
    opt["vort"]["bndry_xin"] = opt["vort"]["bndry_xin"].withDefault("none");
    opt["vort"]["bndry_xout"] = opt["vort"]["bndry_xout"].withDefault("none");
  }
  // y-boundary guard cells only used from vort_aligned
  opt["vort"]["bndry_ydown"] = opt["vort"]["bndry_ydown"].withDefault("none");
  opt["vort"]["bndry_yup"] = opt["vort"]["bndry_yup"].withDefault("none");

  // x-boundary guard cells of vort_aligned not used
  opt["vort_aligned"]["bndry_xin"] = opt["vort_aligned"]["bndry_xin"].withDefault("none");
  opt["vort_aligned"]["bndry_xout"] = opt["vort_aligned"]["bndry_xout"].withDefault("none");
  if (symmetry_plane) {
    opt["vort_aligned"]["bndry_ydown"] = opt["vort_aligned"]["bndry_ydown"].withDefault("neumann_o2");
  } else {
    opt["vort_aligned"]["bndry_ydown"] = opt["vort_aligned"]["bndry_ydown"].withDefault("free_o3");
  }
  opt["vort_aligned"]["bndry_yup"] = opt["vort_aligned"]["bndry_yup"].withDefault("free_o3");

  // electrostatic potential
  // x-boundary guard cells set by Laplacian solver
  opt["phi"]["bndry_xin"] = opt["phi"]["bndry_xin"].withDefault("none");
  opt["phi"]["bndry_xout"] = opt["phi"]["bndry_xout"].withDefault("none");
  // y-boundary guard cells only used from phi_aligned
  opt["phi"]["bndry_ydown"] = opt["phi"]["bndry_ydown"].withDefault("none");
  opt["phi"]["bndry_yup"] = opt["phi"]["bndry_yup"].withDefault("none");

  // x-boundary guard cells of phi_aligned not used
  opt["phi_aligned"]["bndry_xin"] = opt["phi_aligned"]["bndry_xin"].withDefault("none");
  opt["phi_aligned"]["bndry_xout"] = opt["phi_aligned"]["bndry_xout"].withDefault("none");
  if (symmetry_plane) {
    opt["phi_aligned"]["bndry_ydown"] = opt["phi_aligned"]["bndry_ydown"].withDefault("neumann_o2");
  } else {
    opt["phi_aligned"]["bndry_ydown"] = opt["phi_aligned"]["bndry_ydown"].withDefault("free_o3");
  }
  opt["phi_aligned"]["bndry_yup"] = opt["phi_aligned"]["bndry_yup"].withDefault("free_o3");

  if (run_1d) {
    // Not enough x-points to use free_o3 boundary conditions, so use neumann instead
    opt["phi_stag"]["bndry_xin"] = opt["phi_stag"]["bndry_xin"].withDefault("neumann_o2");
    opt["phi_stag"]["bndry_xout"] = opt["phi_stag"]["bndry_xout"].withDefault("neumann_o2");
  } else {
    opt["phi_stag"]["bndry_xin"] = opt["phi_stag"]["bndry_xin"].withDefault("free_o3");
    opt["phi_stag"]["bndry_xout"] = opt["phi_stag"]["bndry_xout"].withDefault("free_o3");
  }
  // y-boundary guard cells not used
  opt["phi_stag"]["bndry_ydown"] = opt["phi_stag"]["bndry_ydown"].withDefault("none");
  opt["phi_stag"]["bndry_yup"] = opt["phi_stag"]["bndry_yup"].withDefault("none");

  // electron temperature
  if (!isothermal) {
    if (!average_radial_boundaries_core_SOL) {
      opt["logp"]["bndry_xin"] = opt["logp"]["bndry_xin"].withDefault("neumann_o2");
      opt["logp"]["bndry_xout"] = opt["logp"]["bndry_xout"].withDefault("neumann_o2");
    } else {
      opt["logp"]["bndry_xin"] = opt["logp"]["bndry_xin"].withDefault("none");
      opt["logp"]["bndry_xout"] = opt["logp"]["bndry_xout"].withDefault("none");
    }
    // y-boundary guard cells only used from T_aligned
    opt["logp"]["bndry_ydown"] = opt["logp"]["bndry_ydown"].withDefault("none");
    opt["logp"]["bndry_yup"] = opt["logp"]["bndry_yup"].withDefault("none");

    // x-boundary guard cells of T_aligned not used
    opt["logT_aligned"]["bndry_xin"] = opt["logT_aligned"]["bndry_xin"].withDefault("none");
    opt["logT_aligned"]["bndry_xout"] = opt["logT_aligned"]["bndry_xout"].withDefault("none");
    if (symmetry_plane) {
      opt["logT_aligned"]["bndry_ydown"] = opt["logT_aligned"]["bndry_ydown"].withDefault("neumann_o2");
    } else {
      opt["logT_aligned"]["bndry_ydown"] = opt["logT_aligned"]["bndry_ydown"].withDefault("free_o3");
    }
    opt["logT_aligned"]["bndry_yup"] = opt["logT_aligned"]["bndry_yup"].withDefault("free_o3");
  }

  // ion velocity
  if (!average_radial_boundaries_core_SOL) {
    opt["U"]["bndry_xin"] = opt["U"]["bndry_xin"].withDefault("neumann_o2");
    opt["U"]["bndry_xout"] = opt["U"]["bndry_xout"].withDefault("neumann_o2");
  } else {
    opt["U"]["bndry_xin"] = opt["U"]["bndry_xin"].withDefault("none");
    opt["U"]["bndry_xout"] = opt["U"]["bndry_xout"].withDefault("none");
  }
  // y-boundary guard cells only used from U_aligned
  opt["U"]["bndry_ydown"] = opt["U"]["bndry_ydown"].withDefault("none");
  opt["U"]["bndry_yup"] = opt["U"]["bndry_yup"].withDefault("none");

  // x-boundary guard cells of U_aligned not used
  opt["U_aligned"]["bndry_xin"] = opt["U_aligned"]["bndry_xin"].withDefault("none");
  opt["U_aligned"]["bndry_xout"] = opt["U_aligned"]["bndry_xout"].withDefault("none");
  if (symmetry_plane) {
    opt["U_aligned"]["bndry_ydown"] = opt["U_aligned"]["bndry_ydown"].withDefault("dirichlet_o2");
  } else {
    // y-boundary guard cells set by sheath boundary conditions function
    opt["U_aligned"]["bndry_ydown"] = opt["U_aligned"]["bndry_ydown"].withDefault("none");
  }
  // y-boundary guard cells set by sheath boundary conditions function
  opt["U_aligned"]["bndry_yup"] = opt["U_aligned"]["bndry_yup"].withDefault("none");

  // electron velocity
  if (!average_radial_boundaries_core_SOL) {
    opt["V"]["bndry_xin"] = opt["V"]["bndry_xin"].withDefault("neumann_o2");
    opt["V"]["bndry_xout"] = opt["V"]["bndry_xout"].withDefault("neumann_o2");
  } else {
    opt["V"]["bndry_xin"] = opt["V"]["bndry_xin"].withDefault("none");
    opt["V"]["bndry_xout"] = opt["V"]["bndry_xout"].withDefault("none");
  }
  // y-boundary guard cells only used from V_aligned
  opt["V"]["bndry_ydown"] = opt["V"]["bndry_ydown"].withDefault("none");
  opt["V"]["bndry_yup"] = opt["V"]["bndry_yup"].withDefault("none");

  // x-boundary guard cells of V_aligned not used
  opt["V_aligned"]["bndry_xin"] = opt["V_aligned"]["bndry_xin"].withDefault("none");
  opt["V_aligned"]["bndry_xout"] = opt["V_aligned"]["bndry_xout"].withDefault("none");
  if (symmetry_plane) {
    opt["V_aligned"]["bndry_ydown"] = opt["V_aligned"]["bndry_ydown"].withDefault("dirichlet_o2");
  } else {
    // y-boundary guard cells set by sheath boundary conditions function
    opt["V_aligned"]["bndry_ydown"] = opt["V_aligned"]["bndry_ydown"].withDefault("none");
  }
  // y-boundary guard cells set by sheath boundary conditions function
  opt["V_aligned"]["bndry_yup"] = opt["V_aligned"]["bndry_yup"].withDefault("none");

  // chiU and chiV do not need boundary conditions - we never take derivatives
  // of either of them, only of psi, U and V
  opt["chiU"]["bndry_xin"] = opt["chiU"]["bndry_xin"].withDefault("none");
  opt["chiU"]["bndry_xout"] = opt["chiU"]["bndry_xout"].withDefault("none");
  opt["chiU"]["bndry_ydown"] = opt["chiU"]["bndry_ydown"].withDefault("none");
  opt["chiU"]["bndry_yup"] = opt["chiU"]["bndry_yup"].withDefault("none");
  opt["chiV"]["bndry_xin"] = opt["chiV"]["bndry_xin"].withDefault("none");
  opt["chiV"]["bndry_xout"] = opt["chiV"]["bndry_xout"].withDefault("none");
  opt["chiV"]["bndry_ydown"] = opt["chiV"]["bndry_ydown"].withDefault("none");
  opt["chiV"]["bndry_yup"] = opt["chiV"]["bndry_yup"].withDefault("none");

  // electron parallel heat flux
  // x-boundaries are not needed
  if (!isothermal) {
    opt["qpar_aligned"]["bndry_xin"] = opt["qpar_aligned"]["bndry_xin"].withDefault("none");
    opt["qpar_aligned"]["bndry_xout"] = opt["qpar_aligned"]["bndry_xout"].withDefault("none");
    if (symmetry_plane) {
      opt["qpar_aligned"]["bndry_ydown"] = opt["qpar_aligned"]["bndry_ydown"].withDefault("dirichlet_o2");
    } else {
      // y-boundary guard cells set by sheath boundary conditions function
      opt["qpar_aligned"]["bndry_ydown"] = opt["qpar_aligned"]["bndry_ydown"].withDefault("none");
    }
    // y-boundary guard cells set by sheath boundary conditions function
    opt["qpar_aligned"]["bndry_yup"] = opt["qpar_aligned"]["bndry_yup"].withDefault("none");
  }
  if (electromagnetic && !isothermal) {
    opt["qpar_centre"]["bndry_xin"] = opt["qpar_centre"]["bndry_xin"].withDefault("neumann_o2");
    opt["qpar_centre"]["bndry_xout"] = opt["qpar_centre"]["bndry_xout"].withDefault("neumann_o2");
    // y-boundaries not needed
    opt["qpar_centre"]["bndry_ydown"] = opt["qpar_centre"]["bndry_ydown"].withDefault("none");
    opt["qpar_centre"]["bndry_yup"] = opt["qpar_centre"]["bndry_yup"].withDefault("none");
  }

  // ExB speed squared
  if (boussinesq == 0) {
    if (run_1d) {
      // Not enough x-points to use free_o3 boundary conditions, so use neumann instead
      opt["uE2"]["bndry_xin"] = opt["uE2"]["bndry_xin"].withDefault("neumann_o2");
      opt["uE2"]["bndry_xout"] = opt["uE2"]["bndry_xout"].withDefault("neumann_o2");
    } else {
      opt["uE2"]["bndry_xin"] = opt["uE2"]["bndry_xin"].withDefault("free_o3");
      opt["uE2"]["bndry_xout"] = opt["uE2"]["bndry_xout"].withDefault("free_o3");
    }
    opt["uE2"]["bndry_ydown"] = opt["uE2"]["bndry_ydown"].withDefault("none");
    opt["uE2"]["bndry_yup"] = opt["uE2"]["bndry_yup"].withDefault("none");
  }

  // parallel component of magnetic vector potential
  if (electromagnetic) {
    // x-boundary guard cells set by Laplacian solver
    opt["psi"]["bndry_xin"] = opt["psi"]["bndry_xin"].withDefault("none");
    opt["psi"]["bndry_xout"] = opt["psi"]["bndry_xout"].withDefault("none");
    // y-boundary guard cells only used from psi_aligned
    opt["psi"]["bndry_ydown"] = opt["psi"]["bndry_ydown"].withDefault("none");
    opt["psi"]["bndry_yup"] = opt["psi"]["bndry_yup"].withDefault("none");

    // x-boundary guard cells of psi_aligned not used
    opt["psi_aligned"]["bndry_xin"] = opt["psi_aligned"]["bndry_xin"].withDefault("none");
    opt["psi_aligned"]["bndry_xout"] = opt["psi_aligned"]["bndry_xout"].withDefault("none");
    if (symmetry_plane) {
      opt["psi_aligned"]["bndry_ydown"] =
        opt["psi_aligned"]["bndry_ydown"].withDefault("dirichlet_o2");
    } else {
      if (!use_psi_boundary_solver) {
        opt["psi_aligned"]["bndry_ydown"] =
          opt["psi_aligned"]["bndry_ydown"].withDefault("free_o3");
      } else {
        // y-boundary guard cells set by sheath boundary conditions function
        opt["psi_aligned"]["bndry_ydown"] =
          opt["psi_aligned"]["bndry_ydown"].withDefault("none");
      }
    }
    if (!use_psi_boundary_solver) {
      opt["psi_aligned"]["bndry_yup"] =
        opt["psi_aligned"]["bndry_yup"].withDefault("free_o3");
    } else {
      // y-boundary guard cells set by sheath boundary conditions function
      opt["psi_aligned"]["bndry_yup"] =
        opt["psi_aligned"]["bndry_yup"].withDefault("none");
    }

    opt["psi_centre"]["bndry_xin"] =
      opt["psi_centre"]["bndry_xin"].withDefault("dirichlet_o2");
    opt["psi_centre"]["bndry_xout"] =
      opt["psi_centre"]["bndry_xout"].withDefault("dirichlet_o2");
    // y-boundary guard cells not needed
    opt["psi_centre"]["bndry_ydown"] =
      opt["psi_centre"]["bndry_ydown"].withDefault("none");
    opt["psi_centre"]["bndry_yup"] = opt["psi_centre"]["bndry_yup"].withDefault("none");
  }
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// Electron Velocity Sheath boundaries 
//////////////////////////////////////////////////////////////////////////////////////////////////

void STORM::Vsheath_yup_staggered(Field3D &var, const FieldPerp &phisheath, const BoutReal phi_wall, const FieldPerp &T_sheath, const BoutReal Vsheath_BC_prefactor){
  ASSERT1(var.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(var.getLocation() == CELL_YLOW);
  ASSERT1(phisheath.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(phisheath.getLocation() == CELL_YLOW);
  ASSERT1(T_sheath.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(T_sheath.getLocation() == CELL_YLOW);

  RangeIterator xrup = mesh->iterateBndryUpperY();
  BoutReal value ;

  j =  mesh->yend+1 ;

  for(xrup.first(); !xrup.isDone(); xrup.next()){
    i = xrup.ind ;
    for(int k=0; k <mesh->LocalNz; k++) {
      value = Vsheath_BC_prefactor*sqrt(T_sheath(i, k))*exp((phi_wall - phisheath(i, k))/T_sheath(i, k))  ;

      if (value > Vsheath_BC_prefactor*sqrt(T_sheath(i, k))) {
        // Limit the maximum electron flux into the sheath.
        // When phisheath < phi_wall the sheath attracts electrons, and the expression we
        // are using, derived for a strongly electron-repelling sheath, is not valid.
        value = Vsheath_BC_prefactor*sqrt(T_sheath(i, k));
      }

      var(i, j, k)   = value ;
      var(i, j+1, k) = 3.0*value - 3.0*var(i, j-1, k) + 1.0*var(i, j-2, k) ;
    }
  }
}

void STORM::Vsheath_ydown_staggered(Field3D &var, const FieldPerp &phisheath, const BoutReal phi_wall, const FieldPerp &T_sheath, const BoutReal Vsheath_BC_prefactor){
  ASSERT1(var.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(var.getLocation() == CELL_YLOW);
  ASSERT1(phisheath.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(phisheath.getLocation() == CELL_YLOW);
  ASSERT1(T_sheath.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(T_sheath.getLocation() == CELL_YLOW);

  RangeIterator xrup = mesh->iterateBndryLowerY();
  BoutReal value ;

  j =  mesh->ystart;
  for(xrup.first(); !xrup.isDone(); xrup.next()){
    i = xrup.ind ;
    for(int k=0; k <mesh->LocalNz; k++) {
      value =  -Vsheath_BC_prefactor*sqrt(T_sheath(i, k))*exp((phi_wall - phisheath(i, k))/T_sheath(i, k))  ;

      if (value < -Vsheath_BC_prefactor*sqrt(T_sheath(i, k))) {
        // Limit the maximum electron flux into the sheath.
        // When phisheath < phi_wall the sheath attracts electrons, and the expression we
        // are using, derived for a strongly electron-repelling sheath, is not valid.
        value = -Vsheath_BC_prefactor*sqrt(T_sheath(i, k));
      }

      var(i, j, k)   = value ;
      var(i, j-1, k) = 3.0*value - 3.0*var(i, j+1, k) + 1.0*var(i, j+2, k) ;
      var(i, j-2, k) = 3.0*var(i, j-1, k) - 3.0*value + 1.0*var(i, j+1, k) ;
    }
  }
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// Ion Velocity Sheath boundaries 
//////////////////////////////////////////////////////////////////////////////////////////////////

void STORM::Usheath_yup_staggered(Field3D &U, const FieldPerp &sqrtT, const BoutReal Usheath_BC_prefactor){
  ASSERT1(U.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(U.getLocation() == CELL_YLOW);
  ASSERT1(sqrtT.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(sqrtT.getLocation() == CELL_YLOW);
  ASSERT1(sqrtT.getIndex() == mesh->yend+1);

  RangeIterator xrup = mesh->iterateBndryUpperY();
  BoutReal Uvalue ;
  BoutReal normalised_sound_speed ;
  j =  mesh->yend+1;
  for(xrup.first(); !xrup.isDone(); xrup.next()){
    i = xrup.ind ;
    for(int k=0; k <mesh->LocalNz; k++) {
      normalised_sound_speed = Usheath_BC_prefactor*sqrtT(i, k) ;
      Uvalue = 3.0*U(i, j-1, k) - 3.0*U(i, j-2, k) + 1.0*U(i, j-3, k) ;
      if (Uvalue < normalised_sound_speed){
        U(i, j, k) = normalised_sound_speed ;
      }else{
        U(i, j, k) = Uvalue ;
      }
      U(i, j+1, k) = 3.0*U(i, j, k) - 3.0*U(i, j-1, k) + 1.0*U(i, j-2, k) ;
    }
  }
}

void STORM::Usheath_ydown_staggered(Field3D &U, const FieldPerp &sqrtT, const BoutReal Usheath_BC_prefactor){
  ASSERT1(U.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(U.getLocation() == CELL_YLOW);
  ASSERT1(sqrtT.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(sqrtT.getLocation() == CELL_YLOW);
  ASSERT1(sqrtT.getIndex() == mesh->ystart);

  RangeIterator xrup = mesh->iterateBndryLowerY();
  BoutReal Uvalue ;
  BoutReal normalised_sound_speed ;
  j =  mesh->ystart;
  for(xrup.first(); !xrup.isDone(); xrup.next()){
    i = xrup.ind ;
    for(int k=0; k <mesh->LocalNz; k++) {
      normalised_sound_speed = Usheath_BC_prefactor*sqrtT(i, k) ;
      Uvalue = 3.0*U(i, j+1, k) - 3.0*U(i, j+2, k) + 1.0*U(i, j+3, k) ;
      if (Uvalue > -normalised_sound_speed){
        U(i, j, k) = -normalised_sound_speed ;
      }else{
        U(i, j, k) = Uvalue ;
      }
      U(i, j-1, k) = 3.0*U(i, j, k) - 3.0*U(i, j+1, k) + 1.0*U(i, j+2, k) ;
      U(i, j-2, k) = 3.0*U(i, j-1, k) - 3.0*U(i, j, k) + 1.0*U(i, j+1, k) ;
    }
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////
// Parallel Heat Flux
//////////////////////////////////////////////////////////////////////////////////////////////////

void STORM::qsheath_yup_staggered(Field3D &var, const FieldPerp &T_sheath, const FieldPerp &n_sheath, const FieldPerp &V_sheath, const BoutReal mu){
  ASSERT1(var.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(var.getLocation() == CELL_YLOW);
  ASSERT1(T_sheath.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(T_sheath.getLocation() == CELL_YLOW);
  ASSERT1(T_sheath.getIndex() == mesh->yend+1);
  ASSERT1(n_sheath.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(n_sheath.getLocation() == CELL_YLOW);
  ASSERT1(n_sheath.getIndex() == mesh->yend+1);
  ASSERT1(V_sheath.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(V_sheath.getLocation() == CELL_YLOW);
  ASSERT1(V_sheath.getIndex() == mesh->yend+1);
  
  RangeIterator xrup = mesh->iterateBndryUpperY();
  BoutReal Vf = 0.5*log(TWOPI/mu) ;
  j =  mesh->yend+1 ;  
  
  for(xrup.first(); !xrup.isDone(); xrup.next()){
    i = xrup.ind ;
    for(int k=0; k <mesh->LocalNz; k++) {
      // q//e = n*T*V*( gamma - 5/2 - 1/(2*mu)*exp(-2*(Vfl + phi/T)) )
      var(i, j, k) = ( (fabs(Vf) - 0.5)*T_sheath(i, j, k) - (0.5/mu)*SQ(V_sheath(i, j, k)) )*n_sheath(i, j, k)*V_sheath(i, j, k);
      var(i, j+1, k) = 3.0*var(i, j, k) - 3.0*var(i, j-1, k) + 1.0*var(i, j-2, k) ;
    }
  }
}

void STORM::qsheath_ydown_staggered(Field3D &var, const FieldPerp &T_sheath, const FieldPerp &n_sheath, const FieldPerp &V_sheath, const BoutReal mu){
  ASSERT1(var.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(var.getLocation() == CELL_YLOW);
  ASSERT1(T_sheath.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(T_sheath.getLocation() == CELL_YLOW);
  ASSERT1(T_sheath.getIndex() == mesh->ystart);
  ASSERT1(n_sheath.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(n_sheath.getLocation() == CELL_YLOW);
  ASSERT1(n_sheath.getIndex() == mesh->ystart);
  ASSERT1(V_sheath.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(V_sheath.getLocation() == CELL_YLOW);
  ASSERT1(V_sheath.getIndex() == mesh->ystart);

  RangeIterator xrup = mesh->iterateBndryLowerY();
  BoutReal Vf = 0.5*log(TWOPI/mu) ;
  j =  mesh->ystart ;
  
  for(xrup.first(); !xrup.isDone(); xrup.next()){
    i = xrup.ind ;
    for(int k=0; k <mesh->LocalNz; k++) {
      // q//e = n*T*V*( gamma - 5/2 - 1/(2*mu)*exp(-2*(Vfl + phi/T)) )
      var(i, j, k) = ( (fabs(Vf) - 0.5)*T_sheath(i, j, k) - (0.5/mu)*SQ(V_sheath(i, j, k)) )*n_sheath(i, j, k)*V_sheath(i, j, k); 
      var(i, j-1, k) = 3.0*var(i, j, k) - 3.0*var(i, j+1, k) + 1.0*var(i, j+2, k) ;
      var(i, j-2, k) = 3.0*var(i, j-1, k) - 3.0*var(i, j, k) + 1.0*var(i, j+1, k) ;      
    }
  }
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// Phi perpendicular boundary condition
//////////////////////////////////////////////////////////////////////////////////////////////////

void STORM::set_xguards(Field3D &f, const BoutReal *g_inner, const BoutReal *g_outer) {
  if (mesh->firstX()) {
    for (i=mesh->xstart-1; i>=0; i--) {
      for (j=0; j < mesh->LocalNy ; j++) {
        int jglobal = mesh->getGlobalYIndex(j);
        for (k=0; k < mesh->LocalNz ; k++) {
          f(i, j, k)   = g_inner[jglobal] ;
        }
      }
    }
  }
  if (mesh->lastX()) {
    for (i=mesh->xend+1; i<mesh->LocalNx; i++) {
      for (j=0; j < mesh->LocalNy ; j++) {
        int jglobal = mesh->getGlobalYIndex(j);
        for (k=0; k < mesh->LocalNz ; k++) {
          f(i, j, k)   = g_outer[jglobal] ;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Extrapolation to the targets
//////////////////////////////////////////////////////////////////////////////////////////////////

FieldPerp STORM::extrap_sheath_upper(const Field3D &var){
  ASSERT1(var.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(var.getLocation() == CELL_CENTRE);

  //3rd order extrapolation of a field to find the sheath value
  FieldPerp result;
  result = 0.375*sliceXZ(var, mesh->yend + 1) + 0.75*sliceXZ(var, mesh->yend) - 0.125*sliceXZ(var, mesh->yend - 1) ;
  result.setIndex(mesh->yend+1);
  result.setLocation(CELL_YLOW);
  result.setDirectionY(YDirectionType::Aligned);
  return result;
}

FieldPerp STORM::extrap_sheath_lower(const Field3D &var){
  ASSERT1(var.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(var.getLocation() == CELL_CENTRE);

  //3rd order extrapolation of a field to find the sheath value
  FieldPerp result;
  result = 0.375*sliceXZ(var, mesh->ystart - 1) + 0.75*sliceXZ(var, mesh->ystart) - 0.125*sliceXZ(var, mesh->ystart + 1) ;
  result.setIndex(mesh->ystart);
  result.setLocation(CELL_YLOW);
  result.setDirectionY(YDirectionType::Aligned);
  return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Time Evolving Radial Boundary Conditions
//////////////////////////////////////////////////////////////////////////////////////////////////

//Routine used to initialise the variables used for the averaged Neumann boundary conditions
void STORM::phi_bc_initialise(bool restarting){
  bool append;
  OPTION(globalOptions, append, false);
  if(restarting && append){
    first_step = false;
  }else{
    first_step = true;
  }

  // If monitor_timestep == false, the evolving boundary conditions are not working
  Options* solver_options = globalOptions->getSection("solver");
  bool monitor_timestep;
  solver_options->get("monitor_timestep", monitor_timestep, false);
  if(!monitor_timestep) throw BoutException("Evolving boundary phi: monitor_timstep = false.");

  restart.add(phi_bc,"phi_bc",0);
  restart.add(time_last_SOL,"time_last_SOL",0);
  restart.add(time_last_PF,"time_last_PF",0);
  restart.add(time_last_core,"time_last_core",0);
  restart.add(cstep_SOL,"cstep_SOL",0);
  restart.add(cstep_PF,"cstep_PF",0);
  restart.add(cstep_core,"cstep_core",0);
  if(!restarting){
    //Start the simulation from scratch
    phi_bc = 0.;
    time_last_SOL = 0.;
    time_last_PF = 0.;
    time_last_core = 0.;
    cstep_SOL = 0.;
    cstep_PF = 0.;
    cstep_core = 0.;
    for(int i = 0; i<mesh->LocalNy; ++i){
      phi_bc(mesh->xstart,i) = phi(mesh->xstart,i,0);
      phi_bc(mesh->xend,i)   = phi(mesh->xend,i,0);
    }
  }
}

//Compute the average of phi(xleft) and phi(xright).
//Update of phi_bc if time - time_last > time_update
void STORM::compute_bndry_phi(const BoutReal time) {
  if(first_step){
    first_step = false;
  }else{
    //Computing the time- and toroidal-average of phi(mesh->xstart,iy) and phi(Nx-mesh->xstart,iy)
    cstep_SOL  += 1.;
    cstep_PF   += 1.;
    cstep_core += 1.;
    if(mesh->firstX()){
      int ix = mesh->xstart;
      for(int iy=0; iy < mesh->LocalNy ; iy++){
        for(int iz=0; iz < mesh->LocalNz; iz++){
          phi_bc(ix-1,iy) += phi(ix,iy,iz)/((double)(mesh->LocalNz));
        }
      }
    }
    if(mesh->lastX()){
      int ix = mesh->xend;
      for(int iy=0; iy < mesh->LocalNy ; iy++){
        for(int iz=0; iz < mesh->LocalNz; iz++){
          phi_bc(ix+1,iy) += phi(ix,iy,iz)/((double)(mesh->LocalNz));
        }
      }
    }
  }

  //Updating the boundary value if time - time_last > time_update
  if(time - time_last_core > time_update_core){
    time_last_core = time;
    if(mesh->firstX() &&
           ((mesh->getGlobalYIndexNoBoundaries(mesh->ystart) > jyseps1_1 && mesh->getGlobalYIndexNoBoundaries(mesh->yend) <= jyseps2_1) ||
            (mesh->getGlobalYIndexNoBoundaries(mesh->ystart) > jyseps1_2 && mesh->getGlobalYIndexNoBoundaries(mesh->yend) <= jyseps2_2) ||
            (jyseps2_1 == jyseps1_2 && mesh->getGlobalYIndexNoBoundaries(mesh->ystart) > jyseps1_1 && mesh->getGlobalYIndexNoBoundaries(mesh->yend) <= jyseps2_2))) {
      for(int iy=0; iy < mesh->LocalNy ; iy++){
        phi_bc(mesh->xstart,iy) = phi_bc(mesh->xstart-1,iy)/cstep_core;
        phi_bc(mesh->xstart-1,iy) = 0.;
      }
    }
    cstep_core = 0.;
  }

  if(time - time_last_PF > time_update_PF){
    time_last_PF = time;
    if(mesh->firstX() &&
           ((mesh->getGlobalYIndexNoBoundaries(mesh->yend) <= jyseps1_1) ||
            (mesh->getGlobalYIndexNoBoundaries(mesh->ystart) > jyseps2_1 && mesh->getGlobalYIndexNoBoundaries(mesh->yend) <= jyseps1_2) ||
            (mesh->getGlobalYIndexNoBoundaries(mesh->ystart) > jyseps2_2))) {
      for(int iy=0; iy < mesh->LocalNy ; iy++){
        phi_bc(mesh->xstart,iy) = phi_bc(mesh->xstart-1,iy)/cstep_PF;
        phi_bc(mesh->xstart-1,iy) = 0.;
      }
    }
    cstep_PF = 0.;
  }

  if(time - time_last_SOL > time_update_SOL){
    time_last_SOL = time;
    if(mesh->lastX()){
      for(int iy=0; iy < mesh->LocalNy ; iy++){
        phi_bc(mesh->xend,iy) = phi_bc(mesh->xend+1,iy)/cstep_SOL;
        phi_bc(mesh->xend+1,iy) = 0.;
      }
    }
    cstep_SOL = 0.;
  }
}

//Apply the values in phi_bc to phi
void STORM::apply_bndry_phi(){
  if(mesh->firstX()){
    int ixs  = mesh->xstart;
    int ixsg = ixs - 1;
    for(int ix=ixsg; ix >= 0 ; ix--){
      for(int iy=0; iy < mesh->LocalNy ; iy++){
        for(int iz=0; iz < mesh->LocalNz ; iz++){
          phi(ix,iy,iz)  = phi_bc(ixs,iy);
        }
      }
    }
  }
  if (mesh->lastX()){
    int ixe = mesh->xend;
    int ixeg = ixe + 1;
    for(int ix=ixeg; ix < mesh->LocalNx ; ix++){
      for(int iy=0; iy < mesh->LocalNy ; iy++){
        for(int iz=0; iz < mesh->LocalNz ; iz++){
          phi(ix,iy,iz)  = phi_bc(ixe,iy);
        }
      }
    }
  }
}

// Apply the boundary conditions on psi
void STORM::applyPsiBoundaries() {

  // Apply x-boundary conditions, and possibly free_o3 y-boundary conditions
  psi.applyBoundary();
  psi_aligned.applyBoundary();

  if (!use_psi_boundary_solver) {
    // Just extrapolate psi at the boundaries, free_o3 y-boundary conditions should have
    // been set

    if (mesh->yend - mesh->ystart + 1 <= mesh->ystart) {
      // Number of y-grid points equals number of y-guard cells, so boundary condition
      // changed a point that needs to be communicated
      mesh->communicate(psi_aligned);
    }
    return;
  }

  if (mesh->hasBndryLowerY() && !symmetry_plane) {
    ASSERT1(U_sheath_lower.getDirectionY() == YDirectionType::Aligned);
    ASSERT1(U_sheath_lower.getLocation() == CELL_YLOW);
    ASSERT1(U_sheath_lower.getIndex() == mesh->ystart);
    ASSERT1(V_sheath_lower.getDirectionY() == YDirectionType::Aligned);
    ASSERT1(V_sheath_lower.getLocation() == CELL_YLOW);
    ASSERT1(V_sheath_lower.getIndex() == mesh->ystart);
    ASSERT1(n_sheath_lower.getDirectionY() == YDirectionType::Aligned);
    ASSERT1(n_sheath_lower.getLocation() == CELL_YLOW);
    ASSERT1(n_sheath_lower.getIndex() == mesh->ystart);

    // Apply lower y-boundary condition at ystart
    int y = mesh->ystart;

    // Solve Delp2(psi) = -Jpar
    // Solve in x-z orthogonal coordinates and transform to field-aligned to apply
    // Note n_sheath_lower, U_sheath_lower and V_sheath_lower are in field-aligned
    // coordinates.
    psi_sheath_lower = psiBoundarySolver->solve(
          -fromFieldAligned(n_sheath_lower*(U_sheath_lower - V_sheath_lower), "RGN_NOX"),
          psi_sheath_lower);
    psi_aligned = toFieldAligned(psi_sheath_lower, "RGN_NOX");

    // Above only works in divertor configurations where whole x-z plane at ystart is a
    // y-boundary, so just assume that's true here and extrapolate into boundary cells
    // everywhere
    for (int x=mesh->xstart; x<=mesh->xend; x++) {
      for (int z=0; z<mesh->LocalNz; z++) {
        psi_aligned(x, y-1, z) = 3.*psi_aligned(x, y, z) - 3.*psi_aligned(x, y+1, z)
                                 + psi_aligned(x, y+2, z);
        psi_aligned(x, y-2, z) = 3.*psi_aligned(x, y-1, z) - 3.*psi_aligned(x, y, z)
                                 + psi_aligned(x, y+1, z);
      }
    }
  }

  if (mesh->hasBndryUpperY()) {
    ASSERT1(U_sheath_upper.getDirectionY() == YDirectionType::Aligned);
    ASSERT1(U_sheath_upper.getLocation() == CELL_YLOW);
    ASSERT1(U_sheath_upper.getIndex() == mesh->yend+1);
    ASSERT1(V_sheath_upper.getDirectionY() == YDirectionType::Aligned);
    ASSERT1(V_sheath_upper.getLocation() == CELL_YLOW);
    ASSERT1(V_sheath_upper.getIndex() == mesh->yend+1);
    ASSERT1(n_sheath_upper.getDirectionY() == YDirectionType::Aligned);
    ASSERT1(n_sheath_upper.getLocation() == CELL_YLOW);
    ASSERT1(n_sheath_upper.getIndex() == mesh->yend+1);

    // Apply upper y-boundary condition at yend+1
    int y = mesh->yend+1;

    // Solve Delp2(psi) = -Jpar
    // Solve in x-z orthogonal coordinates and transform to field-aligned to apply
    // Note n_sheath_upper, U_sheath_upper and V_sheath_upper are in field-aligned
    // coordinates.
    psi_sheath_upper = psiBoundarySolver->solve(
          -fromFieldAligned(n_sheath_upper*(U_sheath_upper - V_sheath_upper), "RGN_NOX"),
          psi_sheath_upper);
    psi_aligned = toFieldAligned(psi_sheath_upper, "RGN_NOX");

    // Above only works in divertor configurations where whole x-z plane at yend+1 is a
    // y-boundary, so just assume that's true here and extrapolate into boundary cells
    // everywhere
    for (int x=mesh->xstart; x<=mesh->xend; x++) {
      for (int z=0; z<mesh->LocalNz; z++) {
        psi_aligned(x, y+1, z) = 3.*psi_aligned(x, y, z) - 3.*psi_aligned(x, y-1, z)
                                 + psi_aligned(x, y-2, z);
      }
    }
  }

  if (mesh->yend - mesh->ystart + 1 <= mesh->ystart) {
    // Number of y-grid points equals number of y-guard cells, so boundary condition
    // changed a point that needs to be communicated
    mesh->communicate(psi_aligned);
  }
}

// Apply neumann boundary condition on the DC component
void STORM::average_Z_bndry(Field3D& f, const bool left, const bool right, const bool periodic){
  // if left == true: apply bc to the left boundary
  // if right == true: apply bc to the right boundary
  // if periodic == true: apply bc only if periodc
  int ys = mesh->ystart;
  int ye = mesh->yend;

  if (mesh->firstX() && left && (!periodic || (periodic && mesh->periodicY(mesh->xstart)) ) ) {
    int ix = mesh->xstart;
    int ixm = ix-1;
    for(int iy=ys; iy<=ye; ++iy){
      BoutReal f_av = 0.;
      for(int iz=0; iz<mesh->LocalNz; ++iz){
        f_av += f(ix,iy,iz);
      }
      f_av /= ((double)(mesh->LocalNz));
      for(int iz=0; iz<mesh->LocalNz; ++iz){
        f(ixm,iy,iz) = 2.*f_av - f(ix,iy,iz);
      }
    }

    for (int ii=ixm; ii>0; --ii) {
      for (int iy=ys; iy<=ye; ++iy) {
        for (int iz=0; iz<mesh->LocalNz; ++iz) {
          f(ii-1,iy,iz) = 3.*f(ii,iy,iz) - 3.*f(ii+1,iy,iz) + f(ii+2,iy,iz);
        }
      }
    }

  }
  if(mesh->lastX() && right && (!periodic || (periodic && mesh->periodicY(mesh->xend)) ) ) {
    int ix = mesh->xend;
    int ixp=ix+1;
    for(int iy=ys; iy<=ye; ++iy){
      BoutReal f_av = 0.;
      for(int iz=0; iz<mesh->LocalNz; ++iz){
        f_av += f(ix,iy,iz);
      }
      f_av /= ((double)(mesh->LocalNz));
      for(int iz=0; iz<mesh->LocalNz; ++iz){
        f(ixp,iy,iz) = 2.*f_av - f(ix,iy,iz);
      }
    }

    for (int ii=ixp; ii<mesh->LocalNx-1; ++ii) {
      for (int iy=ys; iy<=ye; ++iy) {
        for (int iz=0; iz<mesh->LocalNz; ++iz) {
          f(ii+1,iy,iz) = 3.*f(ii,iy,iz) - 3.*f(ii-1,iy,iz) + f(ii-2,iy,iz);
        }
      }
    }

  }
}

// Apply neumann boundary condition on the YZ-averaed component
void STORM::average_YZ_bndry(Field3D& f, const bool left, const bool right, const bool periodic){
  if (mesh->firstX() && left && (!periodic || (periodic && mesh->periodicY(mesh->xstart)) ) ) {
    int ix = mesh->xstart;
    BoutReal f_av = 0.;
    for(int iy=mesh->ystart; iy<=mesh->yend; ++iy)
      for(int iz=0; iz<mesh->LocalNz; ++iz)
        f_av += f(ix,iy,iz);
    f_av /= BoutReal( (mesh->yend - mesh->ystart + 1) * (mesh->LocalNz) );

    BoutReal buff_send = f_av, buff_recv = 0.;
    int nproc;
    MPI_Comm comm = mesh->getYcomm(mesh->xstart);
    MPI_Comm_size(comm,&nproc);
    MPI_Allreduce(&buff_send, &buff_recv, 1, MPI_DOUBLE, MPI_SUM, comm);
    f_av = buff_recv/BoutReal(nproc);

    for(int iy=mesh->ystart; iy<=mesh->yend; ++iy){
      for(int iz=0; iz<mesh->LocalNz; ++iz){
        f(ix-1,iy,iz) = 2.*f_av - f(ix,iy,iz);
      }
    }
    for (int ii=ix-1; ii>0; --ii) {
      for (int iy=mesh->ystart; iy<=mesh->yend; ++iy) {
        for (int iz=0; iz<mesh->LocalNz; ++iz) {
          f(ii-1,iy,iz) = 3.*f(ii,iy,iz) - 3.*f(ii+1,iy,iz) + f(ii+2,iy,iz);
        }
      }
    }
  }

  if (mesh->lastX() && right && (!periodic || (periodic && mesh->periodicY(mesh->xend)) ) ) {
    int ix = mesh->xend;
    BoutReal f_av = 0.;
    for(int iy=mesh->ystart; iy<=mesh->yend; ++iy)
      for(int iz=0; iz<mesh->LocalNz; ++iz)
        f_av += f(ix,iy,iz);
    f_av /= BoutReal( (mesh->yend - mesh->ystart + 1) * (mesh->LocalNz) );

    BoutReal buff_send = f_av, buff_recv = 0.;
    int nproc;
    MPI_Comm comm = mesh->getYcomm(mesh->xend);
    MPI_Comm_size(comm,&nproc);
    MPI_Allreduce(&buff_send, &buff_recv, 1, MPI_DOUBLE, MPI_SUM, comm);
    f_av = buff_recv/BoutReal(nproc);

    for(int iy=mesh->ystart; iy<=mesh->yend; ++iy){
      for(int iz=0; iz<mesh->LocalNz; ++iz){
        f(ix+1,iy,iz) = 2.*f_av - f(ix,iy,iz);
      }
    }
    for (int ii=ix+1; ii<mesh->LocalNx-1; ++ii) {
      for (int iy=mesh->ystart; iy<=mesh->yend; ++iy) {
        for (int iz=0; iz<mesh->LocalNz; ++iz) {
          f(ii+1,iy,iz) = 3.*f(ii,iy,iz) - 3.*f(ii-1,iy,iz) + f(ii-2,iy,iz);
        }
      }
    }
  }
}
