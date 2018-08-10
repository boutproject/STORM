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

//////////////////////////////////////////////////////////////////////////////////////////////////
// Electron Velocity Sheath boundaries 
//////////////////////////////////////////////////////////////////////////////////////////////////

void STORM::Vsheath_yup_staggered(Field3D &var, const FieldPerp &phisheath, const BoutReal phi_wall, const FieldPerp &T_sheath, const BoutReal Vsheath_BC_prefactor){
  RangeIterator xrup = mesh->iterateBndryUpperY();
  BoutReal value ;

  j =  mesh->yend+1 ;

  for(xrup.first(); !xrup.isDone(); xrup.next()){
    i = xrup.ind ;
    for(int k=0; k <mesh->LocalNz; k++) {
      value = Vsheath_BC_prefactor*sqrt(T_sheath(i, k))*exp((phi_wall - phisheath(i, k))/T_sheath(i, k))  ;
      var(i, j, k)   = value ;
      var(i, j+1, k) = 3.0*value - 3.0*var(i, j-1, k) + 1.0*var(i, j-2, k) ;
    }
  }
}

void STORM::Vsheath_ydown_staggered(Field3D &var, const FieldPerp &phisheath, const BoutReal phi_wall, const FieldPerp &T_sheath, const BoutReal Vsheath_BC_prefactor){
  RangeIterator xrup = mesh->iterateBndryLowerY();
  BoutReal value ;

  j =  mesh->ystart;
  for(xrup.first(); !xrup.isDone(); xrup.next()){
    i = xrup.ind ;
    for(int k=0; k <mesh->LocalNz; k++) {
      value =  -Vsheath_BC_prefactor*sqrt(T_sheath(i, k))*exp((phi_wall - phisheath(i, k))/T_sheath(i, k))  ;
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

void STORM::qsheath_yup_staggered(Field3D &var, const FieldPerp &T_sheath, const FieldPerp &n_sheath, const FieldPerp &U_sheath, const FieldPerp &V_sheath, const BoutReal mu){
  
  RangeIterator xrup = mesh->iterateBndryUpperY();
  BoutReal Vf = 0.5*log(TWOPI/mu) ;
  j =  mesh->yend+1 ;  
  
  for(xrup.first(); !xrup.isDone(); xrup.next()){
    i = xrup.ind ;
    for(int k=0; k <mesh->LocalNz; k++) {
      var(i, j, k) = (2. + fabs(Vf))*n_sheath(i, j, k)*T_sheath(i, j, k)*U_sheath(i, j, k) - 2.5*n_sheath(i, j, k)*T_sheath(i, j, k)*V_sheath(i, j, k) - (0.5/mu)*n_sheath(i, j, k)*pow(V_sheath(i, j, k), 3) ;
      var(i, j+1, k) = 3.0*var(i, j, k) - 3.0*var(i, j-1, k) + 1.0*var(i, j-2, k) ;
    }
  }
}

void STORM::qsheath_ydown_staggered(Field3D &var, const FieldPerp &T_sheath, const FieldPerp &n_sheath, const FieldPerp &U_sheath, const FieldPerp &V_sheath, const BoutReal mu){

  RangeIterator xrup = mesh->iterateBndryLowerY();
  BoutReal Vf = 0.5*log(TWOPI/mu) ;
  j =  mesh->ystart ;
  
  for(xrup.first(); !xrup.isDone(); xrup.next()){
    i = xrup.ind ;
    for(int k=0; k <mesh->LocalNz; k++) {
      var(i, j, k) = (2. + fabs(Vf))*n_sheath(i, j, k)*T_sheath(i, j, k)*U_sheath(i, j, k) - 2.5*n_sheath(i, j, k)*T_sheath(i, j, k)*V_sheath(i, j, k) - (0.5/mu)*n_sheath(i, j, k)*pow(V_sheath(i, j, k), 3) ; 
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
        int jglobal = mesh->YGLOBAL(j) + mesh->ystart;
        for (k=0; k < mesh->LocalNz ; k++) {
          f(i, j, k)   = g_inner[jglobal] ;
        }
      }
    }
  }
  if (mesh->lastX()) {
    for (i=mesh->xend+1; i<mesh->LocalNx; i++) {
      for (j=0; j < mesh->LocalNy ; j++) {
        int jglobal = mesh->YGLOBAL(j) + mesh->ystart;
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
  //3rd order extrapolation of a field to find the sheath value
  FieldPerp result;
  result = 0.375*sliceXZ(var, mesh->yend + 1) + 0.75*sliceXZ(var, mesh->yend) - 0.125*sliceXZ(var, mesh->yend - 1) ;
  return result;
}

FieldPerp STORM::extrap_sheath_lower(const Field3D &var){
//3rd order extrapolation of a field to find the sheath value
FieldPerp result;
result = 0.375*sliceXZ(var, mesh->ystart - 1) + 0.75*sliceXZ(var, mesh->ystart) - 0.125*sliceXZ(var, mesh->ystart + 1) ;
return result;
}

