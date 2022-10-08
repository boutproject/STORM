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

// Utility methods for STORM

#include "storm.hxx"
#include <derivs.hxx>
#include <interpolation.hxx>
#include <bout/coordinates.hxx>

#include <iostream>
#include <sstream>
#include <netcdf>

void STORM::phisolver_1d(){
  // Assume adiabaticity to get parallel profile of phi
  // phi = mu/(mu+1.)*log(n) for

  const int nype = mesh->getNYPE();
  int my_ype = mesh->getYProcIndex();
  MPI_Comm ycomm = mesh->getYcomm(0);
  const int nperp = mesh->LocalNx*mesh->LocalNz;

  Coordinates* coords = phi.getCoordinates();

  // Note, Grad_par_phi_stag should have been already computed in STORM::rhs() before phisolver_1d() is called by the output monitor
  if (symmetry_plane) {
    phi(mesh->xstart,mesh->yend,0) = 0;
    for (int i=mesh->xstart; i<=mesh->xend; i++) {
      for (int j=mesh->yend; j>mesh->ystart-1; j--) {
        for (int k=0; k<mesh->LocalNz; k++) {
          phi(i,j-1,k) = phi(i,j,k) - Grad_par_phi_stag(i,j,k)*coords->dy(i,j-1)*sqrt(coords->g_22(i,j-1));
        }
      }
    }

    phi.applyBoundary();
    phi_aligned = toFieldAligned(phi);
    phi_aligned.applyBoundary();

    // Make phi at the wall be as close to sheath potential as possible...
    if (mesh->lastY()) {
      FieldPerp phi_upper = extrap_sheath_upper(phi_aligned);
      FieldPerp sheath_potential = phi_wall + log(Vsheath_BC_prefactor*sqrt(T_sheath_upper)/U_sheath_upper)*T_sheath_upper;
      FieldPerp temp = fromFieldAligned(-phi_upper + sheath_potential, "RGN_NOX");
      // Set location of temp to be compatible with phi
      temp.setLocation(CELL_CENTRE);
      for (int j=mesh->ystart-1; j<=mesh->yend; j++) {
        // Add to ystart-1 as well, because ystart-1 is the integral that will be communicated
        temp.setIndex(j);
        phi = temp + phi; // this looks weird, but 'temp + phi' returns a FieldPerp, which is then assigned to an x-z slice of phi
      }
    }

    // Gather the integrals across each processor to the other processors
    FieldPerp this_integral = sliceXZ(phi,mesh->ystart-1);
    std::vector<BoutReal> local_integrals;
    local_integrals.resize(nype*nperp);
    int status = MPI_Allgather(&(this_integral(0,0)), nperp, MPI_DOUBLE, local_integrals.data(), nperp, MPI_DOUBLE, ycomm);
    ASSERT0(status == MPI_SUCCESS);

    if (my_ype<nype-1) {
      // Add integrals from processors above this one
      FieldPerp temp = 0.;
      for (int proc=nype-1; proc>my_ype; proc--) {
        for (int i=mesh->xstart; i<=mesh->xend; i++) {
          for (int k=0; k<mesh->LocalNz; k++) {
            temp(i,k) += local_integrals[proc*nperp+i*mesh->LocalNz+k];
          }
        }
      }
      for (int j=mesh->ystart; j<=mesh->yend; j++) {
        temp.setIndex(j);
        phi = temp + phi; // this looks weird, but 'temp + phi' returns a FieldPerp, which is then assigned to an x-z slice of phi
      }
    }

  } else {
    // Integrate symmetrically, one pass up and one down

    // Calculate local parts of integrals for each processor
    // only have Grad_par_phi_stag between ystart and yend (inclusive) because we don't communicate,
    // but Grad_par_phi_stag is at CELL_YLOW, so have the value between ystart-1 and ystart cell centres.
    // So start phi_up integrations at ystart-1 and go to yend, but phi_down integrations at yend and go to ystart-1

    Field3D phi_up = 0.;
    for (int i=mesh->xstart; i<=mesh->xend; i++) {
      for (int j=mesh->ystart-1; j<mesh->yend; j++) {
        for (int k=0; k<mesh->LocalNz; k++) {
          phi_up(i,j+1,k) = phi_up(i,j,k) + Grad_par_phi_stag(i,j+1,k)*coords->dy(i,j)*sqrt(coords->g_22(i,j));
        }
      }
    }
    Field3D phi_down = 0.;
    for (int i=mesh->xstart; i<=mesh->xend; i++) {
      for (int j=mesh->yend; j>mesh->ystart-1; j--) {
        for (int k=0; k<mesh->LocalNz; k++) {
          phi_down(i,j-1,k) = phi_down(i,j,k) - Grad_par_phi_stag(i,j,k)*coords->dy(i,j-1)*sqrt(coords->g_22(i,j-1));
        }
      }
    }

    // Gather the integrals across each processor to the other processors
    FieldPerp this_integral_up = sliceXZ(phi_up,mesh->yend);
    FieldPerp this_integral_down = sliceXZ(phi_down,mesh->ystart-1);
    std::vector<BoutReal> local_integrals_up, local_integrals_down;
    local_integrals_up.resize(nype*nperp);
    local_integrals_down.resize(nype*nperp);
    int status = MPI_Allgather(&(this_integral_up(0,0)), nperp, MPI_DOUBLE, local_integrals_up.data(), nperp, MPI_DOUBLE, ycomm);
    ASSERT0(status == MPI_SUCCESS);
    status = MPI_Allgather(&(this_integral_down(0,0)), nperp, MPI_DOUBLE, local_integrals_down.data(), nperp, MPI_DOUBLE, ycomm);
    ASSERT0(status == MPI_SUCCESS);

    if (my_ype>0) {
      FieldPerp temp = 0.;
      // Add integrals from processors below this one
      for (int proc=0; proc<my_ype; proc++) {
        for (int i=mesh->xstart; i<=mesh->xend; i++) {
          for (int k=0; k<mesh->LocalNz; k++) {
            temp(i,k) += local_integrals_up[proc*nperp+i*mesh->LocalNz+k];
          }
        }
      }
      for (int j=mesh->ystart; j<=mesh->yend; j++) {
        temp.setIndex(j);
        phi_up = temp + phi_up; // this looks weird, but 'temp + phi' returns a FieldPerp, which is then assigned to an x-z slice of phi
      }
    }
    if (my_ype<nype-1) {
      FieldPerp temp = 0.;
      // Add integrals from processors above this one
      for (int proc=nype-1; proc>my_ype; proc--) {
        for (int i=mesh->xstart; i<=mesh->xend; i++) {
          for (int k=0; k<mesh->LocalNz; k++) {
            temp(i,k) += local_integrals_down[proc*nperp+i*mesh->LocalNz+k];
          }
        }
      }
      temp = fromFieldAligned(temp, "RGN_NOX");
      for (int j=mesh->ystart; j<=mesh->yend; j++) {
        temp.setIndex(j);
        phi_down = temp + phi_down; // this looks weird, but 'temp + phi' returns a FieldPerp, which is then assigned to an x-z slice of phi
      }
    }

    // take phi as average of integral up and integral down
    phi = 0.5*(phi_up + phi_down);
    phi.applyBoundary();
    phi_aligned = toFieldAligned(phi);
    phi_aligned.applyBoundary();

    // Make phi at the wall be as close to sheath potentials as possible...
    FieldPerp offset_lower, offset_upper;
    offset_lower.allocate();
    offset_upper.allocate();
    offset_lower.setLocation(CELL_YLOW);
    offset_upper.setLocation(CELL_YLOW);
    if (mesh->firstY()) {
      //... for lower boundary
      FieldPerp phi_lower = extrap_sheath_lower(phi_aligned);
      FieldPerp sheath_potential_lower = phi_wall + log(Vsheath_BC_prefactor*sqrt(T_sheath_lower)/(-U_sheath_lower))*T_sheath_lower;
      offset_lower = sheath_potential_lower - phi_lower;
    }
    if (mesh->lastY()) {
      //... for upper boundary
      FieldPerp phi_upper = extrap_sheath_upper(phi_aligned);
      FieldPerp sheath_potential_upper = phi_wall + log(Vsheath_BC_prefactor*sqrt(T_sheath_upper)/U_sheath_upper)*T_sheath_upper;
      offset_upper = sheath_potential_upper - phi_upper;
    }

    // Broadcast the offsets to all processors
    status = MPI_Bcast(&(offset_lower(0,0)), nperp, MPI_DOUBLE, 0, ycomm);
    ASSERT2(status == MPI_SUCCESS);
    status = MPI_Bcast(&(offset_upper(0,0)), nperp, MPI_DOUBLE, nype-1, ycomm);
    ASSERT2(status == MPI_SUCCESS);

    // add average of offsets
    // doesn't matter what 'yindex' of offset_lower and offset_upper is at the moment -
    // just set to be compatible
    offset_lower.setIndex(mesh->ystart);
    offset_upper.setIndex(mesh->ystart);
    FieldPerp average_offset = 0.5*(offset_lower + offset_upper);
    // set location of average_offset to be compatible with phi
    average_offset.setLocation(CELL_CENTRE);
    average_offset = fromFieldAligned(average_offset, "RGN_NOX");
    for (int j=mesh->ystart; j<=mesh->yend; j++) {
      average_offset.setIndex(j);
      phi = average_offset + phi; // this looks weird, but 'FieldPerp + phi' returns a FieldPerp, which is then assigned to an x-z slice of phi
    }
  }
}

void STORM::set_Lx_Ly_Lz() {
  OPTION(options, Lx, -1.);
  if (Lx < 0.) {
    // Check grid spacings and metrics are constant, otherwis this way of calculating is incorrect.
    if (!(
          min(coordinates_centre->dx, true) == max(coordinates_centre->dx, true)
          && min(coordinates_centre->g_11, true) == max(coordinates_centre->g_11, true)
          )) {
      throw BoutException("Error: g_11 or dx is not constant, so cannot calculate Lx here.");
    }
    Lx = (mesh->GlobalNx - 4)*coordinates_centre->dx(mesh->xstart,mesh->ystart)*sqrt(coordinates_centre->g_11(mesh->xstart,mesh->ystart)) ;
    output<<"\tLx input not found, setting to "<<Lx<<endl;
  }
  OPTION(options, Ly, -1.);
  if (Ly < 0.) {
    // Check grid spacings and metrics are constant, otherwis this way of calculating is incorrect.
    if (!(
          min(coordinates_centre->dy, true) == max(coordinates_centre->dy, true)
          && min(coordinates_centre->g_22, true) == max(coordinates_centre->g_22, true)
          )) {
      throw BoutException("Error: g_22 or dy is not constant, so cannot calculate Ly here.");
    }
    Ly = (mesh->GlobalNy - 4)*coordinates_centre->dy(mesh->xstart,mesh->ystart)*sqrt(coordinates_centre->g_22(mesh->xstart,mesh->ystart)) ;
    output<<"\tLy input not found, setting to "<<Ly<<endl;
  }
  OPTION(options, Lz, -1.);
  if (Lz < 0.) {
    // Check grid spacings and metrics are constant, otherwise this way of calculating is incorrect.
    if (!(
          min(coordinates_centre->g_33, true) == max(coordinates_centre->g_33, true)
          )) {
      throw BoutException("Error: g_33 is not constant, so cannot calculate Lz here.");
    }
    Lz = mesh->GlobalNz*coordinates_centre->dz*sqrt(coordinates_centre->g_33(mesh->xstart,mesh->ystart)) ;
    output<<"\tLz input not found, setting to "<<Lz<<endl;
  }
}

Field3D STORM::vort_from_phi(Field3D phi) {
  Field3D result = 0.;
  if (boussinesq > 0) {
    result = Delp2(phi);
  } else {
    result = n*this_Grad_perp2(phi) + this_Grad_perp_dot_Grad_perp(n, phi);
  }
  return result;
}

// Custom version of perpendicular Laplacian operator, which is the inverse of
// the multigrid solver Delp2 uses FFT z-derivatives and Laplace includes
// y-derivatives, so can't use those The function is a copy of Laplace() with
// the y-derivatives deleted
Field3D STORM::this_Grad_perp2(const Field3D &f) {
#if CHECK>0
  // If the x and z derivative schemes are not C2, then this operator is not
  // the inverse of the multigrid Laplacian solver
  std::string first;
  OPTION(globalOptions.getSection("mesh:ddx"), first, "C2");
  ASSERT1(first == "C2");
  std::string second;
  OPTION(globalOptions.getSection("mesh:ddx"), second, "C2");
  ASSERT1(second == "C2");
  OPTION(globalOptions.getSection("mesh:ddz"), first, "FFT");
  ASSERT1(first == "C2");
  OPTION(globalOptions.getSection("mesh:ddz"), second, "FFT");
  ASSERT1(second == "C2");
#endif
  Coordinates* coords = f.getCoordinates();
  Field3D result = coords->G1 * ::DDX(f) +  coords->G3 * ::DDZ(f)
                   + coords->g11 * ::D2DX2(f) + coords->g33 * ::D2DZ2(f)
                   + 2.0 * coords->g13 * ::D2DXDZ(f);

  return result;
}

Field3D STORM::this_Grad_perp_dot_Grad_perp(const Field3D &f, const Field3D &g) {
#if CHECK>0
  // If the x and z derivative schemes are not C2, then this operator is not
  // the inverse of the multigrid Laplacian solver
  std::string first;
  OPTION(globalOptions.getSection("mesh:ddx"), first, "C2");
  ASSERT1(first == "C2");
  std::string second;
  OPTION(globalOptions.getSection("mesh:ddx"), second, "C2");
  ASSERT1(second == "C2");
  OPTION(globalOptions.getSection("mesh:ddz"), first, "FFT");
  ASSERT1(first == "C2");
  OPTION(globalOptions.getSection("mesh:ddz"), second, "FFT");
  ASSERT1(second == "C2");
#endif
  Coordinates* coords = f.getCoordinates();
  Field3D result = coords->g11 * ::DDX(f) * ::DDX(g) + coords->g33 * ::DDZ(f) * ::DDZ(g)
                   + coords->g13 * (DDX(f)*DDZ(g) + DDZ(f)*DDX(g));

  return result;
}

void STORM::set_lower_ddt_zero(Field3D &var){

  // Should only be used on ylow fields
  ASSERT1(var.getLocation() == CELL_YLOW);

  RangeIterator xrup = mesh->iterateBndryLowerY();

  int j =  mesh->ystart;
  for(xrup.first(); !xrup.isDone(); xrup.next()){
    int i = xrup.ind ;
    for(int k=0; k <mesh->LocalNz; k++) {
      ddt(var)(i, j, k) = 0.0 ;
    }
  }
}

void STORM::set_sources_realistic_geometry() {
  S = 0.; S_E = 0.;
  BoutReal x;

  // Core sources
  Options *options = Options::getRoot()->getSection("storm");
  BoutReal S_n0, S_E0, x_Sn, x_SE, w_Sn, w_SE, S_nTar0;
  bool uniform_source_targets, uniform_energy_source_background;
  OPTION(options, S_n0,    0.);
  OPTION(options, S_E0,    0.);
  OPTION(options, x_Sn, 0.417);
  OPTION(options, x_SE, 0.312);
  OPTION(options, w_Sn, 0.025);
  OPTION(options, w_SE, 0.025);
  OPTION(options, S_nTar0, 1.);
  OPTION(options, uniform_source_targets, false);
  OPTION(options, uniform_energy_source_background, false);

  for(int ix = mesh->xstart; ix <= mesh->xend; ++ix){
    for(int iy = mesh->ystart; iy <= mesh->yend; ++iy){
      if((mesh->getGlobalYIndexNoBoundaries(mesh->ystart) > jyseps1_1   && mesh->getGlobalYIndexNoBoundaries(mesh->yend) <= jyseps2_1) ||
         (mesh->getGlobalYIndexNoBoundaries(mesh->ystart) > jyseps1_2-1 && mesh->getGlobalYIndexNoBoundaries(mesh->yend) <= jyseps2_2)){
        // Note: > jyseps1_2-1 to be compatible both with double null and single null
        for(int iz = 0; iz < mesh->LocalNz; ++iz){
          x = ((double)(mesh->getGlobalXIndex(ix))-1.5)/((double)((mesh->GlobalNx - 4)));
          S(ix,iy,iz) = S_n0*exp(-SQ((x-x_Sn)/w_Sn));
          S_E(ix,iy,iz) = S_E0*exp(-SQ((x-x_SE)/w_SE));
        }
      }
    }
  }

  if (sources_realisticgeometry_background) {
    if(realistic_geometry != "singlenull" && realistic_geometry != "doublenull")
      throw BoutException("sources_realisticgeometry_background = true but wrong realistic_geometry.");

    // Option useful to start a simulation, while the turbulence is not yet sustaining the profiles
    BoutReal y, dy;
    if(uniform_energy_source_background) {
        S_E += 1.;
    }
    int ixseps_inner = std::min(ixseps1, ixseps2);
    int ixseps_outer = std::max(ixseps1, ixseps2);
    if(!uniform_source_targets) {
      for(int ix = mesh->xstart; ix <= mesh->xend; ++ix){
        if (mesh->getGlobalXIndex(ix) < ixseps_inner) {
          if (mesh->getGlobalYIndexNoBoundaries(mesh->yend) <= jyseps1_1) {
            // PF region in between 0 and jyseps1_1 (lower inner divertor)
            for(int iy = mesh->ystart; iy <= mesh->yend; ++iy){
              for(int iz = 0; iz < mesh->LocalNz; ++iz){
                dy = 0.5/((double)(jyseps1_1 + 1));
                // 0 < y < 0.5 in this region
                y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy)) + 0.5)*dy;
                S(ix,iy,iz) += S_nTar0*exp(10.*std::abs(y-0.5))/(exp(5.0)-1.0);
              }
            }
          } 
          else if (mesh->getGlobalYIndexNoBoundaries(mesh->ystart) > jyseps2_1 && mesh->getGlobalYIndexNoBoundaries(mesh->yend) < ny_inner) {
            // PF region between jyseps2_1 and ny_inner (upper inner divertor)
            for(int iy = mesh->ystart; iy <= mesh->yend; ++iy){
              for(int iz = 0; iz < mesh->LocalNz; ++iz){
                dy = 0.5/((double)(ny_inner - jyseps2_1 - 1));
                // 0.5 < y < 1.0 in this region
                y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy) - jyseps2_1) - 0.5)*dy + 0.5;
                S(ix,iy,iz) += S_nTar0*exp(10.*std::abs(y-0.5))/(exp(5.0)-1.0);
              }
            }
          }
          else if(mesh->getGlobalYIndexNoBoundaries(mesh->ystart) >= ny_inner &&  mesh->getGlobalYIndexNoBoundaries(mesh->yend) <= jyseps1_2) {
            // PF region between ny_inner and jyseps1_2 (upper outer divertor)
            for(int iy = mesh->ystart; iy <= mesh->yend; ++iy){
              for(int iz = 0; iz < mesh->LocalNz; ++iz){
                dy = 0.5/((double)(jyseps1_2 - ny_inner + 1));
                // 0 < y < 0.5 in this region
                y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy)) - ny_inner + 0.5)*dy;
                S(ix,iy,iz) += S_nTar0*exp(10.*std::abs(y-0.5))/(exp(5.0)-1.0);
              }
            }
          }
          else if(mesh->getGlobalYIndexNoBoundaries(mesh->ystart) > jyseps2_2) {
            // PF region between jyseps2_2 and NyGlobal (lower outer divertor)
            for(int iy = mesh->ystart; iy <= mesh->yend; ++iy){
              for(int iz = 0; iz < mesh->LocalNz; ++iz){
                dy = 0.5/((double)(mesh->GlobalNy -2*mesh->ystart - jyseps2_2 - 1));
                // 0.5 < y < 1.0 in this region
                y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy) - jyseps2_2) - 0.5)*dy + 0.5;
                S(ix,iy,iz) += S_nTar0*exp(10.*std::abs(y-0.5))/(exp(5.0)-1.0);
              }
            }
          }
        }
        else if (mesh->getGlobalXIndex(ix) >= ixseps_outer) {
          if (mesh->getGlobalYIndexNoBoundaries(mesh->yend) < ny_inner) {
            // SOL region on inner side
            for(int iy = mesh->ystart; iy <= mesh->yend; ++iy){
              for(int iz = 0; iz < mesh->LocalNz; ++iz){
                dy = 1./((double)(ny_inner));
                // 0 < y < 1 in this region
                y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy)) + 0.5)*dy;
                S(ix,iy,iz) += S_nTar0*exp(10.*std::abs(y-0.5))/(exp(5.0)-1.0);
              }
            }
          }
          else {
            // SOL region on outer side
            for(int iy = mesh->ystart; iy <= mesh->yend; ++iy){
              for(int iz = 0; iz < mesh->LocalNz; ++iz){
                dy = 1./((double)(mesh->GlobalNy -2*mesh->ystart - ny_inner));
                // 0 < y < 1 in this region
                y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy)) - ny_inner + 0.5)*dy;
                S(ix,iy,iz) += S_nTar0*exp(10.*std::abs(y-0.5))/(exp(5.0)-1.0);
              }
            }
          }
        }
        else if (ixseps1< ixseps2){
          // inter-separatrix for lower double null
          if (mesh->getGlobalYIndexNoBoundaries(mesh->yend) <= jyseps2_1) {
            // inter-separatrix including lower inner divertor leg
            for(int iy = mesh->ystart; iy <= mesh->yend; ++iy){
              for(int iz = 0; iz < mesh->LocalNz; ++iz){
                dy = 0.5/((double)(jyseps2_1 + 1));
                // 0 < y < 0.5 in this region
                y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy)) + 0.5)*dy;
                S(ix,iy,iz) += S_nTar0*exp(10.*std::abs(y-0.5))/(exp(5.0)-1.0);
              }
            }
          }
          else if (mesh->getGlobalYIndexNoBoundaries(mesh->ystart) > jyseps2_1 && mesh->getGlobalYIndexNoBoundaries(mesh->yend) < ny_inner) {
            // inter-separatrix in upper inner divertor leg
            for(int iy = mesh->ystart; iy <= mesh->yend; ++iy){
              for(int iz = 0; iz < mesh->LocalNz; ++iz){
                dy = 0.5/((double)(ny_inner - jyseps2_1 - 1));
                // 0.5 < y < 1.0 in this region
                y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy) - jyseps2_1) - 0.5)*dy + 0.5;
                S(ix,iy,iz) += S_nTar0*exp(10.*std::abs(y-0.5))/(exp(5.0)-1.0);
              }
            }
          }
          else if (mesh->getGlobalYIndexNoBoundaries(mesh->ystart) >= ny_inner &&  mesh->getGlobalYIndexNoBoundaries(mesh->yend) <= jyseps1_2) {
            // inter-separatrix in upper outer divertor leg
            for(int iy = mesh->ystart; iy <= mesh->yend; ++iy){
              for(int iz = 0; iz < mesh->LocalNz; ++iz){
                dy = 0.5/((double)(jyseps1_2 - ny_inner + 1));
                // 0 < y < 0.5 in this region
                y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy)) - ny_inner + 0.5)*dy;
                S(ix,iy,iz) += S_nTar0*exp(10.*std::abs(y-0.5))/(exp(5.0)-1.0);
              }
            }
          }
          else if (mesh->getGlobalYIndexNoBoundaries(mesh->ystart) > jyseps1_2-1) {
            // inter-separatrix including lower outer divertor leg
            for(int iy = mesh->ystart; iy <= mesh->yend; ++iy){
              for(int iz = 0; iz < mesh->LocalNz; ++iz){
                dy = 0.5/((double)(mesh->GlobalNy -2*mesh->ystart - jyseps1_2 - 1));
                // 0.5 < y < 1.0 in this region
                y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy) - jyseps1_2) - 0.5)*dy + 0.5;
                S(ix,iy,iz) += S_nTar0*exp(10.*std::abs(y-0.5))/(exp(5.0)-1.0);
              }
            }
          }
        }
        else if (ixseps1> ixseps2){
          //inter-separatrix for upper double null
          if (mesh->getGlobalYIndexNoBoundaries(mesh->yend) <= jyseps1_1) {
            // inter-separatrix in lower inner divertor leg
            for(int iy = mesh->ystart; iy <= mesh->yend; ++iy){
              for(int iz = 0; iz < mesh->LocalNz; ++iz){
                dy = 0.5/((double)(jyseps1_1 + 1));
                // 0 < y < 0.5 in this region
                y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy)) + 0.5)*dy;
                S(ix,iy,iz) += S_nTar0*exp(10.*std::abs(y-0.5))/(exp(5.0)-1.0);
              }
            }
          }
          else if (mesh->getGlobalYIndexNoBoundaries(mesh->ystart) > jyseps1_1 && mesh->getGlobalYIndexNoBoundaries(mesh->yend) < ny_inner) {
            // inter-separatrix including upper inner divertor leg
            for(int iy = mesh->ystart; iy <= mesh->yend; ++iy){
              for(int iz = 0; iz < mesh->LocalNz; ++iz){
                dy = 0.5/((double)(ny_inner - jyseps1_1 - 1));
                // 0.5 < y < 1.0 in this region
                y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy) - jyseps1_1) - 0.5)*dy + 0.5;
                S(ix,iy,iz) += S_nTar0*exp(10.*std::abs(y-0.5))/(exp(5.0)-1.0);
              }
            }
          }
          else if (mesh->getGlobalYIndexNoBoundaries(mesh->ystart) >= ny_inner &&  mesh->getGlobalYIndexNoBoundaries(mesh->yend) <= jyseps2_2) {
            // inter-separatrix including upper outer divertor leg
            for(int iy = mesh->ystart; iy <= mesh->yend; ++iy){
              for(int iz = 0; iz < mesh->LocalNz; ++iz){
                dy = 0.5/((double)(jyseps2_2 - ny_inner + 1));
                // 0 < y < 0.5 in this region
                y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy)) - ny_inner + 0.5)*dy;
                S(ix,iy,iz) += S_nTar0*exp(10.*std::abs(y-0.5))/(exp(5.0)-1.0);
              }
            }
          }
          else if (mesh->getGlobalYIndexNoBoundaries(mesh->ystart) > jyseps2_2) {
            // inter-separatrix in lower outer divertor leg
            for(int iy = mesh->ystart; iy <= mesh->yend; ++iy){
              for(int iz = 0; iz < mesh->LocalNz; ++iz){
                dy = 0.5/((double)(mesh->GlobalNy -2*mesh->ystart - jyseps2_2 - 1));
                // 0.5 < y < 1.0 in this region
                y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy) - jyseps2_2) - 0.5)*dy + 0.5;
                S(ix,iy,iz) += S_nTar0*exp(10.*std::abs(y-0.5))/(exp(5.0)-1.0);
              }
            }
          }
        }
      }
    }
    else { // uniform S distribution in the divertor legs, same as the PF regions above
      for(int ix = mesh->xstart; ix <= mesh->xend; ++ix){
        if (mesh->getGlobalYIndexNoBoundaries(mesh->yend) <= jyseps1_1) {
          // between 0 and jyseps1_1 (lower inner divertor)
          for(int iy = mesh->ystart; iy <= mesh->yend; ++iy){
            for(int iz = 0; iz < mesh->LocalNz; ++iz){
              dy = 0.5/((double)(jyseps1_1 + 1));
              y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy)) + 0.5)*dy;
              S(ix,iy,iz) += S_nTar0*exp(10.*std::abs(y-0.5))/(exp(5.0)-1.0);
            }
          }
        } 
        else if (mesh->getGlobalYIndexNoBoundaries(mesh->ystart) > jyseps2_1 && mesh->getGlobalYIndexNoBoundaries(mesh->yend) < ny_inner) {
          // between jyseps2_1 and ny_inner (upper inner divertor)
          for(int iy = mesh->ystart; iy <= mesh->yend; ++iy){
            for(int iz = 0; iz < mesh->LocalNz; ++iz){
              dy = 0.5/((double)(ny_inner - jyseps2_1 - 1));
              y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy) - jyseps2_1) - 0.5)*dy + 0.5;
              S(ix,iy,iz) += S_nTar0*exp(10.*std::abs(y-0.5))/(exp(5.0)-1.0);
            }
          }
        }
        else if(mesh->getGlobalYIndexNoBoundaries(mesh->ystart) >= ny_inner &&  mesh->getGlobalYIndexNoBoundaries(mesh->yend) <= jyseps1_2) {
          // between ny_inner and jyseps1_2 (upper outer divertor)
          for(int iy = mesh->ystart; iy <= mesh->yend; ++iy){
            for(int iz = 0; iz < mesh->LocalNz; ++iz){
              dy = 0.5/((double)(jyseps1_2 - ny_inner + 1));
              y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy)) - ny_inner + 0.5)*dy;
              S(ix,iy,iz) += S_nTar0*exp(10.*std::abs(y-0.5))/(exp(5.0)-1.0);
            }
          }
        }
        else if(mesh->getGlobalYIndexNoBoundaries(mesh->ystart) > jyseps2_2) {
          // between jyseps2_2 and NyGlobal (lower outer divertor)
          for(int iy = mesh->ystart; iy <= mesh->yend; ++iy){
            for(int iz = 0; iz < mesh->LocalNz; ++iz){
              dy = 0.5/((double)(mesh->GlobalNy -2*mesh->ystart - jyseps2_2 - 1));
              y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy) - jyseps2_2) - 0.5)*dy + 0.5;
              S(ix,iy,iz) += S_nTar0*exp(10.*std::abs(y-0.5))/(exp(5.0)-1.0);
            }
          }
        }
      }
    }
  }

  S /= 13500.;
  S_E /= 13500.;

  mesh->communicate(S);
  S.applyBoundary("free_o3");
  S_stag = interp_to(S, CELL_YLOW, "RGN_NOBNDRY");

  if (isothermal) S_E = 0.;
}

void STORM::setup_history_tracking(bool restarting, const std::string& storm_git_hash) {
  // track history of inputs over restarts of a simulation

  int restart_first_iteration;
  if (not restarting) {
    restart_first_iteration = 0;

    // Resize here so that std::vector storage does not get reallocated after being passed
    // to restart.addOnce() and dump.addOnce().
    restart_at_iteration_history.resize(1);
    input_files_history.resize(1);
  } else {
    // Open restart file directly with NcFile, because BOUT++'s Datafile interface is too
    // limited. We only read stuff here, and the restart files should not currently be
    // open, so this should not cause any problems.
    std::stringstream restart_file_name;
    restart_file_name << globalOptions["datadir"].as<std::string>() << "/BOUT.restart."
                      << BoutComm::rank() << ".nc";
    netCDF::NcFile restart_reader(restart_file_name.str(), netCDF::NcFile::read);

    {
      auto var = restart_reader.getVar("restart_counter");
      var.getVar(&restart_counter);
    }

    {
      auto var = restart_reader.getVar("hist_hi");
      var.getVar(&restart_first_iteration);
    }

    // Resize here so that std::vector storage does not get reallocated after being passed
    // to restart.addOnce() and dump.addOnce().
    restart_at_iteration_history.resize(restart_counter+1);
    input_files_history.resize(restart_counter+1);

    for (int i=0; i < restart_counter; i++) {
      std::stringstream iteration_history_name;
      iteration_history_name << "restart_at_iteration_history" << i;

      restart.addOnce(restart_at_iteration_history[i], iteration_history_name.str());
      dump.addOnce(restart_at_iteration_history[i], iteration_history_name.str());

      std::stringstream input_history_name;
      input_history_name << "input_files_history" << i;

      // Have to read the input_files_history* variables before adding to `restart` and
      // `dump` because string lengths are not allowed to change after being added.
      {
        auto var = restart_reader.getVar(input_history_name.str());
        auto size = var.getDim(0).getSize();

        if (size > 0) {
          // Have to read as C-style char array
          char* this_input_contents = new char[size];
          var.getVar(this_input_contents);

          input_files_history[i] = std::string(this_input_contents, size);

          delete [] this_input_contents;
        } else {
          // existing value was empty string
          input_files_history[i] = "";
        }
      }

      restart.addOnce(input_files_history[i], input_history_name.str());
      dump.addOnce(input_files_history[i], input_history_name.str());
    }
  }

  restart_at_iteration_history[restart_counter] = restart_first_iteration;
  std::stringstream iteration_history_name;
  iteration_history_name << "restart_at_iteration_history" << restart_counter;
  dump.addOnce(restart_at_iteration_history[restart_counter], iteration_history_name.str());

  // Add current input file
  std::stringstream input_file_name;
  input_file_name << globalOptions["datadir"].as<std::string>() << "/"
                  << globalOptions["optionfile"].as<std::string>() ;
  std::ostringstream input_contents;
  input_contents << "Restart using STORM version " << storm_git_hash << endl
    << "----------------------------------------------------------------------------------------"
    << endl << endl;
  {
    // Implementation borrowed from
    // https://insanecoding.blogspot.com/2011/11/how-to-read-in-file-in-c.html
    std::ifstream in(input_file_name.str(), std::ios::in | std::ios::binary);
    if (in) {
      input_contents << in.rdbuf();
      in.close();
    } else {
      throw BoutException("Failed to open input file %s", input_file_name.str().c_str());
    }
  }
  input_files_history[restart_counter] = input_contents.str();
  std::stringstream input_history_name;
  input_history_name << "input_files_history" << restart_counter;
  dump.addOnce(input_files_history[restart_counter], input_history_name.str());

  restart.addOnce(restart_counter, "restart_counter");
  dump.addOnce(restart_counter, "restart_counter");

}

void STORM::history_tracking_first_rhs() {
  // Add new variables to restart files here, to be called after restart files have been
  // read, to avoid needing to use `restart:init_missing=true` because they were not
  // present in old restart files.
  std::stringstream iteration_history_name;
  iteration_history_name << "restart_at_iteration_history" << restart_counter;
  restart.addOnce(restart_at_iteration_history[restart_counter], iteration_history_name.str());

  std::stringstream input_history_name;
  input_history_name << "input_files_history" << restart_counter;
  restart.addOnce(input_files_history[restart_counter], input_history_name.str());

  restart_counter++;
}

void STORM::enhance_in_radial_buffers(Field2D& coefficient,
                                      const int buffer_size,
                                      const BoutReal peak_enhancement) {
  // Inner radial boundary
  for (int xglobal = mesh->getGlobalXIndex(0);
       xglobal < std::min(buffer_size + mesh->xstart,
                          mesh->getGlobalXIndex(mesh->LocalNx));
       xglobal++) {

    int xlocal = mesh->getLocalXIndex(xglobal);
    for (int y = 0; y < mesh->LocalNy; y++) {
      coefficient(xlocal, y) +=
          BoutReal(buffer_size + mesh->xstart - xglobal)
          / BoutReal(buffer_size)
          * (peak_enhancement - 1.0) * coefficient(xlocal, y);
    }
  }

  // Outer radial boundary
  for (int xglobal = std::max(mesh->GlobalNx-1-buffer_size-mesh->xstart,
                              mesh->getGlobalXIndex(0));
       xglobal < mesh->getGlobalXIndex(mesh->LocalNx); xglobal++) {

    int xlocal = mesh->getLocalXIndex(xglobal);
    for (int y = 0; y < mesh->LocalNy; y++) {
      coefficient(xlocal, y) +=
          BoutReal(buffer_size + xglobal - (mesh->GlobalNx - 1 - mesh->xstart))
          / BoutReal(buffer_size)
          * (peak_enhancement - 1.0) * coefficient(xlocal, y);
    }
  }
}
