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
  if (boussinesq) {
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
  OPTION(globalOptions->getSection("mesh:ddx"), first, "C2");
  ASSERT1(first == "C2");
  std::string second;
  OPTION(globalOptions->getSection("mesh:ddx"), second, "C2");
  ASSERT1(second == "C2");
  OPTION(globalOptions->getSection("mesh:ddz"), first, "FFT");
  ASSERT1(first == "C2");
  OPTION(globalOptions->getSection("mesh:ddz"), second, "FFT");
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
  OPTION(globalOptions->getSection("mesh:ddx"), first, "C2");
  ASSERT1(first == "C2");
  std::string second;
  OPTION(globalOptions->getSection("mesh:ddx"), second, "C2");
  ASSERT1(second == "C2");
  OPTION(globalOptions->getSection("mesh:ddz"), first, "FFT");
  ASSERT1(first == "C2");
  OPTION(globalOptions->getSection("mesh:ddz"), second, "FFT");
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
