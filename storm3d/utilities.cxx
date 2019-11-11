// Utility methods for STORM

#include "storm.hxx"
#include <derivs.hxx>
#include <interpolation.hxx>

void STORM::phisolver_1d(){
  // Assume adiabaticity to get parallel profile of phi
  // phi = mu/(mu+1.)*log(n) for
  if (mesh->getNYPE()>1) throw BoutException("Cannot run in parallel when run_1d=true");
  if (mesh->GlobalNz>1) throw BoutException("Cannot use more than 1 z point when run_1d=true");
  if (mesh->GlobalNx>5) throw BoutException("Cannot use more than 1 x point when run_1d=true");

  // Apply U boundary conditions before solving for phi
  if (!symmetry_plane) {
    Usheath_ydown_staggered(U, sqrt(T_sheath_lower), Usheath_BC_prefactor);
  }
  Usheath_yup_staggered(U, sqrt(T_sheath_upper), Usheath_BC_prefactor);

  Field3D gradpar_phi;
  if(isothermal) {
    gradpar_phi = mu/(mu+1.)*Grad_par(n,CELL_YLOW)/interp_to(n,CELL_YLOW);
  } else {
    gradpar_phi = mu/(mu+1.)*Grad_par(n*T,CELL_YLOW)/interp_to(n,CELL_YLOW) + 0.71*Grad_par(T,CELL_YLOW);
  }
  if (symmetry_plane) {
    phi(mesh->xstart,mesh->yend,0) = 0;
    for (int j=mesh->yend; j>mesh->ystart; j--) {
      phi(mesh->xstart,j-1,0) = phi(mesh->xstart,j,0) - gradpar_phi(mesh->xstart,j,0)*mesh->getCoordinates()->dy(mesh->xstart,j-1)*sqrt(mesh->getCoordinates()->g_22(mesh->xstart,j-1));
    }
    BoutReal phi_upper;
    phi.applyBoundary();
    phi_upper = extrap_sheath_upper(phi)(mesh->xstart,0);

    // Offset phi so that the boundary value is the sheath value
    phi -= phi_upper;
    phi += phi_wall + log(Vsheath_BC_prefactor*sqrt(T_sheath_upper(mesh->xstart,0))/U(mesh->xstart,mesh->yend+1,0))*T_sheath_upper(mesh->xstart,0);
  } else {
    // Integrate symmetrically, one pass up and one down
    Field3D phi_up = 0.;
    Field3D phi_down = 0.;
    phi_up(mesh->xstart,mesh->ystart,0) = 0.;
    for (int j=mesh->ystart; j<mesh->yend; j++) {
      phi_up(mesh->xstart,j+1,0) = phi_up(mesh->xstart,j,0) + gradpar_phi(mesh->xstart,j+1,0)*mesh->getCoordinates()->dy(mesh->xstart,j)*sqrt(mesh->getCoordinates()->g_22(mesh->xstart,j));
    }
    phi_down(mesh->xstart,mesh->yend,0) = 0.;
    for (int j=mesh->yend; j>mesh->ystart; j--) {
      phi_down(mesh->xstart,j-1,0) = phi_down(mesh->xstart,j,0) - gradpar_phi(mesh->xstart,j,0)*mesh->getCoordinates()->dy(mesh->xstart,j-1)*sqrt(mesh->getCoordinates()->g_22(mesh->xstart,j-1));
    }
    phi = (phi_up+phi_down)/2.;
    BoutReal phi_lower, phi_upper;
    phi.applyBoundary();
    phi_lower = extrap_sheath_lower(phi)(mesh->xstart,0);
    phi_upper = extrap_sheath_upper(phi)(mesh->xstart,0);

    // Offset phi so that the average of the sheath boundary values is zero
    phi -= (phi_lower+phi_upper)/2.;

    // Offset phi relative to phi_wall
    phi += 0.5*(log(Vsheath_BC_prefactor*sqrt(T_sheath_upper(mesh->xstart,0))/U(mesh->xstart,mesh->yend+1,0))*T_sheath_upper(mesh->xstart,0)
                + log(Vsheath_BC_prefactor*sqrt(T_sheath_lower(mesh->xstart,0))/(-U(mesh->xstart,mesh->ystart,0)))*T_sheath_lower(mesh->xstart,0))
           + phi_wall;
  }
}

////////////////////////////////////////////////////////////////////////

void STORM::set_Lx_Ly_Lz() {
  OPTION(options, Lx, -1.);
  if (Lx < 0.) {
    // Check grid spacings and metrics are constant, otherwis this way of calculating is incorrect.
    if (!(
          min(mesh->getCoordinates()->dx, true) == max(mesh->getCoordinates()->dx, true)
          && min(mesh->getCoordinates()->g_11, true) == max(mesh->getCoordinates()->g_11, true)
          )) {
      throw BoutException("Error: g_11 or dx is not constant, so cannot calculate Lx here.");
    }
    Lx = (mesh->GlobalNx - 4)*mesh->getCoordinates()->dx(mesh->xstart,mesh->ystart)*sqrt(mesh->getCoordinates()->g_11(mesh->xstart,mesh->ystart));
    output<<"\tLx input not found, setting to "<<Lx<<endl;
  }
  OPTION(options, Ly, -1.);
  if (Ly < 0.) {
    // Check grid spacings and metrics are constant, otherwis this way of calculating is incorrect.
    if (!(
          min(mesh->getCoordinates()->dy, true) == max(mesh->getCoordinates()->dy, true)
          && min(mesh->getCoordinates()->g_22, true) == max(mesh->getCoordinates()->g_22, true)
          )) {
      throw BoutException("Error: g_22 or dy is not constant, so cannot calculate Ly here.");
    }
    Ly = (mesh->GlobalNy - 4)*mesh->getCoordinates()->dy(mesh->xstart,mesh->ystart)*sqrt(mesh->getCoordinates()->g_22(mesh->xstart,mesh->ystart));
    output<<"\tLy input not found, setting to "<<Ly<<endl;
  }
  OPTION(options, Lz, -1.);
  if (Lz < 0.) {
    // Check grid spacings and metrics are constant, otherwis this way of calculating is incorrect.
    if (!(
          min(mesh->getCoordinates()->g_33, true) == max(mesh->getCoordinates()->g_33, true)
          )) {
      throw BoutException("Error: g_33 is not constant, so cannot calculate Lz here.");
    }
    Lz = mesh->GlobalNz*mesh->getCoordinates()->dz*sqrt(mesh->getCoordinates()->g_33(mesh->xstart,mesh->ystart)) ;
    output<<"\tLz input not found, setting to "<<Lz<<endl;
  }
}

