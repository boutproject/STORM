
#include "D-vpar.hxx"

#include <boutexception.hxx>
#include <derivs.hxx>
#include <bout/constants.hxx>

NeutralDVpar::NeutralDVpar(Solver *solver, Options &options, Datafile &dump)
    : NeutralModel(options, dump), minmaxmean_monitor(*this) {

  mesh = bout::globals::mesh;
  coordinates_centre = mesh->getCoordinates(CELL_CENTRE);
  coordinates_stag = mesh->getCoordinates(CELL_YLOW);
  OPTION(options, Tn0,                 3.5);
  OPTION(options, munvn,                -1);
  OPTION(options, dt_update,           0.1);
  OPTION(options, monitor_minmaxmean, true);
  OPTION(options, lambdamax,            -1);
  OPTION(options, bndry_neutrals,        2); // tried different options, this is the one working better
  OPTION(options, Rc,                  0.8);
  OPTION(options, yderivatives,      false);

  updaterates = true;

  if (mesh->xstart > 1) {
    orderupwindscheme = "U2";
  }
  else {
    output << "WARNING, ONLY ONE GUARD CELL IN X => USING U1 SCHEME IN X!" << endl;
    orderupwindscheme = "U1";
  }

  auto& opt = *Options::getRoot();
  opt["lognn"]["bndry_xin"] = opt["lognn"]["bndry_xin"].withDefault("neumann_o2");
  opt["lognn"]["bndry_xout"] = opt["lognn"]["bndry_xout"].withDefault("neumann_o2");
  opt["lognn"]["bndry_ydown"] = opt["lognn"]["bndry_ydown"].withDefault("none");
  opt["lognn"]["bndry_yup"] = opt["lognn"]["bndry_yup"].withDefault("none");

  opt["lognn_aligned"]["bndry_xin"] = opt["lognn_aligned"]["bndry_xin"].withDefault("none");
  opt["lognn_aligned"]["bndry_xout"] = opt["lognn_aligned"]["bndry_xout"].withDefault("none");
  if (bndry_neutrals >= 2) {
    opt["lognn_aligned"]["bndry_ydown"] = opt["lognn_aligned"]["bndry_ydown"].withDefault("neumann_o2");
    opt["lognn_aligned"]["bndry_yup"] = opt["lognn_aligned"]["bndry_yup"].withDefault("neumann_o2");
  }
  else {
    opt["lognn_aligned"]["bndry_ydown"] = opt["lognn_aligned"]["bndry_ydown"].withDefault("none");
    opt["lognn_aligned"]["bndry_yup"] = opt["lognn_aligned"]["bndry_yup"].withDefault("none");
  }

  opt["nvn"]["bndry_xin"] = opt["nvn"]["bndry_xin"].withDefault("neumann_o2");
  opt["nvn"]["bndry_xout"] = opt["nvn"]["bndry_xout"].withDefault("neumann_o2");
  opt["nvn"]["bndry_ydown"] = opt["nvn"]["bndry_ydown"].withDefault("none");
  opt["nvn"]["bndry_yup"] = opt["nvn"]["bndry_yup"].withDefault("none");

  opt["nvn_aligned"]["bndry_xin"] = opt["nvn_aligned"]["bndry_xin"].withDefault("none");
  opt["nvn_aligned"]["bndry_xout"] = opt["nvn_aligned"]["bndry_xout"].withDefault("none");
  if (bndry_neutrals == 2 || bndry_neutrals == 5) {
    opt["nvn_aligned"]["bndry_ydown"] = opt["nvn_aligned"]["bndry_ydown"].withDefault("neumann_o2");
    opt["nvn_aligned"]["bndry_yup"] = opt["nvn_aligned"]["bndry_yup"].withDefault("neumann_o2");
  }
  else {
    opt["nvn_aligned"]["bndry_ydown"] = opt["nvn_aligned"]["bndry_ydown"].withDefault("none");
    opt["nvn_aligned"]["bndry_yup"] = opt["nvn_aligned"]["bndry_yup"].withDefault("none");
  }

  lognn_stag.setLocation(CELL_YLOW);
  nn_stag.setLocation(CELL_YLOW);
  nvn.setLocation(CELL_YLOW);
  Tn_stag.setLocation(CELL_YLOW);
  Snvn.setLocation(CELL_YLOW);

  SOLVE_FOR(lognn);
  SOLVE_FOR(nvn);

  // Initialize fields
  nn = 0.001;
  nvn = 0.;
  lognn = log(nn);
  lognn_stag = log(0.001);
  nn_stag = 0.001;
  Dn = 1.;
  recycled_lower = 0.;
  recycled_upper = 0.;
  vn_lower = 0.;
  vn_upper = 0.;
  Snn = 0.;
  Snvn = 0.;
  vn = nvn/nn_stag;

  BoutReal Snn0;
  OPTION(options, Snn0, 0.);
  for(int ix = mesh->xstart; ix <= mesh->xend; ++ix){
    if (!mesh->periodicY(ix)) {
      for(int iy = mesh->ystart; iy <= mesh->yend; ++iy){
        for(int iz = 0; iz < mesh->LocalNz; ++iz){
          BoutReal x = ((double)(mesh->getGlobalXIndex(ix))-1.5)/((double)((mesh->GlobalNx - 4)));
          Snn(ix,iy,iz) += Snn0*exp(-SQ(x)/0.002)/15000.;
          Snn(ix,iy,iz) += Snn0*exp(-SQ(x-1.)/0.002)/15000.;
        }
      }
    }
    else {
      for(int iy = mesh->ystart; iy <= mesh->yend; ++iy){
        for(int iz = 0; iz < mesh->LocalNz; ++iz){
          BoutReal x = ((double)(mesh->getGlobalXIndex(ix))-1.5)/((double)((mesh->GlobalNx - 4)));
          Snn(ix,iy,iz) += 0.00000002*exp(-SQ(x)/0.008);
        }
      }
    }
  }

  nn_lower = 0., nn_upper = 0.;
  nn_lower.setIndex(mesh->ystart);
  nn_lower.setLocation(CELL_YLOW);
  nn_lower.setDirectionY(YDirectionType::Aligned);
  nn_upper.setIndex(mesh->yend+1);
  nn_upper.setLocation(CELL_YLOW);
  nn_upper.setDirectionY(YDirectionType::Aligned);
  vn_lower.setIndex(mesh->ystart);
  vn_lower.setLocation(CELL_YLOW);
  vn_lower.setDirectionY(YDirectionType::Aligned);
  vn_upper.setIndex(mesh->yend+1);
  vn_upper.setLocation(CELL_YLOW);
  vn_upper.setDirectionY(YDirectionType::Aligned);

  lognn_aligned.setBoundary("lognn_aligned");
  nvn_aligned.setBoundary("nvn_aligned");
  vn.setBoundary("nvn");
  Dn.setBoundary("lognn");

  if (bndry_neutrals != 1 && bndry_neutrals != 2 && bndry_neutrals != 3 && bndry_neutrals != 4 && bndry_neutrals != 5)
    throw BoutException("Unrecognized option bndry_neutrals!");
  if (Rc < 0. || Rc > 1.)
    throw BoutException("Rc must be in [0,1]");
  if (bndry_neutrals != 4 && bndry_neutrals != 5) {
    SAVE_ONCE2(Snn,Snvn);
  }
  else {
    initialiseSource(options);
    SAVE_REPEAT2(Snn,Snvn);
  }
  SAVE_ONCE(Tn);
  SAVE_REPEAT(nn);
  SAVE_REPEAT(vn);
  SAVE_REPEAT(Dn);

  if (yderivatives) {
    DDY_Jogyy_oJ = DDY(coordinates_centre->J/coordinates_centre->g_22)/coordinates_centre->J;
    g_12og_22 = coordinates_centre->g_12/coordinates_centre->g_22;
    g_23og_22 = coordinates_centre->g_23/coordinates_centre->g_22;
  }
  
  if (neutral_dvpar_instances.size() == 0) {
    // No NeutralDVpar instances have been created before, so this is the first one: Need
    // to add the timestep monitor function to the solver.
    solver->addTimestepMonitor(timestepmonitor_func);
  }
  // Add this object to the global list, so it can be found by timestepmonitor_func
  neutral_dvpar_instances.push_back(this);

  if (monitor_minmaxmean) {
    solver->addMonitor(&minmaxmean_monitor, MonitorPosition::BACK);
  }
}

// 3D model, diffusive in X-Z, fluid in Y
void NeutralDVpar::update(const Field3D& n, const Field3D& n_stag, const Field3D& U, 
                          const Field3D& V, const Field3D& T, const Field3D& T_stag,
                          BoutReal time) {
  TRACE("NeutralDVpar::update");
  rhs_counter++;

  Tn = Tn0/T0;
  Tn_stag = Tn0/T0;

  // lognn_aligned, nvn_aligned
  lognn_aligned = toFieldAligned(lognn, "RGN_NOBNDRY");
  nvn_aligned =  toFieldAligned(nvn, "RGN_NOBNDRY");
  mesh->communicate(lognn, nvn, lognn_aligned, nvn_aligned);
  
  // Apply boundary conditions in y to logn_aligend and nvn_aligned, compute vn_lower and vn_upper
  recycleFluxes(time);

  // nn, nn_aligned, nn_stag
  nn = exp(lognn, "RGN_NOY");
  nn_aligned = exp(lognn_aligned, "RGN_NOX");
  lognn_stag = myInterp_to(lognn_aligned, CELL_YLOW);
  nn_stag = exp(lognn_stag, "RGN_NOBNDRY");

  // vn, vn_aligned
  vn = nvn/nn_stag;
  vn_aligned = toFieldAligned(vn, "RGN_NOBNDRY");
  mesh->communicate(vn, vn_aligned);
  vn.applyBoundary(time);
  int i, j = mesh->ystart;
  for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
    i = r.ind;
    for (int k = 0; k < mesh->LocalNz; k++) {
      vn_aligned(i,j,k) = vn_lower(i,k);
      vn_aligned(i,j-1,k) = 3.*vn_aligned(i,j,k) + 3.*vn_aligned(i,j+1,k) - vn_aligned(i,j+2,k);
      vn_aligned(i,j-2,k) = 3.*vn_aligned(i,j-1,k) + 3.*vn_aligned(i,j,k) - vn_aligned(i,j+1,k);
    }
  }
  j = mesh->yend+1;
  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    i = r.ind;
    for (int k = 0; k < mesh->LocalNz; k++) {
      vn_aligned(i,j,k) = vn_upper(i,k);
      vn_aligned(i,j+1,k) = 3.*vn_aligned(i,j,k) + 3.*vn_aligned(i,j-1,k) - vn_aligned(i,j-2,k);
    }
  }
  vn_centre = myInterp_to(vn_aligned, CELL_CENTRE);

  // Update collision frequencies and reaction rates
  neutral_rates(n, n_stag, U, V, T, T_stag, nn, nn_stag, vn, updaterates);
  
  // Compute the diffusion coefficients
  // Use hard-sphere collisions for neutrals, see e.g. https://en.wikipedia.org/wiki/Cross_section_(physics)#Collision_among_gas_particles
  const BoutReal a0 = PI*SQ(2.*5.29e-11); // pi*(2*r0)^2, where r0 is the Bohr radius
  BOUT_FOR(i, Dn.getRegion("RGN_NOBNDRY")) {
    BoutReal lambda_nn = 1./(nn[i]*a0*n0);
    if(lambda_nn > lambdamax && lambdamax > 0.) lambda_nn = lambdamax;
    lambda_nn /= L0;
    Dn[i] = Tn[i]/(sqrt(Tn[i])/lambda_nn + nuiz[i] + nucx[i]);
  }
  mesh->communicate(Dn);
  Dn.applyBoundary(time);

  // Dn*nn*Grad_par(vn)
  Field3D Dn_nn_Gradpar_vn_aligned = toFieldAligned( Dn*nn*myGrad_par(vn_aligned, CELL_CENTRE) , "RGN_NOBNDRY");

  // LaplPerp(nn) and Grad_perp(Dn)*Grad_perp(log(nn))
  Field3D LaplPerpnn = Delp2(nn);
  Field3D GradPerp_Dn_x = DDX(Dn), GradPerp_Dn_z = DDZ(Dn);
  Field3D GradPerp_nn_x = DDX(nn), GradPerp_nn_z = DDZ(nn);
  if (yderivatives) {
    Field3D Dn_aligned =  toFieldAligned(Dn, "RGN_NOBNDRY");
    mesh->communicate(Dn_aligned);
    Dn_aligned.applyBoundary("free_o3");
    Field3D DDY_Dn = fromFieldAligned(DDY(Dn_aligned), "RGN_NOBNDRY"), DDY_nn = fromFieldAligned(DDY(nn_aligned), "RGN_NOBNDRY");

    GradPerp_Dn_x -= g_12og_22*DDY_Dn;
    GradPerp_Dn_z -= g_23og_22*DDY_Dn;
    GradPerp_nn_x -= g_12og_22*DDY_nn;
    GradPerp_nn_z -= g_23og_22*DDY_nn;

    // assumes g12=g13=0
    LaplPerpnn += (coordinates_centre->G2 - DDY_Jogyy_oJ)*DDY_nn 
                + (coordinates_centre->g22 - 1./coordinates_centre->g_22)*fromFieldAligned(D2DY2(nn_aligned, CELL_CENTRE, "C2"), "RGN_NOBNDRY")
                + 2.*coordinates_centre->g23*DDZ(DDY_nn);
  }
  Field3D GradPerpDn_GradPerpnn = coordinates_centre->g11*GradPerp_Dn_x*GradPerp_nn_x 
                                + coordinates_centre->g33*GradPerp_Dn_z*GradPerp_nn_z;

  // Field aligned of Div(vn_perp*nn) = -Div(Dn*Grad_perp(nn)) = -Dn*Lapla_perp(nn) - nn*Grad_per(Dn)*Grad_perp(log(nn))
  Field3D div_vnperpnn_aligned =  -toFieldAligned( Dn*LaplPerpnn + GradPerpDn_GradPerpnn , "RGN_NOBNDRY");

  // vn_perp aligned
  Field3D vnperpx = -Dn/nn*GradPerp_nn_x, vnperpz = -Dn/nn*GradPerp_nn_z;
  Field3D vnperpx_aligned = toFieldAligned(vnperpx, "RGN_NOBNDRY");
  Field3D vnperpz_aligned = toFieldAligned(vnperpz, "RGN_NOBNDRY");

  mesh->communicate(Dn_nn_Gradpar_vn_aligned, div_vnperpnn_aligned, vnperpx_aligned, vnperpz_aligned);
  Dn_nn_Gradpar_vn_aligned.applyBoundary("free_o3");
  div_vnperpnn_aligned.applyBoundary("free_o3");
  vnperpx_aligned.applyBoundary("free_o3");
  vnperpz_aligned.applyBoundary("free_o3");

  // nn * vperp * Grad_perp(vn) on staggered grid
  Field3D vnperpx_stag = myInterp_to(vnperpx_aligned, CELL_YLOW);
  Field3D vnperpz_stag = myInterp_to(vnperpz_aligned, CELL_YLOW);
  Field3D vnperp_GradPerpvn = coordinates_stag->g11*VDDX(vnperpx_stag, vn, CELL_YLOW, orderupwindscheme) +
                              coordinates_stag->g33*vnperpz_stag*DDZ(vn);
  if (yderivatives) {
    vnperp_GradPerpvn -= 
            coordinates_stag->g11*coordinates_stag->g_12/coordinates_stag->g_22*
            fromFieldAligned(VDDY(vnperpx_aligned, vn_aligned, CELL_YLOW, "U2"), "RGN_NOBNDRY")
          + coordinates_stag->g33*coordinates_stag->g_23/coordinates_stag->g_22*
            fromFieldAligned(VDDY(vnperpz_aligned, vn_aligned, CELL_YLOW, "U2"), "RGN_NOBNDRY");
  }
  Field3D nn_vperp_GradPerpvn = nn_stag*vnperp_GradPerpvn;

  neutral_density_equation["parallel_advection"] =
    - vn_centre*myGrad_par(lognn_aligned, CELL_CENTRE)
    - myDiv_par(vn_aligned, CELL_CENTRE);
  neutral_density_equation["perpendicular_diffusion"] =
    Dn*LaplPerpnn/nn + GradPerpDn_GradPerpnn/nn;
  neutral_density_equation["ionisation_recombination"] = S/nn;
  neutral_density_equation["neutral_source"] = Snn/nn;

  neutral_momentum_equation["parallel_advection"] =
    - vn*myGrad_par(nvn_aligned, CELL_YLOW)
    - nvn*myDiv_par(vn_aligned, CELL_YLOW);
  neutral_momentum_equation["density_gradient"] =
    - Tn_stag*myGrad_par(nn_aligned, CELL_YLOW);
  neutral_momentum_equation["perpendicular_diffusion"] =
    - vn*myInterp_to(div_vnperpnn_aligned, CELL_YLOW)
    - nn_vperp_GradPerpvn;
  neutral_momentum_equation["parallel_diffusion"] =
    myDiv_par(Dn_nn_Gradpar_vn_aligned, CELL_YLOW);
  neutral_momentum_equation["ion_friction"] = Fi;
  neutral_momentum_equation["electron_friction"] = Fe/mu;
  neutral_momentum_equation["neutral_momentum_source"] = Snvn;

  if (munvn > 0.) {
    neutral_momentum_equation["neutral_viscosity"] = munvn*Delp2(nvn);
  }

  if (mesh->hasBndryLowerY()) {
    int j = mesh->ystart;
    for (int i=mesh->xstart; i<=mesh->xend; ++i) {
      for (int k=0; k<mesh->LocalNz; ++k) {
        ddt(nvn)(i,j,k) = 0.;
      }
    }
  }
}

// Need to declare neutral_dvpar_instances static variable so that it can be linked
// correctly in timestepmonitor_func
std::vector<NeutralDVpar*> NeutralDVpar::neutral_dvpar_instances;

int NeutralDVpar::timestepmonitor_func(Solver*, BoutReal simtime, BoutReal) {
  for (auto neutral_dvpar : neutral_dvpar_instances) {
    int retval = neutral_dvpar->timestepmonitor(simtime);
    if (retval) {
      // If return value was not zero, there was an error, so return immediately
      return retval;
    }
  }
  return 0;
}

int NeutralDVpar::timestepmonitor(BoutReal simtime) {
  if (simtime - monitor_timelast >= dt_update) {
    monitor_timelast = simtime;
    updaterates = true;
  }
  else {
    updaterates = false;
  }

  return 0;
}

void NeutralDVpar::printMinMaxMean() {
  output.write("min(nn) = %e, max(nn) = %e, mean(nn) = %e, ",
    min(nn,true,"RGN_NOBNDRY"),
    max(nn,true,"RGN_NOBNDRY"),
    mean(nn,true,"RGN_NOBNDRY")
  );
  output.write("min(ddt(lognn)) = %e, max(ddt(lognn)) = %e, mean(ddt(lognn)) = %e\n",
    min(ddt(lognn),true,"RGN_NOBNDRY"),
    max(ddt(lognn),true,"RGN_NOBNDRY"),
    mean(ddt(lognn),true,"RGN_NOBNDRY")
  );
  output.write("min(nvn) = %e, max(nvn) = %e, mean(nvn) = %e, ",
    min(nvn,true,"RGN_NOBNDRY"),
    max(nvn,true,"RGN_NOBNDRY"),
    mean(nvn,true,"RGN_NOBNDRY")
  );
  output.write("min(ddt(nvn)) = %e, max(ddt(nvn)) = %e, mean(ddt(nvn)) = %e\n",
    min(ddt(nvn),true,"RGN_NOBNDRY"),
    max(ddt(nvn),true,"RGN_NOBNDRY"),
    mean(ddt(nvn),true,"RGN_NOBNDRY")
  );
  output.write("min(vn) = %e, max(vn) = %e, mean(vn) = %e\n\n",
    min(vn,true,"RGN_NOBNDRY"),
    max(vn,true,"RGN_NOBNDRY"),
    mean(vn,true,"RGN_NOBNDRY")
  );
}

void NeutralDVpar::precon(BoutReal UNUSED(t), BoutReal gamma, BoutReal UNUSED(delta)) {
  // Neutral gas diffusion
  // Solve (1 - gamma*Dn*Delp2)^{-1} 

  if (precon_firststep) {
    output << endl;
    output << "Creating laplacian for preconditioning neutral equations." << endl;
    precon_firststep = false;
    inv = Laplacian::create(Options::getRoot()->getSection("preconSolver"));

    // Zero value outer boundary   
    inv->setInnerBoundaryFlags(3);
    
    inv->setCoefA(1.0);
  }
  inv->setCoefD(-gamma*Dn);
  ddt(lognn) = inv->solve(ddt(lognn));
}

void NeutralDVpar::recycleFluxes(BoutReal time) {
  // Apply boundary conditions at the target plates:
  // bndry_neutrals = 1: vn = vth, nvn = -Rc*n*U, nn = nvn/vn
  // bndry_neutrals = 2: Neumann for nn and nvn, vn = nvn/nn, source of density and momentum in the first cell
  // bndry_neutrals = 3: Neumann for nn, vn = vth, nvn = nn*vn, source of density and momentum in the first cell
  // bndry_neutrals = 4: Neumann for nn, vn = vth, nvn = nn*vn, source of density and momentum exponentially decaying from targets
  // bndry_neutrals = 5: Neumann for nn and nvn, vn = nvn/nn, source of density and momentum exponentially decaying from targets
  BoutReal vnth = sqrt(Tn0/T0); 
  lognn_aligned.applyBoundary(time);
  nvn_aligned.applyBoundary(time);
  nn_lower = 0.;
  nn_upper = 0.;
  int i, j;

  if (bndry_neutrals == 1) {
    // nn = -Rc*n*U/vth
    j =  mesh->ystart;
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      i = r.ind;
      for(int k=0; k <mesh->LocalNz; ++k){
        nn_lower(i,k) = recycled_lower(i,k)/vnth;
        lognn_aligned(i,j-1,k) = (8.*log(nn_lower(i,k)) - 6.*lognn_aligned(i,j,k) + lognn_aligned(i,j+1,k))/3.;
        lognn_aligned(i,j-2,k) = 3.*lognn_aligned(i,j-1,k) - 3.*lognn_aligned(i,j,k) + lognn_aligned(i,j+1,k);
      }
    }

    j =  mesh->yend;
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      i = r.ind ;
      for(int k=0; k <mesh->LocalNz; ++k){
        nn_upper(i,k) = recycled_upper(i,k)/vnth;
        lognn_aligned(i,j+1,k) = (8.*log(nn_upper(i,k)) - 6.*lognn_aligned(i,j,k) + lognn_aligned(i,j-1,k))/3.;
        lognn_aligned(i,j+2,k) = 3.*lognn_aligned(i,j+1,k) - 3.*lognn_aligned(i,j,k) + lognn_aligned(i,j-1,k);
      }
    }
  }
  else {
    j =  mesh->ystart;
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      i = r.ind;
      for(int k=0; k <mesh->LocalNz; ++k){
        nn_lower(i,k) = exp(lognn_aligned(i,j,k));
      }
    }

    j =  mesh->yend;
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      i = r.ind ;
      for(int k=0; k <mesh->LocalNz; ++k){
        nn_upper(i,k) = exp(lognn_aligned(i,j,k));
      }
    }
  }

  // Boundary conditions for vn and nvn
  if (bndry_neutrals == 1 || bndry_neutrals == 3 || bndry_neutrals == 4) {
    // vn = vth
    vn_lower = vnth;
    vn_upper = -vnth;
    
    // nvn = nn*vn
    j =  mesh->ystart;
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      i = r.ind;
      for(int k=0; k <mesh->LocalNz; ++k){
        nvn_aligned(i,j,k) = vn_lower(i,k)*nn_lower(i,k);
        nvn_aligned(i,j-1,k) = 3.*nvn_aligned(i,j,k) - 3.*nvn_aligned(i,j+1,k) + nvn_aligned(i,j+2,k);
        nvn_aligned(i,j-2,k) = 3.*nvn_aligned(i,j-1,k) - 3.*nvn_aligned(i,j,k) + nvn_aligned(i,j+1,k);
      }
    }

    j =  mesh->yend;
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      i = r.ind ;
      for(int k=0; k <mesh->LocalNz; ++k){
        nvn_aligned(i,j+1,k) = vn_upper(i,k)*nn_upper(i,k);
        nvn_aligned(i,j+2,k) = 3.*nvn_aligned(i,j+1,k) - 3.*nvn_aligned(i,j,k) + nvn_aligned(i,j-1,k);
      }
    }
  }
  else {
    // nvn is Neumann, vn = nvn/nn
    j = mesh->ystart;
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      i = r.ind;
      for (int k = 0; k < mesh->LocalNz; ++k) {
        BoutReal vnl = nvn_aligned(i,j,k)/nn_lower(i,k);
        if (vnl > vnth) vnl = vnth;
        vn_lower(i,k) = vnl;
      }
    }

    j = mesh->yend;
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      i = r.ind;
      for (int k = 0; k < mesh->LocalNz; ++k) {
        BoutReal vnu = nvn_aligned(i,j+1,k)/nn_upper(i,k);
        if (vnu < -vnth) vnu = -vnth;
        vn_upper(i,k) = vnu;
      }
    }
  }
   
  if (bndry_neutrals == 2 || bndry_neutrals == 3) {
    // Sources to recycle density and momentum in the first cell 
    auto neutral_density_recycling =
      neutral_density_equation["recycling_source"].localAccessor();
    auto neutral_momentum_recycling =
      neutral_momentum_equation["recycling_source"].localAccessor();
    j = mesh->ystart;
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      i = r.ind;
      BoutReal dd = 1./( sqrt(coordinates_stag->g_22(i,j)) * coordinates_stag->dy(i,j) );
      for (int k = 0; k < mesh->LocalNz; k++) {
        BoutReal dnndt = recycled_lower(i,k)*dd;
        neutral_density_recycling(i,j,k) += dnndt/nn(i,j,k);
        neutral_momentum_recycling(i,j+1,k) += dnndt * vnth;
      }
    }

    j = mesh->yend;
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      i = r.ind;
      BoutReal dd = 1./( sqrt(coordinates_stag->g_22(i,j+1)) * coordinates_stag->dy(i,j+1) );
      for (int k = 0; k < mesh->LocalNz; k++) {
        BoutReal dnndt = recycled_upper(i,k)*dd;
        neutral_density_recycling(i,j,k) += dnndt/nn(i,j,k);
        neutral_momentum_recycling(i,j,k) -= dnndt * vnth;
      }
    }
  }
  
  if (bndry_neutrals == 4 || bndry_neutrals == 5) {
    // Recycle fluxes at the sheats by modifying the source terms (exponential decay)
    MPI_Comm comm;
    Snn = 0.;
    Snvn = 0.;
    int Nx_local = mesh->xend - mesh->xstart + 1;
    BoutReal* particlesdt = new BoutReal[Nx_local];
    j =  mesh->ystart;
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      i = r.ind;
      particlesdt[i-mesh->xstart] = 0.;
      for (int k=0; k <mesh->LocalNz; ++k)
        particlesdt[i-mesh->xstart] += recycled_lower(i,k);
      
      particlesdt[i-mesh->xstart] *= coordinates_stag->J(i,j)/( sqrt(coordinates_stag->g_22(i,j)) * BoutReal(mesh->LocalNz) );
    }
    for (int ix=mesh->xstart; ix<=mesh->xend; ++ix) {
      if (!mesh->periodicY(ix)) {
        comm = mesh->getYcomm(ix);
        MPI_Bcast(&particlesdt[ix-mesh->xstart],1,MPI_DOUBLE,0,comm);
        for (int iy=mesh->ystart; iy<=mesh->yend; ++iy) {
          Snn(ix,iy)  += particlesdt[ix-mesh->xstart]*profileSn_lower(ix,iy);
          Snvn(ix,iy) += vnth*particlesdt[ix-mesh->xstart]*profileSn_lower_stag(ix,iy);
        }
      }
    }

    j =  mesh->yend;
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      i = r.ind ;
      particlesdt[i-mesh->xstart] = 0.;
      for (int k=0; k <mesh->LocalNz; ++k)
        particlesdt[i-mesh->xstart] += recycled_upper(i,k);
      
      particlesdt[i-mesh->xstart] *= coordinates_stag->J(i,j+1)/( sqrt(coordinates_stag->g_22(i,j+1)) * BoutReal(mesh->LocalNz) );
    }
    for (int ix=mesh->xstart; ix<=mesh->xend; ++ix) {
      if (!mesh->periodicY(ix)) {
        comm = mesh->getYcomm(ix);
        int nproc;
        MPI_Comm_size(comm,&nproc);
        MPI_Bcast(&particlesdt[ix-mesh->xstart],1,MPI_DOUBLE,nproc-1,comm);
        for (int iy=mesh->ystart; iy<=mesh->yend; ++iy) {
          Snn(ix,iy)  += particlesdt[ix-mesh->xstart]*profileSn_upper(ix,iy);
          Snvn(ix,iy) -= vnth*particlesdt[ix-mesh->xstart]*profileSn_upper_stag(ix,iy);
        }
      }
    }

    delete [] particlesdt;
  }

}

// Compute the profile of the neutral sources, assuming an exponential decay in the poloidal plane along y from the targets
void NeutralDVpar::initialiseSource(Options &options) {
  OPTION(options, lambdaSnn, 0.1);
  profileSn_lower = 0.;
  profileSn_upper = 0.;
  profileSn_lower_stag = 0.;
  profileSn_upper_stag = 0.;
  profileSn_lower_stag.setLocation(CELL_YLOW);
  profileSn_upper_stag.setLocation(CELL_YLOW);
  Field2D hthe;
  GRID_LOAD(hthe);
  MPI_Comm comm;
  int Ny_local = mesh->yend - mesh->ystart + 1;
  BoutReal* sendbuffer = new BoutReal[Ny_local];
  for (int ix=mesh->xstart; ix<=mesh->xend; ++ix) {
    if (!mesh->periodicY(ix)) {
      comm = mesh->getYcomm(ix);
      int me,nproc;
      MPI_Comm_size(comm,&nproc);
      MPI_Comm_rank(comm,&me);
      int Ny = nproc*Ny_local;
      int globalystart = me*Ny_local;
      int globalyend = (me+1)*Ny_local-1;

      // Compute the distance along a flux surface in the poloidal plane (in m) from the target 
      BoutReal* recvbuffer = new BoutReal[Ny];
      for (int iy=mesh->ystart; iy<=mesh->yend; ++iy) {
        sendbuffer[iy-mesh->ystart] = coordinates_centre->dy(ix,iy)*hthe(ix,iy);
      }
      MPI_Allgather(sendbuffer,Ny_local,MPI_DOUBLE,recvbuffer,Ny,MPI_DOUBLE,comm);
      BoutReal* ydistance_lower  = new BoutReal[Ny];
      BoutReal* ydistance_upper  = new BoutReal[Ny];
      ydistance_lower[0] = 0.5*recvbuffer[0];
      ydistance_upper[Ny-1] = 0.5*recvbuffer[Ny-1];
      for (int iy=1; iy<Ny; ++iy) {
        ydistance_lower[iy] = ydistance_lower[iy-1] + 0.5*(recvbuffer[iy-1] + recvbuffer[iy]);
        ydistance_upper[Ny-1-iy] = ydistance_upper[Ny-iy] + 0.5*(recvbuffer[Ny-iy] + recvbuffer[Ny-1-iy]);
      }
      
      // Compute the two exponentials
      BoutReal* profile_lower = new BoutReal[Ny];
      BoutReal* profile_upper = new BoutReal[Ny];
      for (int iy=0; iy<Ny; ++iy) {
        profile_lower[iy] = exp(-ydistance_lower[iy]/lambdaSnn);
        profile_upper[iy] = exp(-ydistance_upper[iy]/lambdaSnn);
      }
      
      // Normalize the profiles to 1: sum_y(profile[y]*dy[y]*J[y]) = 1
      for (int iy=mesh->ystart; iy<=mesh->yend; ++iy) {
        sendbuffer[iy-mesh->ystart] = coordinates_centre->dy(ix,iy)*coordinates_centre->J(ix,iy);
      }
      MPI_Allgather(sendbuffer,Ny_local,MPI_DOUBLE,recvbuffer,Ny,MPI_DOUBLE,comm);
      BoutReal intlower = 0., intupper = 0.;
      for (int iy=0; iy<Ny; ++iy) {
        intlower += profile_lower[iy]*recvbuffer[iy];
        intupper += profile_upper[iy]*recvbuffer[iy];
      }

      for (int iy=0; iy<Ny; ++iy) {
        if(iy >= globalystart && iy <= globalyend) {
          int localyindex = iy - globalystart + mesh->ystart;
          
          profileSn_lower(ix,localyindex) = profile_lower[iy]/intlower;
          profileSn_upper(ix,localyindex) = profile_upper[iy]/intupper;
        }
      }
      delete [] profile_upper;
      delete [] profile_lower;
      delete [] ydistance_upper;
      delete [] ydistance_lower;
      delete [] recvbuffer;
    }
  }
  delete [] sendbuffer;
  mesh->communicate(profileSn_lower);
  mesh->communicate(profileSn_upper);
  profileSn_lower.applyBoundary("free_o3");
  profileSn_upper.applyBoundary("free_o3");
  for (int ix=mesh->xstart; ix<=mesh->xend; ++ix) {
    for( int iy=mesh->ystart; iy<=mesh->yend; ++iy) {
      profileSn_lower_stag(ix,iy) = 0.5*(profileSn_lower(ix,iy-1) + profileSn_lower(ix,iy));
      profileSn_upper_stag(ix,iy) = 0.5*(profileSn_upper(ix,iy-1) + profileSn_upper(ix,iy));
    }
  }
  mesh->communicate(profileSn_lower_stag);
  mesh->communicate(profileSn_upper_stag);
  profileSn_lower_stag.applyBoundary("free_o3");
  profileSn_upper_stag.applyBoundary("free_o3");
  SAVE_ONCE(profileSn_lower);
  SAVE_ONCE(profileSn_upper);
  SAVE_ONCE(profileSn_lower_stag);
  SAVE_ONCE(profileSn_upper_stag);
}

