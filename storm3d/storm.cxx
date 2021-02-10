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

#include <derivs.hxx>
#include <field_factory.hxx>
#include <invert_parderiv.hxx>
#include <bout/assert.hxx>
#include <bout/constants.hxx>

#include <vector>
#include <numeric>

class STORM;

BOUTMAIN(STORM);

int STORM::init(bool restarting) {
  
  globalOptions = Options::getRoot();
  options = Options::getRoot()->getSection("storm");
  coordinates_centre = mesh->getCoordinates(CELL_CENTRE);
  
  OPTION(options, mu_n0,                  -1.0) ;
  OPTION(options, mu_vort0,               -1.0) ;
  OPTION(options, mu,                     -1.0) ;
  OPTION(options, nu_parallel0,           -1.0) ;
  OPTION(options, g0,                     -1.0) ;  
  OPTION(options, kappa0,                 -1.0) ;
  OPTION(options, kappa0_perp,            -1.0) ;
  OPTION(options, diff_perp_U,             -1.) ;
  OPTION(options, diff_perp_V,             -1.) ;
  OPTION(options, beta0,                  -1.0) ;
  OPTION(options, phi_wall,                  0) ;
  OPTION(options, old_phi_wall_value,    false) ;
  
  OPTION(options, isothermal,            false) ;
  OPTION(options, boussinesq,             true) ;
  OPTION(options, electromagnetic,       false) ;
  OPTION(options, uniform_diss_paras,    false) ;
  OPTION(options, run_1d,                false) ;
  OPTION(options, hydrodynamic,          false) ;
  OPTION(options, symmetry_plane,         true) ;
  OPTION(options, normalise_lengths,      true) ;
  OPTION(options, add_blob,               true) ;
  OPTION(options, two_blobs,             false) ;
  OPTION(options, initial_noise,         false) ;
  OPTION(options, initial_perturbation,  false) ;
  OPTION(options, verbose,               false) ;
  OPTION(options, save_aligned_fields,   false) ;
  
  OPTION(options, B_0,                     0.5) ;   // Tesla 
  OPTION(options, T_e0,                     40) ;   // eV
  OPTION(options, T_i0,                     40) ;   // eV
  OPTION(options, m_i,                       2) ;   // Atomic Units
  OPTION(options, q,                         7) ;   // Dimensionless
  OPTION(options, R_c,                     1.5) ;   // m
  OPTION(options, n_0,                  0.8e19) ;   // m^-3
  OPTION(options, loglambda,                -1) ;   // Dimensionless
  
  OPTION(options, equilibrium_source, "1d_profiles"); // possible values: 1d_profiles, input_file, profiles_file
  OPTION(options, equilibrium_file_path, "");
  if (equilibrium_file_path == "") {
    // The user has not explicitly set equilibrium_file_path, use path relative to simulation output directory instead
    std::string equilibrium_directory;
    OPTION(options, equilibrium_directory, "equilibrium");
    std::string data_dir;
    Options::getRoot()->get("datadir", data_dir, "data");
    equilibrium_file_path = data_dir + "/" + equilibrium_directory;
  }
  OPTION(options, equilibrium_data_file, "");

  OPTION(options, average_radial_boundaries_core_SOL,            false) ;

  // Set default values for boundary conditions for all variables.
  // Can be overridden in input file, but should never need to be.
  // Needs to be called before adding anything to the time-solver or using setBoundary().
  setBoundaryConditionsOptions();

  // Check if we're using ShiftedMetric to run simulation in x-z orthogonal coordinates
  auto& mesh_options = (*globalOptions)["mesh"];
  if (mesh_options["paralleltransform"] == "shifted") {
    // As we use staggered grids, we must transform to/from field-aligned coordinates
    // rather than using 'parallel slices'. Therefore we don't want to calculate parallel
    // slices in 'mesh->communicate' as it would just be wasted computation, so set this
    // option must be set to false to disable it.
    if (mesh_options["calcParallelSlices_on_communicate"].as<bool>() == true) {
      throw BoutException("Option mesh:calcParallelSlices_on_communicate is set to true. "
          "We don't use parallel slices, so calculating them would be a waste of time. "
          "This option should be set to false.");
    }
  }

  // Load in ixseps
  mesh->get(ixseps1, "ixseps1") ; 
  mesh->get(ixseps2, "ixseps2") ; 
  mesh->get(jyseps1_1, "jyseps1_1");
  mesh->get(jyseps2_1, "jyseps2_1");
  mesh->get(jyseps1_2, "jyseps1_2");
  mesh->get(jyseps2_2, "jyseps2_2");
  mesh->get(ny_inner, "ny_inner");

  // Variables for evolving radial phi boundary conditions
  OPTION(options, evolving_bcs,          false);
  OPTION(options, time_update_SOL,         10.); // Length of time interval for outer boundary averaging (SOL) 
  OPTION(options, time_update_PF,          10.); // Length of time interval for inner boundary averaging (PF)
  OPTION(options, time_update_core,        10.); // Length of time interval for inner boundary averaging (core)
  if (evolving_bcs) {
    // check that jyseps* were set explicitly, if not set the default values (as defaults calculated in mesh initialization are not picked up here)
    if ( jyseps1_1 == 0 && jyseps2_1 == 0 && jyseps1_2 == 0 && jyseps2_2 == 0 ){
      int ny;
      mesh->get(ny, "ny");
      jyseps1_1 = -1;
      jyseps1_2 = ny / 2;
      jyseps2_1 = jyseps1_2;
      jyseps2_2 = ny - 1;
      ny_inner = jyseps2_1;
    }
  }

  //********* Parameters ***********
  
  // Specify Constants
  u = 1.66053892e-27 ;              // kg
  e = 1.602176565e-19 ;             // C
  m_e = 9.10938291e-31;             // Kg
  epsilon_0 = 8.854187817e-12 ; 
  mu_0 = 4*PI*1e-7 ;                
  
  // Convert Parameters to SI units 
  T_e0 = e*T_e0 ;                   // Joules
  T_i0 = e*T_i0 ;                   // Joules
  
  m_i = m_i*u ;                     // kg 
  
  // Calculate secondary parameters
  c_s = sqrt(T_e0/m_i) ;
  Omega_i = e*B_0/m_i ; 
  rho_s = c_s/Omega_i ;  
  
  V_the = sqrt(T_e0/m_e) ; 
  V_thi = sqrt(T_i0/m_i) ; 
  
  Omega_e = e*B_0/m_e ; 
  rho_e = V_the/Omega_e ; 
  rho_i = V_thi/Omega_i ; 
  
  if (loglambda < 0){
    loglambda = 18.0 - log(sqrt(n_0/1.0e19)*pow(T_e0/(1000.0*e), -1.5)) ; 
  } 
  
  nu_ei = n_0*pow(e, 4)*loglambda/(sqrt(m_e)*SQ(epsilon_0)*3.0*pow(TWOPI*T_e0, 1.5)) ; // electron - ion collision frequency
  nu_ii = n_0*pow(e, 4)*loglambda/(sqrt(m_i)*SQ(epsilon_0)*3.0*pow(TWOPI*T_i0, 1.5)*sqrt(2)) ; // ion - ion collision frequency 
  
  // Simulation parameters
  if (mu < 0){
    mu = m_i/m_e ;
  } 
  if (nu_parallel0 < 0 ){
    nu_parallel0 = 0.51*nu_ei/Omega_i ; 
  }
  if(mu_n0 < 0){
    mu_n0 = (1.0 + 1.3*pow(q, 2))*(1.0 + T_i0/T_e0)*pow(rho_e, 2)*nu_ei ; // Neo-Classical Particle diffusion, m^2s-1
    mu_n0 = mu_n0/(SQ(rho_s)*Omega_i) ; 
  }
  if(mu_vort0 < 0){
    mu_vort0 = (1.0 + 1.6*pow(q, 2))*(6.0/8.0)*pow(rho_i, 2)*nu_ii ;      // Neo-Classical Ion viscosity, m^2s-1
    mu_vort0 = mu_vort0/(SQ(rho_s)*Omega_i) ; 
  }
  if(g0 < 0){
    g0 = 2.0*rho_s/R_c ; 
  }
  
  if(kappa0 < 0 ){
    kappa0 = 3.16*(SQ(V_the)/nu_ei)/(SQ(rho_s)*Omega_i);
  }

  if(kappa0_perp < 0 ){
    kappa0_perp = (1.0 + 1.6*SQ(q))*(4.66)*SQ(rho_e)*nu_ei/(SQ(rho_s)*Omega_i);
  }

  if(beta0 < 0) {
    // ratio of nominal plasma pressure to nominal magnetic pressure
    beta0 = 2. * mu_0 * n_0 * T_e0 / SQ(B_0) ;
  }

  Vsheath_BC_prefactor = sqrt((mu/TWOPI)/(1.0+(1.0/mu))) ;
  Usheath_BC_prefactor = sqrt(1.0/(1.0+(1.0/mu))) ;
  
  // Get from options or calculate Lx, Ly and Lz
  set_Lx_Ly_Lz();

  //********* Metric, Bracket scheme and staggering of fields ***********

  // Normalise grid mesh.   
  if (normalise_lengths){
    coordinates_centre->g_11 /= SQ(rho_s) ;
    coordinates_centre->g_22 /= SQ(rho_s) ;
    coordinates_centre->g_33 /= SQ(rho_s) ;
    coordinates_centre->g_12 /= SQ(rho_s) ;
    coordinates_centre->g_13 /= SQ(rho_s) ;
    coordinates_centre->g_23 /= SQ(rho_s) ;
    coordinates_centre->g11 *= SQ(rho_s) ;
    coordinates_centre->g22 *= SQ(rho_s) ;
    coordinates_centre->g33 *= SQ(rho_s) ;
    coordinates_centre->g12 *= SQ(rho_s) ;
    coordinates_centre->g13 *= SQ(rho_s) ;
    coordinates_centre->g23 *= SQ(rho_s) ;
    //mesh->dx /= rho_s ;
    //mesh->dz /= rho_s ;
    //mesh->dy /= rho_s ;
    //mesh->zlength/= rho_s ;
    coordinates_centre->geometry();
  }
  int bracket; 
  OPTION(options, bracket, 2);
  
  switch(bracket) {
    case 0: {
      bm = BRACKET_STD;
      output << "\tBrackets: default differencing\n";
      break;
    }
    case 1: {
      bm = BRACKET_SIMPLE; 
      output << "\tBrackets: simplified operator\n";
      break;
    }
    case 2: {
      bm = BRACKET_ARAKAWA; 
      output << "\tBrackets: Arakawa scheme\n";
      break;
    }
    case 3: {
      bm = BRACKET_CTU; 
      output << "\tBrackets: Corner Transport Upwind method\n";
      break;
    }
    default:
      output << "ERROR: Invalid choice of bracket method. Must be 0-3\n";
      return 1;
  }
  
  // Stagger Fields
  V.setLocation(CELL_YLOW) ;
  U.setLocation(CELL_YLOW) ;
  chiU.setLocation(CELL_YLOW);
  chiV.setLocation(CELL_YLOW);
  psi.setLocation(CELL_YLOW);
  logn_stag.setLocation(CELL_YLOW) ;
  n_stag.setLocation(CELL_YLOW) ;
  S_stag.setLocation(CELL_YLOW) ;
  phi_stag.setLocation(CELL_YLOW) ;  
  nu_parallel.setLocation(CELL_YLOW) ;
  qpar_aligned.setLocation(CELL_YLOW);
  logT_stag.setLocation(CELL_YLOW);
  T_stag.setLocation(CELL_YLOW);
  UmV.setLocation(CELL_YLOW);
  Grad_par_phi_stag.setLocation(CELL_YLOW);
  Grad_par_T_stag.setLocation(CELL_YLOW);
  powT_1_5_stag.setLocation(CELL_YLOW);

  // ******** Initialise constant fields ********
  
  S = FieldFactory::get()->create3D("S:function", Options::getRoot(), mesh, CELL_CENTRE, 0);
  S_stag = FieldFactory::get()->create3D("S:function", Options::getRoot(), mesh, CELL_YLOW, 0);
  if (!isothermal) {
    S_E = FieldFactory::get()->create3D("S_E:function", Options::getRoot(), mesh, CELL_CENTRE, 0);
  }
 
  // Normalise S and S_stag by Ly, so that we do not need Ly in the input file
  S /= Ly;
  S_stag /= Ly;
  if (!isothermal) {
    S_E /= Ly;
  }

  if(normalise_lengths){
    S = S*rho_s ; 
    S_stag = S_stag*rho_s ; 
    if (!isothermal) {
      S_E = S_E*rho_s ;
    }
  }
    
  //************* Definition of the fileds to solve and Laplacian solver *************
  
  SOLVE_FOR2(logn, vort);
  n = exp(logn);
  SAVE_REPEAT(n);
  SAVE_REPEAT(phi);
  SAVE_ONCE(S);
  if (isothermal) {
    T = 1.;
    logT_aligned = toFieldAligned(logT, "RGN_NOX");
    logT_aligned = 0.;
    T_aligned = exp(logT_aligned, "RGN_NOX");
    logT_stag = stormInterp_to(logT_aligned, CELL_YLOW);
    logT_stag = 0.;
    T_stag = exp(logT_stag, "RGN_NOY");
    sqrtT = 1.;
    T_sheath_upper = exp(extrap_sheath_upper(logT_aligned));
    if (!symmetry_plane) {
      T_sheath_lower = exp(extrap_sheath_lower(logT_aligned));
    }
    SAVE_ONCE(T);
  } else {
    SOLVE_FOR(logT);
    T = exp(logT);
    SAVE_REPEAT(T);
    SAVE_ONCE(S_E);
  }
  if (electromagnetic) {
    SOLVE_FOR2(chiU, chiV);
    SAVE_REPEAT3(U, V, psi);
  } else {
    solver->add(chiU, "U"); // chiU equivalent to U when psi is neglected
    solver->add(chiV, "V"); // chiV equivalent to V when psi is neglected
    U = chiU;
    V = chiV;
  }

  if (save_aligned_fields) {
    SAVE_REPEAT5(n_aligned, U_aligned, V_aligned, vort_aligned, phi_aligned);
    if (!isothermal) {
      SAVE_REPEAT2(T_aligned, qpar_aligned);
    }
    if (electromagnetic) {
      SAVE_REPEAT(psi_aligned);
    }
  }

  Options* phiOptions;
  if (globalOptions->isSection("phiSolver")) {
    phiOptions = globalOptions->getSection("phiSolver");
  } else {
    phiOptions = globalOptions->getSection("laplace");
  }
  phiSolver = Laplacian::create(phiOptions);
  phi = phi_wall ;  // Sets phi field to 0 for inital 'guess' of Laplace solver
  phi.setBoundary("phi") ;
  phi_aligned.setBoundary("phi_aligned") ;
  phi_stag.setBoundary("phi_stag") ;
  if (!boussinesq) {
    // Check that the Laplacian solver uses 3D coefficients
    ASSERT0(phiSolver->uses3DCoefs());
  }

  if (electromagnetic) {
    auto psi_options = globalOptions->getSection("psiSolver");
    psiSolver = Laplacian::create(psi_options, CELL_YLOW);

    // Check that the Laplacian solver uses 3D coefficients
    ASSERT0(psiSolver->uses3DCoefs());

    psiSolver->setCoefD(-1.0);
    psi = 0.; // make sure all guard cells are set to something
    psi.setBoundary("psi");
    psi_aligned.setBoundary("psi_aligned");
    psi_centre.setBoundary("psi_centre");

    use_psi_boundary_solver = (*psi_options)["use_psi_boundary_solver"].withDefault(true);
    if (use_psi_boundary_solver) {
      // Check there isn't a limiter-type boundary where we can't use a solver on the
      // boundary point
      // Must be: - either a divertor topology with non-zero length legs (so includes PFR
      //            region
      //          - or SOL only simulation
      //          - or core only simulation
      ASSERT0( (jyseps1_1 >= 0 and jyseps2_2 <= mesh->GlobalNy - 2*mesh->ystart)
            or (ixseps1 <= 0 and ixseps2 <= 0)
            or (ixseps1 >= mesh->GlobalNx and ixseps2 >= mesh->GlobalNx));

      std::string psi_solver_type = (*globalOptions)["psiSolver"]["type"];
      if (psi_solver_type == "naulin") {
        // Actually only need the FFT sub-solver as at the boundary we solve Delp2(psi) =
        // n*(U-V) which has constant coefficients
        psiBoundarySolver = Laplacian::create(psi_options->getSection("delp2solver"), CELL_YLOW);

        // Set boundary flags from psiSolver options, like LaplaceNaulin does in its
        // constructor
        psiBoundarySolver->setInnerBoundaryFlags(
            (*psi_options)["inner_boundary_flags"].as<int>());
        psiBoundarySolver->setOuterBoundaryFlags(
            (*psi_options)["outer_boundary_flags"].as<int>());
      } else {
        psiBoundarySolver = Laplacian::create(psi_options, CELL_YLOW);
      }

      if (mesh->hasBndryLowerY() and not symmetry_plane) {
        psi_sheath_lower = 0.;
        psi_sheath_lower.setLocation(CELL_YLOW).setIndex(mesh->ystart);
      }
      if (mesh->hasBndryUpperY()) {
        psi_sheath_upper = 0.;
        psi_sheath_upper.setLocation(CELL_YLOW).setIndex(mesh->yend+1);
      }
    }
  }
 
  U.setBoundary("U");
  U_aligned.setBoundary("U_aligned");
  U_centre.setBoundary("U"); // Use same boundary conditions for U_centre as for U (NB don't need parallel boundary conditions)
  V.setBoundary("V");
  V_aligned.setBoundary("V_aligned");
  V_centre.setBoundary("V"); // Use same boundary conditions for V_centre as for V (NB don't need parallel boundary conditions)
  check_U_V_x_boundary_conditions();
  UmV_centre.setBoundary("V"); // Use same boundary conditions for UmV as for V (NB don't need parallel boundary conditions)

  logn_aligned.setBoundary("n_aligned");
  n_aligned.setBoundary("n_aligned");
  logn_stag.setBoundary("n");
  n_stag.setBoundary("n");
  vort_aligned.setBoundary("vort_aligned");
  if (!isothermal) {
    logT_aligned.setBoundary("T_aligned");
    T_aligned.setBoundary("T_aligned");
    logT_stag.setBoundary("T");
    T_stag.setBoundary("T");
  }

  comms.add(logn);
  comms.add(vort) ;

  if (!isothermal) {
    qpar_aligned = 0.0;
    qpar_aligned.setBoundary("qpar_aligned");
    if (electromagnetic) {
      qpar_centre.setBoundary("qpar_centre");
    }
    comms.add(logT) ;
  }
  if(!boussinesq){
    uE2 = 0.0;
    uE2.setBoundary("uE2");
  }
  
  //**************** Initialisation of the fields *******************
  
  //Initialise backgrounds (if necessary)
  initialise_background() ;
  
  if (((mesh->firstX() && phi_array_inner == NULL)
        || (mesh->lastX() && phi_array_outer == NULL))
      && !evolving_bcs) {
    throw BoutException("Error: phi_array_inner or phi_array_outer has not been set.");
  }else if(evolving_bcs){
    //Intialise the variables for the averaged Neumann boundary conditions
    phi_bc_initialise(restarting);
  }
    
  //*************** Filament initialisation, if not in background ***************

  if (!restarting) {
    if (add_blob){
      initialise_blob("blob");
      // if two filaments are initialised, add the section blob_1 to BOUT.inp
      if (two_blobs){
        initialise_blob("blob_1");
      }
    }
    if (initial_noise){
      add_noise();
    }
    if (initial_perturbation) {
      add_perturbation(n, "n_pert");
      add_perturbation(U, "U_pert");
      add_perturbation(V, "V_pert");
      add_perturbation(phi, "phi_pert");
      if (electromagnetic) {
        add_perturbation(psi, "psi_pert");
      }
      vort = vort_from_phi(phi);
      if (!isothermal)
        add_perturbation(T, "T_pert");
    }
  }

  // set initial values of logn and logT
  logn = log(n);
  logT = log(T);

  // set initial values of chiU, chiV
  if (!restarting) {
    if (electromagnetic) {
      chiU = U + (beta0/2.0) * psi;
      chiV = V - mu * (beta0/2.0) * psi;
    } else {
      chiU = U;
      chiV = V;
    }
  }

  // Add phi and psi to restart files so we save the initial guess for
  // iterative solvers
  restart.addOnce(phi, "phi");
  if (electromagnetic) {
    restart.addOnce(psi, "psi");
    if (use_psi_boundary_solver) {
      SAVE_REPEAT(psi_sheath_upper);
      restart.addOnce(psi_sheath_upper, "psi_sheath_upper");
      if (not symmetry_plane) {
        SAVE_REPEAT(psi_sheath_lower);
        restart.addOnce(psi_sheath_lower, "psi_sheath_lower");
      }
    }
  }

  if (run_1d) {
    OPTION(options, run_1d_T_slowdown, 20.) ;
    hydrodynamic = true; // run_1d uses all the 'hydrodynamic' options, plus run_1d_T_slowdown
    if (mesh->xend > mesh->xstart) {
      throw BoutException("run_1d is true, but number of x-grid-points is >1");
    }
    if (mesh->LocalNz > 1) {
      throw BoutException("run_1d is true, but number of z-grid-points is >1");
    }
  }
   
  //************** Output of Summary **********************
  
  output.write("\n*******************************************************************") ; 
  output.write("\nGit Version of STORM: %s", STORM_VERSION) ;
  output.write("\n*******************************************************************") ; 
  output.write("\nCalculated Parameters") ; 
  output.write("\n*******************************************************************") ; 
  output.write("\n\tT_e0         = %e J", T_e0) ; 
  output.write("\n\tm_i          = %e kg", m_i) ; 
  output.write("\n\tc_s          = %e m/s", c_s) ; 
  output.write("\n\trho_s        = %e m", rho_s) ; 
  output.write("\n\tOmega_i      = %e s^-1", Omega_i) ; 
  output.write("\n\tOmega_e      = %e s-1", Omega_e) ; 
  output.write("\n\tnu_ei        = %e Hz", nu_ei) ; 
  output.write("\n\tnu_ii        = %e Hz", nu_ii) ;
  output.write("\n\tV_the        = %e m/s", V_the) ; 
  output.write("\n\tV_thi        = %e m/s", V_thi) ;
  output.write("\n\trho_e        = %e m", rho_e) ; 
  output.write("\n\trho_i        = %e m", rho_i) ; 
  output.write("\n\tloglambda    = %e ", loglambda) ; 
  output.write("\n\tbeta0        = %e ", beta0) ;
  output.write("\nDimensionless Parameters:") ;  
  output.write("\n\tmu_n0        = %e ", mu_n0) ; 
  output.write("\n\tmu_vort0     = %e ", mu_vort0) ; 
  output.write("\n\tdiff_perp_U  = %e ", diff_perp_U);
  output.write("\n\tdiff_perp_V  = %e ", diff_perp_V);
  output.write("\n\tnu_parallel0 = %e ", nu_parallel0) ; 
  output.write("\n\tg0           = %e ", g0) ; 
  output.write("\n\tmu           = %e ", mu) ;
  output.write("\n\tkappa_par    = %e ", kappa0);
  output.write("\n\tkappa_perp   = %e ", kappa0_perp);
  output.write("\nDimensionless Lengths:") ;  
  output.write("\n\tdx           = %e ", coordinates_centre->dx(mesh->xstart,mesh->ystart)) ;
  output.write("\n\tdy           = %e ", coordinates_centre->dy(mesh->xstart,mesh->ystart)) ;
  output.write("\n\tdz           = %e ", coordinates_centre->dz) ;
  output.write("\n\tzlength      = %e ", coordinates_centre->zlength()) ;
  output.write("\n\tLx           = %e ", Lx) ; 
  output.write("\n\tLy           = %e ", Ly) ; 
  output.write("\n\tLz           = %e ", Lz) ; 
  output.write("\nDimensional Lengths:") ;  
  output.write("\n\tLx           = %e m", Lx*rho_s) ; 
  output.write("\n\tLy           = %e m", Ly*rho_s) ; 
  output.write("\n\tLz           = %e m", Lz*rho_s) ; 

  if (run_1d) {
    output.write("\n\tWarning, run_1d is assumed to be finding an equilibrium"
        "\n\tEvolution of T is artificially slowed down to allow larger timesteps"
        "\n\tTime evolution with run_1d is therefore not physical, unless run_1d_T_slowdown is 1."
        "\n\trun_1d_T_slowdown   = %e", run_1d_T_slowdown);
  }
  output.write("\n*******************************************************************\n") ; 

  // Set up output for synthetic Langmuir probe trace
  if (fast_output.enabled) {

    // Add monitor if necessary
    if (fast_output.enable_monitor) {
      solver->addMonitor(&fast_output);
    }

    // Add points from the input file
    int i = 0;
    BoutReal xpos, ypos, zpos;
    int ix, iy, iz;
    Options* fast_output_options = globalOptions->getSection("fast_output");
    while (true) {
      // Add more points if explicitly set in input file
      fast_output_options->get("xpos"+std::to_string(i), xpos, -1.);
      fast_output_options->get("ypos"+std::to_string(i), ypos, -1.);
      fast_output_options->get("zpos"+std::to_string(i), zpos, -1.);
      if (xpos<0. || ypos<0. || zpos<0.) {
        output.write("\tAdded %i fast_output points\n", i);
        break;
      }
      ix = int(xpos*mesh->GlobalNx);
      iy = int(ypos*mesh->GlobalNy);
      iz = int(zpos*mesh->GlobalNz);
      fast_output.add("n"+std::to_string(i), n, ix, iy, iz);
      fast_output.add("T"+std::to_string(i), T, ix, iy, iz);
      i++;
    }
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////

int STORM::rhs(BoutReal time) { 

  if(verbose){ output<<"\r\t\t\t\t\t\t\t t = "<<time<<std::flush;} //Stream simulation time to screen

  mesh->communicate(comms) ;
  if (average_radial_boundaries_core_SOL) {
    average_Z_bndry(logn);
    average_YZ_bndry(logn, true, false, true);
    average_Z_bndry(vort);
    if (!isothermal) {
      average_Z_bndry(logT);
      average_YZ_bndry(logT, true, false, true);
    }
  }

  n = exp(logn, "RGN_NOY");
  if (!isothermal) T = exp(logT, "RGN_NOY");

  // Calculate field-aligned versions of fields
  logn_aligned = toFieldAligned(logn, "RGN_NOBNDRY");
  vort_aligned = toFieldAligned(vort, "RGN_NOBNDRY");
  if (!isothermal) {
    logT_aligned = toFieldAligned(logT, "RGN_NOBNDRY");
    mesh->communicate(logn_aligned, vort_aligned, logT_aligned);
    logT_aligned.applyBoundary(time);
    T_aligned = exp(logT_aligned, "RGN_NOX");
  }
  else {
    mesh->communicate(logn_aligned, vort_aligned);
  }
  logn_aligned.applyBoundary(time);
  vort_aligned.applyBoundary(time);
  n_aligned = exp(logn_aligned, "RGN_NOX");

  //*********** Calculate staggered interpolation for n and T, evaluate parameters depending on n and T ************
  logn_stag = stormInterp_to(logn_aligned, CELL_YLOW);
  if (not isothermal) {
    logT_stag = stormInterp_to(logT_aligned, CELL_YLOW);
    mesh->communicate(logn_stag, logT_stag);
    logT_stag.applyBoundary(time);
    if (average_radial_boundaries_core_SOL) {
      average_Z_bndry(logT_stag);
      average_YZ_bndry(logT_stag, true, false, true);
    }
    T_stag = exp(logT_stag, "RGN_NOY");
  }
  else {
    mesh->communicate(logn_stag);
  }
  logn_stag.applyBoundary(time);
  if (average_radial_boundaries_core_SOL) {
    average_Z_bndry(logn_stag);
    average_YZ_bndry(logn_stag, true, false, true);
  }
  n_stag = exp(logn_stag, "RGN_NOY");

  if(!isothermal){
    p = n*T;
    sqrtT = sqrt(T, "RGN_NOY"); // Need x boundary guard cells, because we will take x derivatives of sqrtT
  }else{
    p = n;
  }
  
  //********** Extrapolate BC for density and temperature to calculate BC on V, U and q_par **********
  // Extrapolate in field-aligned coordinates since these are used for parallel boundary
  // conditions
  n_sheath_upper = exp(extrap_sheath_upper(logn_aligned));
  if(!symmetry_plane){
    n_sheath_lower = exp(extrap_sheath_lower(logn_aligned));
  }
  if (!isothermal) {
    T_sheath_upper = exp(extrap_sheath_upper(logT_aligned));
    if(!symmetry_plane){
      T_sheath_lower = exp(extrap_sheath_lower(logT_aligned));
    }
  }

  //********* Laplace inversion: calculation of phi, including value at sheath *********

  if (!boussinesq){
    phiSolver->setCoefD(n);
    phiSolver->setCoefC2(n);
  } 

  if (hydrodynamic) {
    // in 1d or 2d this does nothing, but don't include ExB drifts if
    // hydrodynamic mode is used in 3d
    phi = 0.;
  }
  else{
    if(!evolving_bcs){
      set_xguards(phi, phi_array_inner, phi_array_outer);
    }else{
      apply_bndry_phi();
    }
    phi = phiSolver->solve(vort, phi) ;
  }
	
  phi_aligned = toFieldAligned(phi, "RGN_NOBNDRY");
  mesh->communicate(phi, phi_aligned);
  phi.applyBoundary(time);
  phi_aligned.applyBoundary(time);
  
  phi_stag = stormInterp_to(phi_aligned, CELL_YLOW);
  mesh->communicate(phi_stag);
  phi_stag.applyBoundary(time);
  // Note phi_stag is at CELL_YLOW, but we never take y-derivatives or y-interpolate it,
  // so don't ever need to communicate again to communicate ystart values set by boundary
  // condition

  // define uE2, set boundaries and communicate it (check staggered buisiness)
  if(!boussinesq){
    Vector3D grad_perp_phi = Grad_perp(phi);
    uE2 = grad_perp_phi*grad_perp_phi; // Note, cannot use SQ() here because SQ() is a template that returns the same type as its argument, whereas (Vector3D)*(Vector3D) is a dot product that returns Field3D
    mesh->communicate(uE2);
    uE2.applyBoundary(time);
  }  

  // Extrapolate phi at boundary for V boundary condition
  phi_sheath_lower = extrap_sheath_lower(phi_aligned) ;
  phi_sheath_upper = extrap_sheath_upper(phi_aligned) ;

  if (electromagnetic) {
    //********* Laplace inversion: calculation of psi, U and V *********
    // BC:
    //chi variables do not need explicit BCs as no gradients taken
    //instead, separately apply BCs to psi and U/V

    //***********PSI FIRST THEN U/V
    //***********Setup laplacian inversion to solve:
    //- delp2 psi + n * (mu + 1) * beta0/2 * psi = n * (chiU - chiV)
    //***********this originates from:
    //J = - Delp2 psi
    //J = n(U-V)
    //U = chiU - beta0/2 psi
    //V = chiV + mu beta0/2 psi

    //U = chiU - beta0/2 psi
    psiSolver->setCoefA((beta0/2.0)*(mu+1.0)*n_stag);

    //solve for A_par, U and V
    //implicit interp_to(n,CELL_YLOW)
    //need to pass psi in as parameter for multigrid - initial guess for each timestep
    psi = psiSolver->solve((chiU - chiV)*n_stag, psi);
    psi_aligned = toFieldAligned(psi, "RGN_NOBNDRY");
    // psi_aligned boundary conditions handled below in applyPsiBoundaries()

    U = chiU - (beta0/2.0) * psi;
    V = chiV + mu * (beta0/2.0) * psi;

    U_aligned = toFieldAligned(U, "RGN_NOBNDRY");
    V_aligned = toFieldAligned(V, "RGN_NOBNDRY");

    mesh->communicate(psi, psi_aligned, U, V, U_aligned, V_aligned);
    U.applyBoundary(time);
    V.applyBoundary(time);
  } else {
    U = chiU;
    V = chiV;
    U_aligned = toFieldAligned(U, "RGN_NOBNDRY");
    V_aligned = toFieldAligned(V, "RGN_NOBNDRY");
    mesh->communicate(U, V, U_aligned, V_aligned);
    // As chiU and chiV were added to the solver as U and V, the (non-sheath) boundary
    // conditions for U and V have already been applied to them
  }
  if(average_radial_boundaries_core_SOL) {
    average_Z_bndry(U);
    average_Z_bndry(V);
  }
  U_aligned.applyBoundary(time);
  V_aligned.applyBoundary(time);
  
  //********** Sheath BC on the ion and electron velocity, U and V ********
  if (!symmetry_plane) {
    Vsheath_ydown_staggered(V_aligned, phi_sheath_lower, phi_wall, T_sheath_lower, Vsheath_BC_prefactor) ;
    Usheath_ydown_staggered(U_aligned, sqrt(T_sheath_lower), Usheath_BC_prefactor);
  }
  Vsheath_yup_staggered(V_aligned, phi_sheath_upper, phi_wall, T_sheath_upper, Vsheath_BC_prefactor);
  Usheath_yup_staggered(U_aligned, sqrt(T_sheath_upper), Usheath_BC_prefactor);
  if (mesh->yend - mesh->ystart + 1 < 3) {
    // Need to communicate mesh->ystart point that has just been set by boundary conditions
    mesh->communicate(U_aligned, V_aligned);
  }

  if (hydrodynamic) {
    V = U; // Set them equal even though we already set ddt(V)=ddt(U) so that the boundary values are exactly equal
    chiV = chiU; // These are the values that will be written to the output files, so set them equal too
    V_aligned = U_aligned;
  }

  if (not isothermal or electromagnetic or hydrodynamic) {
    if (not symmetry_plane) {
      V_sheath_lower = sliceXZ(V_aligned, mesh->ystart);
      U_sheath_lower = sliceXZ(U_aligned, mesh->ystart);
    }
    V_sheath_upper = sliceXZ(V_aligned, mesh->yend+1);
    U_sheath_upper = sliceXZ(U_aligned, mesh->yend+1);
  }

  if (electromagnetic) {
    // Now we have boundary conditions set on U and V, can use the current at the
    // sheath edge to set a boundary condition on psi
    applyPsiBoundaries();

    // Interpolate to psi_centre after applying boundary conditions
    psi_centre = stormInterp_to(psi_aligned, CELL_CENTRE);
    mesh->communicate(psi_centre);
    psi_centre.applyBoundary(time);
  }

  U_centre = stormInterp_to(U_aligned, CELL_CENTRE);
  V_centre = stormInterp_to(V_aligned, CELL_CENTRE);
  if (electromagnetic) {
    mesh->communicate(V_centre, U_centre);
    U_centre.applyBoundary(time);
    V_centre.applyBoundary(time);
    if(average_radial_boundaries_core_SOL) {
      average_Z_bndry(U_centre);
      average_Z_bndry(V_centre);
    }
  }
  UmV = U-V;
  UmV_centre = U_centre - V_centre;
  Field3D UmV_aligned = U_aligned - V_aligned;

  // Calculate coefficients which are dependent on density
  powT_1_5_stag = pow(T_stag, 1.5, "RGN_NOBNDRY");
  nu_parallel = n_stag*nu_parallel0/powT_1_5_stag ;
  
  if(!uniform_diss_paras){
    n_on_sqrt_T = n/sqrtT ;
    mu_n = mu_n0*n_on_sqrt_T;
    mu_vort = mu_vort0*n_on_sqrt_T;
    kappa_perp = n*kappa0_perp*n_on_sqrt_T;
  }

  // Precompute the curvature terms
  Curv_n = Curv(n);
  Curv_phi = Curv(phi);
  if(isothermal){
    Curv_p = Curv_n;
  }else{
    Curv_T = Curv(T);
    Curv_p = T*Curv_n + n*Curv_T; // According to Ben has a better stability than Curv(p)
  }
  
  //****************** Equations of the model **************************
  
  // ***Vorticity Equation***
  if (hydrodynamic) {
    ddt(vort) = 0.;
  }
  else {
    if(boussinesq){
      ddt(vort) = - bracket(phi, vort, bm, CELL_CENTRE)
                  + Div_par_EM(UmV_aligned, UmV_centre, psi_centre, CELL_CENTRE)
                  + UmV_centre*Grad_par_EM(logn_aligned, logn, psi_centre, CELL_CENTRE)
                  - Vpar_Grad_par_EM(U_aligned, vort_aligned, U_centre, vort,
                      psi_centre, CELL_CENTRE);
    }
    else{
      ddt(vort) = - bracket(phi, vort, bm, CELL_CENTRE)
                  + Div_par_EM(UmV_aligned, UmV_centre, psi_centre, CELL_CENTRE)*n
                  + UmV_centre*Grad_par_EM(n_aligned, n, psi_centre, CELL_CENTRE)
                  - Vpar_Grad_par_EM(U_aligned, vort_aligned, U_centre, vort,
                      psi_centre, CELL_CENTRE)
                  - bracket(uE2/2.,n,bm); //check normalisation of last term
    }

    if(uniform_diss_paras){
      ddt(vort) += mu_vort0*Delp2(vort) ;
    }
    else{
     ddt(vort) +=  mu_vort*Delp2(vort) + Grad_perp_dot_Grad_perp(mu_vort,vort) ;
    }

    // Curvature terms for vorticity
    if (boussinesq){
      ddt(vort) += Curv_p/n;
    }
    else {
      ddt(vort) += Curv_p;
    }
  }
    
  Field3D div_par_V = Div_par_EM(V_aligned, V_centre, psi_centre, CELL_CENTRE);
  // ***Density Equation***
  ddt(logn) = - bracket(phi, logn, bm, CELL_CENTRE)
           - div_par_V
           - Vpar_Grad_par_EM(V_aligned, logn_aligned, V_centre, logn, psi_centre, CELL_CENTRE)
           + S/n;
  
  if(uniform_diss_paras){
    ddt(logn) += mu_n0*Delp2(n)/n;
  }
  else{
    ddt(logn) += mu_n*Delp2(n)/n + Grad_perp_dot_Grad_perp(mu_n, logn) ;
  }

  // Curvature terms for density 
  ddt(logn) += -Curv_phi + Curv_p/n;
  
  // ***Ion parallel velocity Equation*** 
  Grad_par_T_stag = Grad_par_EM(T_aligned, T_stag, psi, CELL_YLOW);
  if (hydrodynamic) {
    Grad_par_phi_stag = mu/(mu+1.)*(Grad_par_T_stag
                        + T_stag*fromFieldAligned(Grad_par(logn_aligned, CELL_YLOW),
                                                  "RGN_NOBNDRY"))
                        + 0.71*Grad_par_T_stag;
  } else {
    Grad_par_phi_stag = Grad_par_EM(phi_aligned, phi_stag, psi, CELL_YLOW);
  }
  ddt(chiU) = - bracket(phi_stag, U, bm, CELL_YLOW)
              - Vpar_Grad_par_EM(U_aligned, U_aligned, U, U, psi, CELL_YLOW)
              - Grad_par_phi_stag
              - (nu_parallel/mu)*UmV
              - (U*S_stag)/n_stag
              + 0.71*Grad_par_T_stag;

  if(diff_perp_U > 0.) {
    ddt(chiU) += diff_perp_U*Delp2(U);
  }

  set_lower_ddt_zero(chiU) ;

  // ***Electron parallel velocity Equation***
  if (hydrodynamic) {
    ddt(chiV) = ddt(chiU);
  }
  else {
    ddt(chiV) = - bracket(phi_stag, V, bm, CELL_YLOW)
                - Vpar_Grad_par_EM(V_aligned, V_aligned, V, V, psi, CELL_YLOW)
                + mu*Grad_par_phi_stag
                + nu_parallel*UmV
                - mu*T_stag*Grad_par_EM(logn_aligned, logn_stag, psi, CELL_YLOW)
                - 1.71*mu*Grad_par_T_stag
                - (V*S_stag)/n_stag ;

    if(diff_perp_V > 0.) {
      ddt(chiV) += diff_perp_V*Delp2(V);
    }

    set_lower_ddt_zero(chiV) ;
  }

  if (!isothermal){

    //Calculate parallel heat flux
    //Equivalent to qpar = -2.0*kappa0*Grad_par(pow(T,7.0/2.0),CELL_YLOW)/7.0 - 0.71*n_stag*T_stag*UmV;
    //This has been rewritten as follows to speed up the code.
    qpar_aligned = toFieldAligned(
        T_stag*((-kappa0)*powT_1_5_stag*Grad_par_T_stag - 0.71*n_stag*UmV), "RGN_NOBNDRY");

    //Set heat flux boundary condition
    mesh->communicate(qpar_aligned);
    qpar_aligned.applyBoundary(time);
    if (!symmetry_plane){
      qsheath_ydown_staggered(qpar_aligned, T_sheath_lower, n_sheath_lower, V_sheath_lower, mu);
    }
    qsheath_yup_staggered(qpar_aligned, T_sheath_upper, n_sheath_upper, V_sheath_upper, mu);
    if (mesh->yend - mesh->ystart + 1 < 3) {
      // Need to communicate mesh->ystart point that has just been set by boundary conditions
      mesh->communicate(qpar_aligned);
    }
    if (electromagnetic) {
      qpar_centre = stormInterp_to(qpar_aligned, CELL_CENTRE);
      mesh->communicate(qpar_centre);
      qpar_centre.applyBoundary(time);
    }

    //***Electron temperature Equation***
    Tcoef = (2.0/3.0)/p;
    ddt(logT) = - bracket(phi,logT,bm,CELL_CENTRE)
             - Vpar_Grad_par_EM(V_aligned, logT_aligned, V_centre, logT, psi_centre,
                   CELL_CENTRE)
             - Tcoef*Div_par_EM(qpar_aligned, qpar_centre, psi_centre,CELL_CENTRE)
             - 2./3.*0.71*UmV_centre*Grad_par_EM(logT_aligned, logT, psi_centre, CELL_CENTRE)
             - (2./3.)*div_par_V
             + 2.0/(3.0*mu)*nu_parallel0*(pow(T, -2.5, "RGN_NOBNDRY"))*n*SQ(UmV_centre)
             + Tcoef*S_E
             + S*SQ(V_centre)/(3.*mu*p)
             - S/n ;

    if(uniform_diss_paras){
      ddt(logT) += Tcoef*kappa0_perp*Delp2(T);
    }
    else{
      ddt(logT) += Tcoef*(kappa_perp*Delp2(T) + Grad_perp_dot_Grad_perp(kappa_perp, T));
    }

    //Curvature terms for electron temperature
    ddt(logT) += -2./3.*Curv_phi //Divergence of ExB
      + 2./3.*(Curv_p/n + (5./2.)*Curv_T) //Divergence of diamagnetic flow and diamagnetic heat flux
      - SQ(V_centre)*Curv_p/(3.*mu*p); //gyro-viscous energy transfer term

    if (run_1d) {
      // assume we are only setting up equilibrium:
      // slow down the T evolution so we can take longer timesteps
      ddt(logT) /= run_1d_T_slowdown;
    }
  }

  return 0;
}

