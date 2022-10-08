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
#include "storm_version.hxx"

#include <vector>
#include <numeric>

class STORM;

BOUTMAIN(STORM);

int STORM::init(bool restarting) {

  // Add version information to output files
  dump.setAttribute("", "storm_git_hash", storm_git_hash);
  dump.setAttribute("", "storm_git_diff", storm_git_diff);
  dump.setAttribute("", "storm_configs_git_hash", storm_configs_git_hash);
  dump.setAttribute("", "storm_configs_git_diff", storm_configs_git_diff);
  dump.setAttribute("", "bout_configs_git_diff", bout_configs_git_diff);
  SAVE_ONCE(storm_cmake_cache);
  SAVE_ONCE(module_list);
  
  coordinates_centre = mesh->getCoordinates(CELL_CENTRE);
  
  OPTION(options, mu_n0_scalar,           -1.0) ;
  OPTION(options, mu_vort0_scalar,        -1.0) ;
  OPTION(options, mu,                     -1.0) ;
  OPTION(options, nu_parallel0_scalar,    -1.0) ;
  OPTION(options, g0,                     -1.0) ;  
  OPTION(options, kappa0,                 -1.0) ;
  OPTION(options, kappa0_perp_scalar,     -1.0) ;
  OPTION(options, diff_perp_U,             -1.) ;
  OPTION(options, diff_perp_V,             -1.) ;
  OPTION(options, beta0,                  -1.0) ;
  OPTION(options, phi_wall,                  0) ;
  OPTION(options, old_phi_wall_value,    false) ;
  
  OPTION(options, isothermal,            false) ;
  OPTION(options, boussinesq,                1) ; // 0 = no Boussinesq, 1 = standard STORM, 
                                                  // 2 = Hermes like, 3 = STORM + source term, 4 = STORM - vort*d(logn)/dt
  OPTION(options, split_n0,              false) ; // split n = 0 component and solve it with laplaceXY
  OPTION(options, electromagnetic,       false) ;
  OPTION(options, uniform_diss_paras,    false) ;
  OPTION(options, S_in_peq,               true) ;
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
  
  OPTION(options, equilibrium_source, "1d_profiles"); // possible values: 1d_profiles, input_file, profiles_file, single null, double null
  OPTION(options, equilibrium_file_path, "");
  if (equilibrium_file_path == "") {
    // The user has not explicitly set equilibrium_file_path, use path relative to simulation output directory instead
    std::string equilibrium_directory;
    OPTION(options, equilibrium_directory, "equilibrium");
    std::string data_dir = globalOptions["datadir"].withDefault("data");
    equilibrium_file_path = data_dir + "/" + equilibrium_directory;
  }
  OPTION(options, equilibrium_data_file, "");

  OPTION(options, average_radial_boundaries_core_SOL, false) ;
  OPTION(options, realistic_geometry,                "none") ; // options: none (slab), salpha (circular), singlenull, doublenull
  OPTION(options, normalise_all,                      false) ;
  OPTION(options, monitor_minmaxmean,                 false) ;
  OPTION(options, sources_realistic_geometry,         false) ;
  OPTION(options, sources_realisticgeometry_background, false) ;
  OPTION(options, increased_dissipation_xbndries,        -1) ;
  OPTION(options, xbndry_dissipation_factor,            10.) ;
  OPTION(options, increased_resistivity_xbndries,     false) ;
  OPTION(options, increased_resistivity_core,         false) ;
  OPTION(options, normalise_sources, realistic_geometry == "none") ;

  // Set default values for boundary conditions for all variables.
  // Can be overridden in input file, but should never need to be.
  // Needs to be called before adding anything to the time-solver or using setBoundary().
  setBoundaryConditionsOptions();

  // Check if we're using ShiftedMetric to run simulation in x-z orthogonal coordinates
  auto& mesh_options = globalOptions["mesh"];
  isshifted = false;
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
    isshifted = true;
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
  if (nu_parallel0_scalar < 0 ){
    nu_parallel0_scalar = 0.51*nu_ei/Omega_i;
  }
  if(mu_n0_scalar < 0){
    mu_n0_scalar = (1.0 + 1.3*pow(q, 2))*(1.0 + T_i0/T_e0)*pow(rho_e, 2)*nu_ei ; // Neo-Classical Particle diffusion, m^2s-1
    mu_n0_scalar = mu_n0_scalar/(SQ(rho_s)*Omega_i);
  }
  if(mu_vort0_scalar < 0){
    mu_vort0_scalar = (1.0 + 1.6*pow(q, 2))*(6.0/8.0)*pow(rho_i, 2)*nu_ii ;      // Neo-Classical Ion viscosity, m^2s-1
    mu_vort0_scalar = mu_vort0_scalar/(SQ(rho_s)*Omega_i);
  }
  if(g0 < 0){
    g0 = 2.0*rho_s/R_c ; 
  }
  
  if(kappa0 < 0 ){
    kappa0 = 3.16*(SQ(V_the)/nu_ei)/(SQ(rho_s)*Omega_i);
  }

  if(kappa0_perp_scalar < 0 ){
    kappa0_perp_scalar = (1.0 + 1.6*SQ(q))*(4.66)*SQ(rho_e)*nu_ei/(SQ(rho_s)*Omega_i);
  }

  if(beta0 < 0) {
    // ratio of nominal plasma pressure to nominal magnetic pressure
    beta0 = 2. * mu_0 * n_0 * T_e0 / SQ(B_0) ;
  }

  Vsheath_BC_prefactor = sqrt((mu/TWOPI)/(1.0+(1.0/mu))) ;
  Usheath_BC_prefactor = sqrt(1.0/(1.0+(1.0/mu))) ;
  
  // Get from options or calculate Lx, Ly and Lz
  if (realistic_geometry == "none") {
    set_Lx_Ly_Lz();
  } else {
    Lx = NAN;
    Ly = NAN;
    Lz = NAN;
  }

  //********* Metric, Bracket scheme and staggering of fields ***********
  if (isshifted && mesh->LocalNz % 2 == 0 && mesh->yend - mesh->ystart + 1 < 3) {
    throw BoutException("Shifted coordinates with an even Nz number of points and mesh->yend - mesh->ystart + 1 < 3 is not currently implemented.");
  }

  if (  realistic_geometry == "salpha"
     || realistic_geometry == "singlenull"
     || realistic_geometry == "doublenull") {
    std::string coordinates_type = "";
    if (!mesh->get(coordinates_type, "coordinates_type")) {
      if (coordinates_type != "orthogonal") {
        throw BoutException("Incorrect coordinate system type for realistic geometry, must be orthogonal.");
      }
    } 
    else {
      if (!mesh->get(coordinates_type, "parallel_transform")) {
        if (coordinates_type != "shiftedmetric") {
          throw BoutException("Incorrect parallel transform type for realistic geometry, must be shiftedmetric.");
        }
      }
      else {
        throw BoutException("Realistic geometry, but the grid doesn't have coordinates_type nor parallel_transform.");
      }
    }

    if (!isshifted) {
      throw BoutException("Realistic geometry requires shifted parallel transform!.");
    }

    bxcv.covariant = false; // Read contravariant components
    GRID_LOAD(bxcv);        // Specified components of b0 x kappaa
    Curlb_B = 2.*bxcv/coordinates_centre->Bxy;
    B2 = SQ(coordinates_centre->Bxy);
    G3 = coordinates_centre->G3;
  }
  else if (realistic_geometry != "none") {
    throw BoutException("realistic_geometry must be either none, salpha, singlenull, or doublenull!");
  }

  // Normalise grid mesh.   
  coordinates_stag = mesh->getCoordinates(CELL_YLOW);

  if (normalise_lengths){
    for (auto* coords : {coordinates_centre, coordinates_stag}) {
      coords->g_11 /= SQ(rho_s);
      coords->g_22 /= SQ(rho_s);
      coords->g_33 /= SQ(rho_s);
      coords->g_12 /= SQ(rho_s);
      coords->g_13 /= SQ(rho_s);
      coords->g_23 /= SQ(rho_s);
      coords->g11 *= SQ(rho_s);
      coords->g22 *= SQ(rho_s);
      coords->g33 *= SQ(rho_s);
      coords->g12 *= SQ(rho_s);
      coords->g13 *= SQ(rho_s);
      coords->g23 *= SQ(rho_s);

      // Pass 'false' argument so we use recalculate_staggered=false, since we
      // do not want to interpolate the staggered Coordinates from the
      // CELL_CENTRE ones, in case they were loaded from a grid file.
      coords->geometry(false);
    }
    if (realistic_geometry != "none") {
      G3 = coordinates_centre->G3;
    }
  }

  if (normalise_all) {
    if (normalise_lengths) {
      throw BoutException("Error: not possible to use normalise_lengths = true and normalise_all = true.");
    }
    for (auto* coords : {coordinates_centre, coordinates_stag}) {
      coords->Bxy /= B_0;
      coords->g_11 *= SQ(rho_s*B_0);
      coords->g_22 /= SQ(rho_s);
      coords->g_33 /= SQ(rho_s);
      coords->g_12 *= B_0;
      coords->g_13 *= B_0;
      coords->g_23 /= SQ(rho_s);
      coords->g11 /= SQ(rho_s*B_0);
      coords->g22 *= SQ(rho_s);
      coords->g33 *= SQ(rho_s);
      coords->g12 /= B_0;
      coords->g13 /= B_0;
      coords->g23 *= SQ(rho_s);
      coords->dx /= SQ(rho_s)*B_0;
      coords->J *= B_0/rho_s;

      // Pass 'false' argument so we use recalculate_staggered=false, since we
      // do not want to interpolate the staggered Coordinates from the
      // CELL_CENTRE ones, in case they were loaded from a grid file.
      coords->geometry(false);
    }
    B2 = SQ(coordinates_centre->Bxy);
    if (realistic_geometry != "none") {
      // Load curvature operator
      bxcv.x /= B_0;
      bxcv.y *= SQ(rho_s);
      bxcv.z *= SQ(rho_s);
      Curlb_B = 2.*bxcv/coordinates_centre->Bxy;
      G3 = coordinates_centre->G3;
    }
  }
  else {
    B2 = 1.;
  }
  coordinates_stag->outputVars(dump);

  // Save additional coordinate objects
  // Note: Rxy, Zxy and psixy are *not* normalised
  if (mesh->get(Rxy, "Rxy") == 0) {
    mesh->communicate(Rxy);
    SAVE_ONCE(Rxy);
  }
  if (mesh->get(Zxy, "Zxy") == 0) {
    mesh->communicate(Zxy);
    SAVE_ONCE(Zxy);
  }
  if (mesh->get(psixy, "psixy") == 0) {
    mesh->communicate(psixy);
    SAVE_ONCE(psixy);
  }
  if (realistic_geometry != "none") {
    dump.addOnce(bxcv.x, "bxcvx");
    dump.addOnce(bxcv.y, "bxcvy");
    dump.addOnce(bxcv.z, "bxcvz");
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
  if (!sources_realistic_geometry) {
    S = FieldFactory::get()->create3D("S:function", &globalOptions, mesh, CELL_CENTRE, 0);
    S_stag = FieldFactory::get()->create3D("S:function", &globalOptions, mesh, CELL_YLOW, 0);
    if (!isothermal) {
      S_E = FieldFactory::get()->create3D("S_E:function", &globalOptions, mesh, CELL_CENTRE, 0);
    }
    
    if (normalise_sources) {
      // Normalise S and S_stag by Ly, so that we do not need Ly in the input file
      S /= Ly;
      S_stag /= Ly;
      if (!isothermal) {
        S_E /= Ly;
      }
    }
  }
  else {
    set_sources_realistic_geometry();
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
  SAVE_ONCE(S_stag);
  if (isothermal) {
    T = 1.;
    logT_aligned = toFieldAligned(logT, "RGN_NOX");
    logT_aligned = 0.;
    T_aligned = exp(logT_aligned, "RGN_NOX");
    logT_stag = stormInterp_to(logT_aligned, CELL_YLOW);
    logT_stag = 0.;
    T_stag = exp(logT_stag, "RGN_NOY");
    sqrtT = 1.;
    if (mesh->hasBndryUpperY()) {
      T_sheath_upper = exp(extrap_sheath_upper(logT_aligned));
    }
    if (!symmetry_plane and mesh->hasBndryLowerY()) {
      T_sheath_lower = exp(extrap_sheath_lower(logT_aligned));
    }
    SAVE_ONCE(T);
  } else {
    SOLVE_FOR(logp);
    logT = logp - logn;
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

  SAVE_ONCE(mu_n0, mu_vort0, nu_parallel0, kappa0_perp, diff_perp_U, diff_perp_V);

  if (save_aligned_fields) {
    SAVE_REPEAT5(n_aligned, U_aligned, V_aligned, vort_aligned, phi_aligned);
    if (!isothermal) {
      SAVE_REPEAT2(T_aligned, qpar_aligned);
    }
    if (electromagnetic) {
      SAVE_REPEAT(psi_aligned);
    }
  }

  Options& phiOptions = globalOptions.isSection("phiSolver")
                        ? *globalOptions.getSection("phiSolver")
                        : *globalOptions.getSection("laplace");
  if (split_n0) {
    if (!evolving_bcs) {
      throw BoutException("split_n0 not implemented for evolving_bcs=false!");
    }
    if (boussinesq == 0) {
      throw BoutException("split_n0 not implemented for boussinesq == 0!");
    }

    // Set DC part to zero - this part will be solved by phiSolverxy
    phiOptions["global_flags"] = phiOptions["global_flags"].withDefault(1);
    // keep boundary flags set to 0 (zero-value Dirichlet), which is the BOUT++ default

    auto phixyOptions = globalOptions["laplacexy"];
    // hypre's algebraic multigrid method is usually the best preconditioner
    phixyOptions["pctype"] = phixyOptions["pctype"].withDefault("hypre");
    // Use finite-difference discretisation, instead of default finite-volume one
    phixyOptions["finite_volume"] = phixyOptions["finite_volume"].withDefault(false);

    phiSolverxy = new LaplaceXY(mesh, &phixyOptions);
    if (realistic_geometry != "none") {
      phiSolverxy->setCoefs(1./B2, 0.);
    }
    phi2D = 0.;
    restart.addOnce(phi2D,"phi2D");
    SAVE_REPEAT(phi2D);
  } else {
    // Set the radial boundary conditions to non-zero Dirichlet, with the value passed
    // in the boundary cells of the 'initial guess'
    phiOptions["inner_boundary_flags"]
      = phiOptions["inner_boundary_flags"].withDefault(16);
    phiOptions["outer_boundary_flags"]
      = phiOptions["outer_boundary_flags"].withDefault(16);
  }
  phiSolver = Laplacian::create(&phiOptions);
  phi = phi_wall ;  // Sets phi field to 0 for inital 'guess' of Laplace solver
  phi.setBoundary("phi") ;
  phi_aligned.setBoundary("phi_aligned") ;
  phi_stag.setBoundary("phi_stag") ;
  if (boussinesq == 0) {
    // Check that the Laplacian solver uses 3D coefficients
    ASSERT0(phiSolver->uses3DCoefs());
  }
  if (realistic_geometry != "none" && boussinesq > 0) {
    phiSolver->setCoefD(1./B2);
    phiSolver->setCoefC2(1./B2);
  }

  if (electromagnetic) {
    Options& psi_options = *globalOptions.getSection("psiSolver");
    psiSolver = Laplacian::create(&psi_options, CELL_YLOW);

    // Check that the Laplacian solver uses 3D coefficients
    ASSERT0(psiSolver->uses3DCoefs());

    psiSolver->setCoefD(-1.0);
    psi = 0.; // make sure all guard cells are set to something
    psi.setBoundary("psi");
    psi_aligned.setBoundary("psi_aligned");
    psi_centre.setBoundary("psi_centre");

    use_psi_boundary_solver = psi_options["use_psi_boundary_solver"].withDefault(true);
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

      std::string psi_solver_type = globalOptions["psiSolver"]["type"];
      if (psi_solver_type == "naulin") {
        // Actually only need the FFT sub-solver as at the boundary we solve Delp2(psi) =
        // n*(U-V) which has constant coefficients
        psiBoundarySolver = Laplacian::create(psi_options.getSection("delp2solver"), CELL_YLOW);

        // Set boundary flags from psiSolver options, like LaplaceNaulin does in its
        // constructor
        psiBoundarySolver->setInnerBoundaryFlags(
            psi_options["inner_boundary_flags"].as<int>());
        psiBoundarySolver->setOuterBoundaryFlags(
            psi_options["outer_boundary_flags"].as<int>());
      } else {
        psiBoundarySolver = Laplacian::create(&psi_options, CELL_YLOW);
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

  logn_aligned.setBoundary("logn_aligned");
  logn_stag.setBoundary("logn");
  vort_aligned.setBoundary("vort_aligned");
  if (!isothermal) {
    logT.setBoundary("logp");
    logT_aligned.setBoundary("logT_aligned");
    logT_stag.setBoundary("logp");
  }

  comms.add(logn);
  comms.add(vort) ;
  if (!isothermal) {
    qpar_aligned.setBoundary("qpar_aligned");
    if (electromagnetic) {
      qpar_centre.setBoundary("qpar_centre");
    }
    comms.add(logp) ;
  }
  if(boussinesq == 0){
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

  // set initial values of logn, logT, and logp
  logn = log(n);
  logT = log(T);
  logp = logn + logT;

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
  output.write("\nGit Version of STORM: %s", storm_git_hash) ;
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
  output.write("\n\tmu_n0        = %e ", mu_n0_scalar);
  output.write("\n\tmu_vort0     = %e ", mu_vort0_scalar);
  output.write("\n\tnu_parallel0 = %e ", nu_parallel0_scalar);
  output.write("\n\tg0           = %e ", g0) ; 
  output.write("\n\tmu           = %e ", mu) ;
  output.write("\n\tkappa_par    = %e ", kappa0);
  output.write("\n\tkappa_perp   = %e ", kappa0_perp_scalar);
  output.write("\n\tdiff_perp_U  = %e ", diff_perp_U);
  output.write("\n\tdiff_perp_V  = %e ", diff_perp_V);
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
    Options& fast_output_options = *globalOptions.getSection("fast_output");
    while (true) {
      // Add more points if explicitly set in input file
      xpos = fast_output_options["xpos"+std::to_string(i)].withDefault(-1.);
      ypos = fast_output_options["ypos"+std::to_string(i)].withDefault(-1.);
      zpos = fast_output_options["zpos"+std::to_string(i)].withDefault(-1.);
      if (xpos<0. || ypos<0. || zpos<0.) {
        output.write("\tAdded %i fast_output points\n", i);
        break;
      }
      ix = int(xpos*mesh->GlobalNx);
      iy = int(ypos*mesh->GlobalNy);
      iz = int(zpos*mesh->GlobalNz);
      fast_output.add("n"+std::to_string(i), n, ix, iy, iz);
      fast_output.add("T"+std::to_string(i), T, ix, iy, iz);
      fast_output.add("phi"+std::to_string(i), phi, ix, iy, iz);
      i++;
    }
  }

  // Increase dissipation near boundaries
  mu_n0 = mu_n0_scalar;
  mu_vort0 = mu_vort0_scalar;
  kappa0_perp = kappa0_perp_scalar;
  nu_parallel0 = nu_parallel0_scalar;
  if (increased_dissipation_xbndries > 0) {
    enhance_in_radial_buffers(mu_n0, increased_dissipation_xbndries,
                              xbndry_dissipation_factor);
    enhance_in_radial_buffers(mu_vort0, increased_dissipation_xbndries,
                              xbndry_dissipation_factor);
    enhance_in_radial_buffers(kappa0_perp, increased_dissipation_xbndries,
                              xbndry_dissipation_factor);
    if (increased_resistivity_xbndries) {
      enhance_in_radial_buffers(nu_parallel0, increased_dissipation_xbndries,
                                xbndry_dissipation_factor);
    }
  }

  if (mesh->getGlobalXIndex(mesh->xend) < 60 && mesh->periodicY(mesh->xend) && increased_resistivity_core) {
    nu_parallel0 = 0.1;
  }

  /////////////////////////////////////////////////////////
  // Neutral models
  
  TRACE("Initialising neutral models");
  neutrals = NeutralModel::create(solver, globalOptions["neutral"], dump);
  
  // Set normalisations
  if (neutrals != NULL) {
    neutrals->InitialiseNeutrals(T_e0/e, n_0, B_0, rho_s, Omega_i, mu);
    ionflux_lower = 0.;
    ionflux_upper = 0.;
  }

  // Preconditioner
  setPrecon((preconfunc)&STORM::precon);

  setup_history_tracking(restarting, storm_git_hash);

  return 0;
}

////////////////////////////////////////////////////////////////////////

int STORM::rhs(BoutReal time) { 

  if (first_rhs_evaluation) {
    first_rhs_evaluation = false;

    history_tracking_first_rhs();
  }

  if(verbose){ output<<"\r\t\t\t\t\t\t\t t = "<<time<<std::flush;} //Stream simulation time to screen

  rhs_counter++;

  mesh->communicate(comms) ;
  if (average_radial_boundaries_core_SOL) {
    average_Z_bndry(logn);
    average_YZ_bndry(logn, true, false, true);
    average_Z_bndry(vort);
    if (!isothermal) {
      average_Z_bndry(logp);
      average_YZ_bndry(logp, true, false, true);
    }
  }

  n = exp(logn, "RGN_NOY");
  if (!isothermal) {
    logT = logp - logn;
    T = exp(logT, "RGN_NOY");
  }

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
  if (!isothermal) {
    logT_stag = stormInterp_to(logT_aligned, CELL_YLOW);
  }
  if (electromagnetic) {
    // Derivatives of n_stag and T_stag are used only for electromagnetic runs
    if (isothermal) {
      mesh->communicate(logn_stag);
    }
    else {
      mesh->communicate(logn_stag, logT_stag);
      logT_stag.applyBoundary(time);
      if (average_radial_boundaries_core_SOL) {
        average_Z_bndry(logT_stag);
        average_YZ_bndry(logT_stag, true, false, true);
      }
      T_stag = exp(logT_stag, "RGN_NOY");
    }
    logn_stag.applyBoundary(time);
    if (average_radial_boundaries_core_SOL) {
      average_Z_bndry(logn_stag);
      average_YZ_bndry(logn_stag, true, false, true);
    }
    n_stag = exp(logn_stag, "RGN_NOY");
  }
  else {
    n_stag = exp(logn_stag, "RGN_NOBNDRY");
    if (!isothermal) {
      T_stag = exp(logT_stag, "RGN_NOBNDRY");
    }
  }

  if(!isothermal){
    p = n*T;
    sqrtT = sqrt(T, "RGN_NOY"); // Need x boundary guard cells, because we will take x derivatives of sqrtT
  }else{
    p = n;
  }
  
  //********** Extrapolate BC for density and temperature to calculate BC on V, U and q_par **********
  // Extrapolate in field-aligned coordinates since these are used for parallel boundary
  // conditions
  if (mesh->hasBndryUpperY()) {
    n_sheath_upper = exp(extrap_sheath_upper(logn_aligned));
  }
  if(!symmetry_plane and mesh->hasBndryLowerY()){
    n_sheath_lower = exp(extrap_sheath_lower(logn_aligned));
  }
  if (!isothermal) {
    if (mesh->hasBndryUpperY()) {
      T_sheath_upper = exp(extrap_sheath_upper(logT_aligned));
    }
    if(!symmetry_plane and mesh->hasBndryLowerY()){
      T_sheath_lower = exp(extrap_sheath_lower(logT_aligned));
    }
  }

  //********* Laplace inversion: calculation of phi, including value at sheath *********

  if (boussinesq == 0){
    if (realistic_geometry == "none") {
      phiSolver->setCoefD(n);
      phiSolver->setCoefC2(n);
    }
    else {
      phiSolver->setCoefD(n/B2);
      phiSolver->setCoefC2(n/B2);
    }
  } 

  if (hydrodynamic) {
    // in 1d or 2d this does nothing, but don't include ExB drifts if
    // hydrodynamic mode is used in 3d
    phi = 0.;
  }
  else{
    if (!split_n0) {
      if(!evolving_bcs){
        set_xguards(phi, phi_array_inner, phi_array_outer);
      }else{
        apply_bndry_phi();
      }
      if (realistic_geometry == "none") {
        phi = phiSolver->solve(vort, phi) ;
      }
      else {
        // Set G3 = 0. to avoid instabilities
        coordinates_centre->G3 = 0.;
        phi = phiSolver->solve(vort, phi) ;
        coordinates_centre->G3 = G3;
      }
    }
    else {
      // Solve for the n=0 component
      Field2D vort2D = DC(vort);
      if (mesh->firstX()) {
        int ixs  = mesh->xstart;
        int ixsg = ixs - 1;
        for (int ix=ixsg; ix>=0 ; --ix) {
          for (int iy=mesh->ystart; iy<=mesh->yend ; ++iy) {
            phi2D(ix,iy)  = phi_bc(ixs,iy);
          }
        }
      }
      if (mesh->lastX()) {
        int ixe = mesh->xend;
        int ixeg = ixe + 1;
        for (int ix=ixeg; ix<mesh->LocalNx ; ++ix) {
          for (int iy=mesh->ystart; iy <= mesh->yend ; ++iy) {
            phi2D(ix,iy)  = phi_bc(ixe,iy);
          }
        }
      }
      phi2D = phiSolverxy->solve(vort2D, phi2D);

      // Solve for the non-axisymmetric component
      if (realistic_geometry == "none") {
        phi = phiSolver->solve(vort-vort2D, phi) ;
      }
      else {
        // Set G3 = 0. to avoid instabilities
        coordinates_centre->G3 = 0.;
        phi = phiSolver->solve(vort-vort2D, phi) ;
        coordinates_centre->G3 = G3;
      }
      phi += phi2D;
    }
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
  if(boussinesq == 0){
    if (realistic_geometry != "none") {
      uE2 = Grad_perp_dot_Grad_perp(phi)/B2;
    } else {
      uE2 = Grad_perp_dot_Grad_perp(phi);
    }
    mesh->communicate(uE2);
    uE2.applyBoundary(time);
  }  

  // Extrapolate phi at boundary for V boundary condition
  if (not symmetry_plane and mesh->hasBndryLowerY()) {
    phi_sheath_lower = extrap_sheath_lower(phi_aligned) ;
  }
  if (mesh->hasBndryUpperY()) {
    phi_sheath_upper = extrap_sheath_upper(phi_aligned) ;
  }

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
  if(realistic_geometry == "salpha" && mesh->firstY() && mesh->getGlobalXIndex(mesh->xend+1) == ixseps1){
    // fix x-guard cells near the limiter for s-alpha geometry
    int ix = mesh->xend, iy = mesh->ystart;
    int ixp = ix+1;
    for(int iz = 0; iz<mesh->LocalNz; ++iz){
      U(ixp,iy,iz) = U(ix,iy,iz);
      V(ixp,iy,iz) = V(ix,iy,iz);
    }
  }
  U_aligned.applyBoundary(time);
  V_aligned.applyBoundary(time);
  
  //********** Sheath BC on the ion and electron velocity, U and V ********
  if (!symmetry_plane and mesh->hasBndryLowerY()) {
    Vsheath_ydown_staggered(V_aligned, phi_sheath_lower, phi_wall, T_sheath_lower, Vsheath_BC_prefactor) ;
    Usheath_ydown_staggered(U_aligned, sqrt(T_sheath_lower), Usheath_BC_prefactor);
  }
  if (mesh->hasBndryUpperY()) {
    Vsheath_yup_staggered(V_aligned, phi_sheath_upper, phi_wall, T_sheath_upper, Vsheath_BC_prefactor);
    Usheath_yup_staggered(U_aligned, sqrt(T_sheath_upper), Usheath_BC_prefactor);
  }
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
    if (realistic_geometry != "none") {
      mu_n /= B2;
      mu_vort /= B2;
      kappa_perp /= B2;
    }
  }

  // Precompute the curvature terms
  Curv_n = Curv(n, n_aligned);
  Curv_phi = Curv(phi, phi_aligned);
  if(isothermal){
    Curv_p = Curv_n;
  }else{
    Curv_T = Curv(T, T_aligned);
    Curv_p = T*Curv_n + n*Curv_T; // According to Ben has a better stability than Curv(p)
  }
  
  //****************** Equations of the model **************************
  
  // ***Vorticity Equation***
  if (hydrodynamic) {
    vorticity_equation["all"] = 0.;
  }
  else {
    vorticity_equation["ExB_advection"] = - bracket(phi, vort, bm, CELL_CENTRE);
    vorticity_equation["parallel_advection"] =
      - Vpar_Grad_par_EM(U_aligned, vort_aligned, U_centre, vort, psi_centre,
                         CELL_CENTRE);
    if (boussinesq == 1 || boussinesq == 3 || boussinesq == 4) {
      // Grad(n*Grad_perp(phi)) = n*Delp2(phi)
      vorticity_equation["parallel_current"] =
        Div_par_EM(UmV_aligned, UmV_centre, psi_centre, CELL_CENTRE)
        + UmV_centre*Grad_par_EM(logn_aligned, logn, psi_centre, CELL_CENTRE);
    }
    else if (boussinesq == 0 || boussinesq == 2) {
      // Grad(n*Grad_perp(phi)) = Delp2(phi) or non-boussinesq
       vorticity_equation["parallel_current"] =
         n*( Div_par_EM(UmV_aligned, UmV_centre, psi_centre, CELL_CENTRE)
             + UmV_centre*Grad_par_EM(logn_aligned, logn, psi_centre, CELL_CENTRE) );
    }
    else {
      throw BoutException("Boussinesq option unrecognized.");
    }

    if (boussinesq == 0) {
      vorticity_equation["bracket_uE2"] = -bracket(uE2/2.,n,bm); //check normalisation of last term
    }

    if (boussinesq == 3) {
      // assumes S axisymmetric
      vorticity_equation["density_source"] = - S*vort
                                             - coordinates_centre->g11*DDX(S)*DDX(phi)/B2;
    }

    if(uniform_diss_paras){
      vorticity_equation["perpendicular_diffusion"] = mu_vort0*Delp2(vort);
    }
    else{
     vorticity_equation["perpendicular_diffusion"] =
       mu_vort*Delp2(vort) + Grad_perp_dot_Grad_perp(mu_vort, vort);
    }

    // Curvature terms for vorticity
    if (boussinesq == 1 || boussinesq == 3 || boussinesq == 4) {
      vorticity_equation["curvature"] = Curv_p/n;
    }
    else {
      vorticity_equation["curvature"] = Curv_p;
    }
  }

  Field3D div_par_V = Div_par_EM(V_aligned, V_centre, psi_centre, CELL_CENTRE);
  // ***Density Equation***
  density_equation["ExB_advection"] = - bracket(phi, logn, bm, CELL_CENTRE);
  density_equation["parallel_advection"] =
    - div_par_V
    - Vpar_Grad_par_EM(V_aligned, logn_aligned, V_centre, logn, psi_centre, CELL_CENTRE);
  density_equation["density_source"] = S/n;
  
  if (uniform_diss_paras) {
    density_equation["perpendicular_diffusion"] = mu_n0*Delp2(n)/n;
  }
  else {
    // Include y derivatives??
    density_equation["perpendicular_diffusion"] = mu_n*Delp2(n)/n
                                                    + Grad_perp_dot_Grad_perp(mu_n, logn);
  }

  // Curvature terms for density 
  density_equation["ExB_curvature"] = -Curv_phi;
  density_equation["diamagnetic_curvature"] = Curv_p/n;

  if (boussinesq == 4) {
    vorticity_equation["vort*dndt"] = -vort*ddt(logn);
  }

  // ***Ion parallel velocity Equation*** 
  Grad_par_T_stag = Grad_par_EM(T_aligned, T_stag, psi, CELL_YLOW);
  if (hydrodynamic) {
    Grad_par_phi_stag = mu/(mu+1.)*(Grad_par_T_stag
                                    + T_stag*fromFieldAligned(
                                               Grad_par(logn_aligned, CELL_YLOW),
                                               "RGN_NOBNDRY")
                                   )
                        + 0.71*Grad_par_T_stag;
  } else {
    Grad_par_phi_stag = Grad_par_EM(phi_aligned, phi_stag, psi, CELL_YLOW);
  }
  ion_momentum_equation["ExB_advection"] = - bracket(phi_stag, U, bm, CELL_YLOW);
  ion_momentum_equation["parallel_advection"] =
    - Vpar_Grad_par_EM(U_aligned, U_aligned, U, U, psi, CELL_YLOW);
  ion_momentum_equation["grad_phi"] = - Grad_par_phi_stag;
  ion_momentum_equation["resistivity"] = - (nu_parallel/mu)*UmV;
  ion_momentum_equation["thermal_force"] = 0.71*Grad_par_T_stag;
  ion_momentum_equation["density_source"] = - (U*S_stag)/n_stag;

  if (diff_perp_U > 0.) {
    ion_momentum_equation["perpendicular_diffusion"] = diff_perp_U*Delp2(U);
  }

  // ***Electron parallel velocity Equation***
  if (hydrodynamic) {
    electron_momentum_equation["all"] = ddt(chiU);
  }
  else {
    electron_momentum_equation["ExB_advection"] = - bracket(phi_stag, V, bm, CELL_YLOW);
    electron_momentum_equation["parallel_advection"] =
      - Vpar_Grad_par_EM(V_aligned, V_aligned, V, V, psi, CELL_YLOW);
    electron_momentum_equation["grad_phi"] = mu*Grad_par_phi_stag;
    electron_momentum_equation["resistivity"] = nu_parallel*UmV;
    electron_momentum_equation["density_gradient"] =
      - mu*T_stag*Grad_par_EM(logn_aligned, logn_stag, psi, CELL_YLOW);
    electron_momentum_equation["temperature_gradient"] = - (1.71*mu)*Grad_par_T_stag;
    electron_momentum_equation["density_source"] = - (V*S_stag)/n_stag;

    if (diff_perp_V > 0.) {
      electron_momentum_equation["perpendicular_diffusion"] = diff_perp_V*Delp2(V);
    }

    set_lower_ddt_zero(chiV) ;
  }

  if (!isothermal){

    //Calculate parallel heat flux
    //Equivalent to qpar = -2.0*kappa0*Grad_par(pow(T,7.0/2.0),CELL_YLOW)/7.0 - 0.71*n_stag*T_stag*UmV;
    //This has been rewritten as follows to speed up the code.
    qpar_aligned = toFieldAligned(
        T_stag*((-kappa0)*powT_1_5_stag*Grad_par_T_stag - 0.71*n_stag*UmV), "RGN_NOBNDRY");
    mesh->communicate(qpar_aligned);

    //Set heat flux boundary condition
    qpar_aligned.applyBoundary(time);
    if (!symmetry_plane and mesh->hasBndryLowerY()){
      qsheath_ydown_staggered(qpar_aligned, T_sheath_lower, n_sheath_lower, V_sheath_lower, mu);
    }
    if (mesh->hasBndryUpperY()) {
      qsheath_yup_staggered(qpar_aligned, T_sheath_upper, n_sheath_upper, V_sheath_upper, mu);
    }
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
    electron_pressure_equation["ExB_advection"] = - bracket(phi,logp,bm,CELL_CENTRE);
    electron_pressure_equation["parallel_advection"] =
      - Vpar_Grad_par_EM(V_aligned, logn_aligned + logT_aligned, V_centre, logp,
                         psi_centre, CELL_CENTRE);
    electron_pressure_equation["parallel_conduction"] =
      - Tcoef*Div_par_EM(qpar_aligned, qpar_centre, psi_centre,CELL_CENTRE);
    electron_pressure_equation["thermal_force"] =
      - 2./3.*0.71*UmV_centre*Grad_par_EM(logT_aligned, logT, psi_centre, CELL_CENTRE);
    electron_pressure_equation["compression"] = - (5./3.)*div_par_V;
    // Resistive heating uses nu_parallel0_scalar so that we do not have extra resistive
    // heating in any radial buffer regions, which is added to nu_parallel0.
    electron_pressure_equation["resistive_heating"] =
      2.0/(3.0*mu)*nu_parallel0_scalar*(pow(T, -2.5, "RGN_NOBNDRY"))*n*SQ(UmV_centre);
    electron_pressure_equation["energy_source"] = Tcoef*S_E;

    if(uniform_diss_paras){
      electron_pressure_equation["perpendicular_diffusion"] = Tcoef*kappa0_perp*Delp2(T);
    }
    else{
      electron_pressure_equation["perpendicular_diffusion"] =
        Tcoef*(kappa_perp*Delp2(T) + Grad_perp_dot_Grad_perp(kappa_perp, T));
    }

    if (S_in_peq) {
      electron_pressure_equation["S*V^2/3/mu/p"] = S*SQ(V_centre)/(3.*mu*p);
    }

    //Curvature terms for electron temperature
    // Divergence of ExB
    electron_pressure_equation["ExB_curvature"] = -5./3.*Curv_phi;
    // Divergence of diamagnetic flow and diamagnetic heat flux
    electron_pressure_equation["diamagnetic_curvature"] = 5./3.*(Curv_p/n + Curv_T);
    // gyro-viscous energy transfer term
    electron_pressure_equation["gyroviscous_curvature"] = - SQ(V_centre)*Curv_p/(3.*mu*p);

    if (run_1d) {
      // assume we are only setting up equilibrium:
      // slow down the T evolution so we can take longer timesteps
      // assume that we solve for ddt(logp) = 0
      ddt(logp) = (ddt(logp) - ddt(logn))/run_1d_T_slowdown + ddt(logn);
    }
  }

  // Neutral gas
  if (neutrals) {
    TRACE("Neutral gas model");
   
    // Fluxes at targets for boundary conditions
    if (mesh->hasBndryLowerY()) {
      ionflux_lower = abs(n_sheath_lower*U_sheath_lower);
    }
    if (mesh->hasBndryUpperY()) {
      ionflux_upper = abs(n_sheath_upper*U_sheath_upper);
    }
    neutrals->setRecycledFlux(ionflux_lower, ionflux_upper);
 
    // Update neutral gas model
    neutrals->update(n, n_stag, U, V, T, T_stag, time);

    // Add sources/sinks to plasma equations
    density_equation["ionisation_recombination"] = -neutrals->S/n;

    Field3D Son_stag = neutrals->S_stag/n_stag;
    ion_momentum_equation["neutral_friction"] = -neutrals->Fi/n_stag;
    ion_momentum_equation["neutral_source"] = Son_stag*U;

    electron_momentum_equation["neutral_friction"] = -neutrals->Fe/n_stag;
    electron_momentum_equation["neutral_source"] = Son_stag*V;

    if (boussinesq == 2) {
      vorticity_equation["neutral_friction"] =
        - neutrals->Fperp*vort - Grad_perp_dot_Grad_perp(neutrals->Fperp, phi)/B2;
    }
    else if (boussinesq == 1 || boussinesq == 3 || boussinesq == 4) {
      vorticity_equation["neutral_friction"] =
        - (neutrals->Fperp*vort + Grad_perp_dot_Grad_perp(neutrals->Fperp, phi)/B2)/n;
    }
    else {
      vorticity_equation["neutral_friction"] =
        - neutrals->Fperp*vort/n - Grad_perp_dot_Grad_perp(neutrals->Fperp/n, phi)*n/B2;
    }

    if (!isothermal) {
      electron_pressure_equation["neutral_radiation"] = - Tcoef*neutrals->Rp;
    }

  }

  set_lower_ddt_zero(chiU);

  set_lower_ddt_zero(chiV);

  return 0;
}

int STORM::precon(BoutReal t, BoutReal gamma, BoutReal delta) {
  // Neutral gas preconditioning
  if (neutrals != NULL) neutrals->precon(t, gamma, delta);
  return 0;
}
