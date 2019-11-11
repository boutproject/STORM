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
#include <interpolation.hxx>
#include <field_factory.hxx>
#include <invert_parderiv.hxx>
#include <bout/constants.hxx>

BOUTMAIN(STORM);

int STORM::init(bool restarting) {
  
  globalOptions = Options::getRoot();
  options = Options::getRoot()->getSection("storm");
  
  OPTION(options, mu_n0,                  -1.0) ;
  OPTION(options, mu_vort0,               -1.0) ;
  OPTION(options, mu,                     -1.0) ;
  OPTION(options, nu_parallel0,           -1.0) ;
  OPTION(options, g0,                     -1.0) ;  
  OPTION(options, kappa0,                 -1.0) ;
  OPTION(options, kappa0_perp,            -1.0) ;
  OPTION(options, phi_wall,                  0) ;
  
  OPTION(options, isothermal,            false) ;
  OPTION(options, uniform_diss_paras,    false) ;
  OPTION(options, run_1d,                false) ;
  OPTION(options, symmetry_plane,         true) ;
  OPTION(options, add_blob,               true) ;
  OPTION(options, verbose,               false) ;
  
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

  Vsheath_BC_prefactor = sqrt((mu/TWOPI)/(1.0+(1.0/mu))) ;
  Usheath_BC_prefactor = sqrt(1.0/(1.0+(1.0/mu))) ;
  
  // Get from options or calculate Lx, Ly and Lz
  set_Lx_Ly_Lz();

  //********* Metric, Bracket scheme and staggering of fields ***********

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
  n_stag.setLocation(CELL_YLOW) ;
  S_stag.setLocation(CELL_YLOW) ;
  phi_stag.setLocation(CELL_YLOW) ;  
  nu_parallel.setLocation(CELL_YLOW) ;
  qpar.setLocation(CELL_YLOW);
  T_stag.setLocation(CELL_YLOW);
  
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

  //************* Definition of the fileds to solve and Laplacian solver *************
  
  if(isothermal){
    SOLVE_FOR4(n, U, V, vort) ;
    SAVE_ONCE(T) ;
  }
  else{ SOLVE_FOR5(n,U,V,vort,T) ;
        SAVE_REPEAT(qpar);
  }
	  
  SAVE_REPEAT(phi);
  SAVE_ONCE(S) ;
  if (!isothermal) {
    SAVE_ONCE(S_E);
  }
  	
  phiSolver = Laplacian::create();
  phi = phi_wall ;  // Sets phi field to phi_wall for inital 'guess' of Laplace solver
  phi.setBoundary("phi") ;
  phi_stag.setBoundary("phi_stag") ;
 
  comms.add(n) ;  
  comms.add(vort) ; 
  qpar = 0.0;
  qpar.setBoundary("qpar");
  comms.add(T) ;
  
  //**************** Initialisation of the fields *******************
  
  //Initialise backgrounds (if necessary)
  initialise_background(restarting) ;
  
  if (((mesh->firstX() && phi_array_inner == NULL)
        || (mesh->lastX() && phi_array_outer == NULL))) {
    throw BoutException("Error: phi_array_inner or phi_array_outer has not been set.");
  }
    
  //*************** Filament initialisation, if not in background ***************

  if (!restarting) {
    if (add_blob){
      initialise_blob("blob");
    }
  }
       
  if(!isothermal){
    p = n*T;
  }else{
    p = n;
  }
  UmV = U-V;

  if (run_1d) {
    OPTION(options, run_1d_T_slowdown, 20.) ;
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
  output.write("\nDimensionless Parameters:") ;  
  output.write("\n\tmu_n0        = %e ", mu_n0) ; 
  output.write("\n\tmu_vort0     = %e ", mu_vort0) ; 
  output.write("\n\tnu_parallel0 = %e ", nu_parallel0) ; 
  output.write("\n\tg0           = %e ", g0) ; 
  output.write("\n\tmu           = %e ", mu) ;
  output.write("\n\tkappa_par    = %e ", kappa0);
  output.write("\n\tkappa_perp   = %e ", kappa0_perp);
  output.write("\nDimensionless Lengths:") ;  
  output.write("\n\tdx           = %e ", mesh->getCoordinates()->dx(0,0));
  output.write("\n\tdy           = %e ", mesh->getCoordinates()->dy(0,0));
  output.write("\n\tdz           = %e ", mesh->getCoordinates()->dz);
  output.write("\n\tzlength      = %e ", mesh->getCoordinates()->zlength());
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

  return 0;
}

////////////////////////////////////////////////////////////////////////

int STORM::rhs(BoutReal time) { 

  if (isothermal){
    T = 1.0 ;
  }

  if(verbose){ output<<"\r\t\t\t\t\t\t\t t = "<<time<<std::flush;} //Stream simulation time to screen

  if (mesh->yend - mesh->ystart + 1 < 3) {
    // Need additional communication before applying 'free_o3' boundary
    // conditions because 3 grid points are needed to apply them
    mesh->communicate(comms);
    // Re-apply boundary conditions explicitly after communicating
    n.applyBoundary(time);
    vort.applyBoundary(time);
    if (!isothermal) {
      T.applyBoundary(time);
    }
  }
  mesh->communicate(comms) ;

  if(!isothermal){
    p = n*T;
  }else{
    p = n;
  }
  sqrtT = sqrt(T); // Need boundary guard cells, because we will take derivatives of sqrtT
  
  //********** Extrapolate BC for density and temperature to calculate BC on V, U and q_par **********

  n_sheath_upper = extrap_sheath_upper(n);
  T_sheath_upper = extrap_sheath_upper(T);
  if(!symmetry_plane){
    n_sheath_lower = extrap_sheath_lower(n);
    T_sheath_lower = extrap_sheath_lower(T);
  }

  //********* Laplace inversion: calculation of phi, including value at sheath *********

  if (run_1d) {
    phisolver_1d() ; 
  }
  else{
    set_xguards(phi, phi_array_inner, phi_array_outer);
    phi = phiSolver->solve(vort, phi) ;
  }
	
  if (mesh->yend - mesh->ystart + 1 < 3) {
    // need to communicate before applying 'free_o3' boundary conditions because 3 grid points are needed to apply them
    mesh->communicate(phi);
  }
  phi.applyBoundary(time);
  mesh->communicate(phi);  
  
  phi_stag = interp_to(phi, CELL_YLOW, "RGN_NOBNDRY");
  // Communicate before applying boundary condition in case there are not
  // enough grid points in the x-direction for free boundary conditions.
  // Don't need y-derivatives of phi_stag so don't need to communicate
  // afterward to set yup/ydown
  mesh->communicate(phi_stag) ;
  phi_stag.applyBoundary(time) ;

  // Extrapolate phi at boundary for V boundary condition
  phi_sheath_lower = extrap_sheath_lower(phi) ;
  phi_sheath_upper = extrap_sheath_upper(phi) ;
  
  //********** Sheath BC on the ion and electron velocity, U and V ********
  
  if (mesh->yend - mesh->ystart + 1 < 3) {
    // Need to communicate so lower boundary point gets communicated to next processor for interpolation in y-direction
    // Need to communicate before calling Vsheath_* Usheath_* functions in case there are fewer than 3 grid points.
    mesh->communicate(U, V);
  }

  if (!symmetry_plane) {
    Vsheath_ydown_staggered(V, phi_sheath_lower, phi_wall, T_sheath_lower, Vsheath_BC_prefactor) ;
    Usheath_ydown_staggered(U, sqrt(T_sheath_lower), Usheath_BC_prefactor);
  }
  Vsheath_yup_staggered(V, phi_sheath_upper, phi_wall, T_sheath_upper, Vsheath_BC_prefactor);
  Usheath_yup_staggered(U, sqrt(T_sheath_upper), Usheath_BC_prefactor);
  if (run_1d) {
    V = U; // Set them equal even though we already set ddt(V)=ddt(U) so that the boundary values are exactly equal
  }

  mesh->communicate(U, V);

  V_centre = interp_to(V, CELL_CENTRE, "RGN_NOBNDRY");
  UmV = U-V;
  UmV_centre = interp_to(UmV, CELL_CENTRE, "RGN_NOBNDRY");
  
  //*********** Calculate staggered interpolation for n and T, evaluate parameters depending on n and T ************
  n_stag = interp_to(n, CELL_YLOW, "RGN_NOBNDRY");
  T_stag = interp_to(T, CELL_YLOW, "RGN_NOBNDRY");

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
    Curv_p = T*Curv_n + n*Curv_T;
  }
  
  //****************** Equations of the model **************************
  
  // ***Vorticity Equation***
  ddt(vort) = - bracket(phi, vort, bm, CELL_CENTRE) 
            + Div_par(UmV, CELL_CENTRE) 
            + UmV_centre*Grad_par(n)/n
            - Vpar_Grad_par(U, vort) ;

  if(uniform_diss_paras){
    ddt(vort) += mu_vort0*Delp2(vort) ;
  }
  else{
   ddt(vort) +=  mu_vort*Delp2(vort) + Grad_perp(mu_vort)*Grad_perp(vort) ;
  }

  // Curvature terms for vorticity
  ddt(vort) += Curv_p/n;

  // ***Density Equation***
  ddt(n) = - bracket(phi, n, bm, CELL_CENTRE) 
  	   - n*Div_par(V, CELL_CENTRE) 
  	   - Vpar_Grad_par(V,n)
  	   + S ;
  
  if(uniform_diss_paras){
    ddt(n) += mu_n0*Delp2(n) ;
  }
  else{
    ddt(n) += mu_n*Delp2(n) + Grad_perp(mu_n)*Grad_perp(n) ;
  }

  // Curvature terms for density 
  ddt(n) += -n*Curv_phi + Curv_p;
  
  // ***Ion parallel velocity Equation*** 
  Grad_par_phi = Grad_par(phi, CELL_YLOW);
  Grad_par_T   = Grad_par(T,CELL_YLOW);
  ddt(U) = - bracket(phi_stag, U, bm, CELL_YLOW)
           - Vpar_Grad_par(U, U) 
           - Grad_par_phi
           - (nu_parallel/mu)*UmV
           - (U*S_stag)/n_stag
           + 0.71*Grad_par_T;

  // ***Electron parallel velocity Equation***
  if (run_1d) {
    ddt(V) = ddt(U);
  }
  else {
    ddt(V) = - bracket(phi_stag, V, bm, CELL_YLOW)
             - Vpar_Grad_par(V, V)
             + mu*Grad_par_phi
             + nu_parallel*UmV
             - (mu/n_stag)*Grad_par(p, CELL_YLOW)
             - 0.71*mu*Grad_par_T
             - (V*S_stag)/n_stag ;
  }

  if (!isothermal){

    //Calculate parallel heat flux
    //Equivalent to qpar = -2.0*kappa0*Grad_par(pow(T,7.0/2.0),CELL_YLOW)/7.0 - 0.71*n_stag*T_stag*UmV;
    //This has been rewritten as follows to speed up the code.
    qpar = T_stag*((-kappa0)*powT_1_5_stag*Grad_par_T - 0.71*n_stag*UmV);

    //Set heat flux boundary condition
    V_sheath_lower = sliceXZ(V, mesh->ystart) ;
    V_sheath_upper = sliceXZ(V, mesh->yend+1) ;
    U_sheath_lower = sliceXZ(U, mesh->ystart) ;
    U_sheath_upper = sliceXZ(U, mesh->yend+1) ;

    if (mesh->yend - mesh->ystart + 1 < 3) {
      // Need to communicate so lower boundary point gets communicated to next processor for interpolation in y-direction
      // Need to communicate before calling qsheath_* functions in case there are fewer than 3 grid points.
      mesh->communicate(qpar);
    }
    if (!symmetry_plane){
      qsheath_ydown_staggered(qpar, T_sheath_lower, n_sheath_lower, U_sheath_lower, V_sheath_lower, mu);
    }
    qsheath_yup_staggered(qpar, T_sheath_upper, n_sheath_upper, U_sheath_upper, V_sheath_upper, mu);
    qpar.applyBoundary(time);
    mesh->communicate(qpar);

    //***Electron temperature Equation***
    Tcoef = 2.0/(3.0*n);
    ddt(T) = - bracket(phi,T,bm,CELL_CENTRE)
             - Vpar_Grad_par(V,T)
             - Tcoef*Div_par(qpar,CELL_CENTRE)
             - 2./3.*0.71*UmV_centre*Grad_par(T)
             - (2./3.)*T*Div_par(V,CELL_CENTRE)
             + 2.0/(3.0*mu)*nu_parallel0*(pow(T, -1.5, "RGN_NOBNDRY"))*n*SQ(UmV_centre)
             + Tcoef*S_E
             + S*SQ(V_centre)/(3.*mu*n)
             - T*S/n ;

    if(uniform_diss_paras){
      ddt(T) += Tcoef*kappa0_perp*Delp2(T) ;
    }
    else{
      ddt(T) += Tcoef*(kappa_perp*Delp2(T) + Grad_perp(kappa_perp)*Grad_perp(T));
    }

    //Curvature terms for electron temperature
    ddt(T) += -2./3.*T*Curv_phi //Divergence of ExB
      + 2./3.*T*(Curv_p/n + (5./2.)*Curv_T) //Divergence of diamagnetic flow and diamagnetic heat flux
      + Tcoef*(0.5/mu)*SQ(V_centre)*Curv_p; //gyro-viscous energy transfer term

    if (run_1d) {
      // assume we are only setting up equilibrium:
      // slow down the T evolution so we can take longer timesteps
      ddt(T) /= run_1d_T_slowdown;
    }
  }

  return 0;
}


