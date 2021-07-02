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

#include <bout/physicsmodel.hxx>
#include <smoothing.hxx>
#include <invert_laplace.hxx>
#include <initialprofiles.hxx>
#include <derivs.hxx>
#include <bout/constants.hxx>

class STORM2D : public PhysicsModel {
public:
  virtual ~STORM2D() {}
protected:
  int init(bool restarting);
  int rhs(BoutReal t);
private:
  Field3D n, T, vort ;                  // Evolving density, temp and vorticity
  Field3D phi ;
  Field3D S, S_E;                       // Density and energy source terms
  Field3D n_bg_source, T_bg_source ;    // Density and temperature sources to sustain a background of (n_bg, T_bg)
  Field2D sigma_n, sigma_T, sigma_vort; // ESEL model dissipation terms

  FieldGroup comms ;                    // Communicated variables

  BRACKET_METHOD bm;                    // Bracket method for advection terms

  BoutReal L ;                          // parallel connection length from midplane to target
  BoutReal lambda ;                     // sheath dissipation parameter = 1/L

  BoutReal mu_n0, mu_vort0, g0, kappa0_perp, kappa0 ; // mu
  BoutReal B_0, T_e0, T_i0, m_i, q, R_c, n_0, Z, loglambda ;
  BoutReal u, e, m_e, epsilon_0, mu_0;
  BoutReal c_s, rho_s, Omega_i, V_the, V_thi, Omega_e, rho_i, rho_e, nu_ei, nu_ii, nu_ee ;

  // Variables and functions used for parallel loss terms
  BoutReal mu, V_float ;                // floating potential for the nonisothermal sheath loss terms
  BoutReal V_sheath_prefactor, U_sheath_prefactor ;
  std::string SOL_closure ;                  // switch for 2d closure to use, either sheath_diss, vort_adv, or ESEL_like.
  bool sheath_linear ;                  // switch for linearised sheath conditions when using sheath dissipation closure
  bool initial_noise ;                  // switch to add random noise to density and vorticity to trigger instabilities quicker
  BoutReal n_bg, T_bg, phi_bg ;                 // background value to which fluctuations decay, as a fraction of the normalisation values
  const Field3D V_sheath(const Field3D &phi, const Field3D &T) ;
  const Field3D loss_n(const BoutReal &lambda, const Field3D &n, const Field3D &phi, const Field3D &T);
  const Field3D loss_vort(const BoutReal &lambda, const Field3D &vort, const Field3D &n, const Field3D &phi, const Field3D &T);
  const Field3D loss_T(const BoutReal &lambda, const Field3D &n, const Field3D &phi, const Field3D &T);

  const Field3D curv_op(const Field3D &f);

  bool isothermal ;                     // switch for isothermal simulations

  class Laplacian* phiSolver;           // Laplacian solver for vort -> phi
};

int STORM2D::init(bool UNUSED(restarting)) {
  
  Options *options = Options::getRoot()->getSection("storm");
  
  OPTION(options, mu_n0,                  -1.0) ;
  OPTION(options, mu_vort0,               -1.0) ;
  OPTION(options, g0,                     -1.0) ;  
  OPTION(options, kappa0_perp,            -1.0) ;
  OPTION(options, kappa0,                 -1.0) ;
  
  OPTION(options, B_0,                     0.5) ;   // Tesla 
  OPTION(options, T_e0,                     40) ;   // eV
  OPTION(options, T_i0,                     40) ;   // eV
  OPTION(options, m_i,                       2) ;   // Atomic Units
  OPTION(options, q,                         7) ;   // Dimensionless
  OPTION(options, R_c,                     1.5) ;   // m
  OPTION(options, n_0,                  0.8e19) ;   // m^-3
  OPTION(options, Z,                         1) ;   // Dimensionless
  OPTION(options, loglambda,                -1) ;   // Dimensionless
  
  Options::getRoot()->getSection("mesh")->get("Ly", L, 5500.);

  int bracket; 
  OPTION(options, bracket,                   2) ;
  OPTION(options, isothermal,            false) ;
  OPTION(options, SOL_closure,   "sheath_diss") ;
  OPTION(options, sheath_linear,          true) ;
  OPTION(options, initial_noise,         false) ;
  
  OPTION(options, n_bg,                    1.0) ;
  OPTION(options, T_bg,                    1.0) ;   

  // Specify Constants
  u = 1.66053892e-27 ;              // kg
  e = 1.602176565e-19 ;             // C
  m_e = 9.10938291e-31;             // kg
  epsilon_0 = 8.854187817e-12 ; 
  mu_0 = 4.*PI*1.e-7 ;
  
  // Convert Parameters to SI units 
  T_e0 = e*T_e0 ;                   // Joules
  T_i0 = e*T_i0 ;                   // Joules
  
  m_i = m_i*u ;                     // kg 
  
  // Calculate secondary parameters
  c_s = sqrt(T_e0/m_i) ;
  Omega_i = Z*e*B_0/m_i ; 
  rho_s = c_s/Omega_i ;  
  
  V_the = sqrt(T_e0/m_e) ; 
  V_thi = sqrt(T_i0/m_i) ; 
  
  Omega_e = e*B_0/m_e ; 
  rho_e = V_the/Omega_e ; 
  rho_i = V_thi/Omega_i ; 
  
  if (loglambda < 0){
    loglambda = 18.0 - log(sqrt(n_0/1.0e19)*pow(T_e0/(1000.0*e), -1.5)) ; 
  } 
  nu_ee = n_0*pow(e,4)*loglambda/(sqrt(m_e)*pow(epsilon_0, 2)*3.0*pow(TWOPI*T_e0, 1.5)) ;          // electron - electron collision frequency
  nu_ei = n_0*pow(Z,2)*pow(e,4)*loglambda/(sqrt(m_e)*pow(epsilon_0, 2)*3.0*pow(TWOPI*T_e0, 1.5)) ; // electron - ion collision frequency
  nu_ii = n_0*pow(Z*e,       4)*loglambda/(sqrt(m_i)*pow(epsilon_0, 2)*3.0*pow(TWOPI*T_i0, 1.5)*sqrt(2)) ; // ion - ion collision frequency 

  if(mu_n0 < 0){
    mu_n0 = (1.0 + 1.3*pow(q, 2))*(1.0 + T_i0/T_e0)*pow(rho_e, 2)*nu_ei ; // Neo-Classical Particle diffusion, m^2s-1
    mu_n0 = mu_n0/(rho_s*rho_s*Omega_i) ; 
  }
  if(mu_vort0 < 0){
    mu_vort0 = (1.0 + 1.6*pow(q, 2))*(6.0/8.0)*pow(rho_i, 2)*nu_ii ;      // Neo-Classical Ion viscosity, m^2s-1
    mu_vort0 = mu_vort0/(rho_s*rho_s*Omega_i) ; 
  }
  if(g0 < 0){
    g0 = 2.0*rho_s/R_c ; 
  }
  if(kappa0_perp < 0 ){
    kappa0_perp = (1.0 + 1.6*pow(q,2))*(4.66)*pow(rho_e, 2)*nu_ee/(rho_s*rho_s*Omega_i);
  }
  if(kappa0 < 0 ){
    kappa0 = 3.16*(pow(V_the, 2)/nu_ei)/(rho_s*rho_s*Omega_i);
  }
  
  // Floating potential (only used for nonisothermal cases)
  mu = m_i/m_e ;
  V_float = 0.5*log(TWOPI/mu) ;
  V_sheath_prefactor = sqrt((mu/TWOPI) / (1.0+(1.0/mu))) ;
  U_sheath_prefactor = sqrt(1.0 / (1.0+(1.0/mu))) ; 

  if(SOL_closure == "sheath_diss") {
    output << "\tClosure: sheath dissipation\n";
  }
  else if(SOL_closure == "vort_adv") { 
    output << "\tClosure: vorticity advection\n";
  }
  else if(SOL_closure == "ESEL_like") {
    output << "\tClosure: ESEL-like loss terms\n";  
  }
  else {
    output << "ERROR: Invalid choice of 2D closure method. Must be either sheath_diss, vort_adv, or ESEL_like.\n";
    return 1;
  }

  if (sheath_linear && !isothermal) {
    throw BoutException("Error: sheath_linear is not supported when isothermal=false");
  }

  if (SOL_closure == "ESEL_like") {
     // Set sinks from input file
    output << "\tSetting initial source and sink profiles\n";

    initial_profile("sigma_n", sigma_n);
    initial_profile("sigma_T", sigma_T);
    initial_profile("sigma_vort", sigma_vort);
    lambda = 0.0 ;
  }
  else{
    lambda = 1.0/L ;

    sigma_n = 0.0 ;
    sigma_vort = 0.0 ;
    sigma_T = 0.0 ;
  }
  
  output << "\tSetting initial source profiles\n";
  initial_profile("S", S);
  initial_profile("S_E", S_E);

  // Poisson brackets: b_hat x Grad(f) dot Grad(g) / B = [f, g]
  // Method to use: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE
  // Choose method to use for Poisson bracket advection terms
  
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
  
  if (isothermal){
    SOLVE_FOR2(n,vort) ;
  }
  else{
    SOLVE_FOR3(n, T, vort) ;
  }
  SAVE_REPEAT(phi);
  SAVE_ONCE(sigma_n) ;
  SAVE_ONCE(sigma_vort) ;
  SAVE_ONCE(sigma_T) ;
  
  phiSolver = Laplacian::create();
  // Starting phi
  if (!isothermal){
    // Laplace inversion flag must be = 16 for this to be used instead of 0
    phi = -V_float ;
    phi_bg = -V_float * T_bg ;
  }
  else{
    phi = 0.0 ;
    phi_bg = 0.0 ;
    T_bg = 1.0 ; // (to ensure background source term for density is correct even if user specifies T_bg != 1.0 for an isothermal simulation)
  }


  // Compute the background source terms
  if (n_bg > 0.0){
    n_bg_source = loss_n(lambda, n_bg, phi_bg, T_bg) ;
  }
  if (T_bg > 0.0){
    T_bg_source = loss_T(lambda, n_bg, phi_bg, T_bg) ;
  }

  // Initialise the fields
  initial_profile("n", n);
  if (isothermal) {
    T = 1.0;
  } else {
    initial_profile("T", T);
  }

  comms.add(n) ;
  if (!isothermal) comms.add(T) ;
  comms.add(vort) ; 

  // Seed turbulence with random noise
  if (initial_noise){
    output << "\tSeeding random noise for triggering turbulent instabilities\n";
    int MYPE = BoutComm::rank();
    srand (MYPE);
    for(int i=0; i < mesh->LocalNx ; i++){
      for(int k=0; k < mesh->LocalNz; k++){
         n(i,0,k)    += 2.*(((double) rand()/(RAND_MAX)) - 0.5)*0.0001;
         vort(i,0,k) += 2.*(((double) rand()/(RAND_MAX)) - 0.5)*0.00001;
         if (!isothermal){
           T(i,0,k)  += 2.*(((double) rand()/(RAND_MAX)) - 0.5)*0.0001;
         }
      }
    }
  }

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
  output.write("\n\tg0           = %e ", g0) ; 
  output.write("\n\tkappa_perp   = %e ", kappa0_perp) ;
  output.write("\n\tkappa        = %e ", kappa0) ;
  output.write("\n\tlambda       = %e ", lambda) ;
  output.write("\n\tV_float      = %e ", V_float) ;

  return 0;
}

int STORM2D::rhs(BoutReal UNUSED(time)) {

  // Communicate variables
  mesh->communicate(comms);

  // Solve for potential
  if (!isothermal){
    // Set the boundary conditions to the floating potential at each timestep using flag =16
    phi = phiSolver->solve(vort, -V_float);  
  }
  else{
    phi = phiSolver->solve(vort);
  }

  // Communicate phi
  mesh->communicate(phi);

  // Continuity equation:  
  ddt(n) = - bracket(phi,n,bm) - n*curv_op(phi) + curv_op(n*T) + mu_n0*Delp2(n) + S;
  
  // Choice of parallel loss terms for density
  ddt(n) -= loss_n(lambda, n, phi, T) ;

  // Source term to sustain background
  if (n_bg > 0.0){
    ddt(n) += n_bg_source ;
  }


  // Energy equation:
  if (!isothermal){
    ddt(T) = - bracket(phi,T,bm) - (2.0/3.0)*T*curv_op(phi) + (7.0/3.0)*T*curv_op(T) + (2.0/3.0)*(SQ(T)/n)*curv_op(n) + kappa0_perp*Delp2(T) + (2.0/3.0)*S_E/n - T*S/n;

    // Choice of parallel loss terms for temperature
    ddt(T) -= loss_T(lambda, n, phi, T) ;

    // Source term to sustain background
    // Background terms are divided by (n/n_bg) because they represent a source added to the energy equation, not the temperature equation
    if (T_bg > 0.0){
      ddt(T) += T_bg_source*(n_bg/n) ;
    }
  }
  

  // Vorticity equation:
  ddt(vort) = - bracket(phi, vort, bm) + curv_op(n*T)/n + mu_vort0*Delp2(vort);
  
  // Choice of parallel loss terms for vorticity
  ddt(vort) -= loss_vort(lambda, vort, n, phi, T) ;

  return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Loss operators and sheath velocity for different parallel closures.
// Fields (n, phi, vort, T) are supplied as arguments explicitly, while run parameters
// (SOL_closure, sheath_linear, unit_bg, sigma_*, isothermal) are supplied as global constants.
//
// "Sheath_dissipation" closure assumes no parallel gradients in (n,T,phi,vort), and that the electron and ion velocities integrated from mid-plane to target equal the Debye
// sheath expressions, with stagnation at the mid-plane. Parallel currents are allowed, and are set by the electron and ion velocities. 2D equations are same as 3D but averaged 
// in the parallel direction.
//
// "Vorticity_advection" closure assumes that there are no parallel currents, and that (n,vort,T) parallel advection terms are approximated by c_s/L_parallel * field. 2D equations
// are similar to 3D equations evaluated at the mid-plane.
//
// "ESEL_like" closure assumes (n,T,vort) decay away at a rate directly proportional to their value, and is included for comparison against earlier studies using ESEL and STORM2D.
//
// It will cause decay back to zero if phi = 0.0 (isothermal runs) or phi = -V_float (thermal runs).
// The extra terms imcluded with this option equal -loss_field(n=1,T=1,vort=0,phi=(0,-V_float)), so that the background term exactly cancels the loss. 

const Field3D STORM2D::V_sheath(const Field3D &phi, const Field3D &T){
  // Calculates the electron velocity at the sheath, which is needed in most of the closures.
  if (!isothermal){
    return V_sheath_prefactor*sqrt(T)*exp(-phi/T) ;
  }
  else{
    if (sheath_linear){
      return 1.0 - phi ;
    }
    else{
      return exp(-phi) ;
    }
  }
}


const Field3D STORM2D::loss_n(const BoutReal &lambda, const Field3D &n, const Field3D &phi, const Field3D &T){
// Parallel loss term for density equation

  Field3D n_loss ;
  n_loss.allocate() ;

  if(SOL_closure == "sheath_diss") {
    n_loss = lambda*n*V_sheath(phi, T) ;
  }
  else if(SOL_closure == "vort_adv") { 
    n_loss = lambda*n*sqrt(T) ;
  }
  else { // SOL_closure == "ESEL_like"
    n_loss = sigma_n*n ;
  }

  return n_loss ;
}


const Field3D STORM2D::loss_vort(const BoutReal &lambda, const Field3D &vort, const Field3D &UNUSED(n), const Field3D &phi, const Field3D &T){
// Parallel loss term for vorticity equation

  Field3D vort_loss ;
  vort_loss.allocate() ;

  if(SOL_closure == "sheath_diss") {
    // Sheath dissipation closure
    vort_loss = lambda*(V_sheath(phi, T) - sqrt(T)) ;
  }
  else if(SOL_closure == "vort_adv") { 
    // Vorticity advection closure
    vort_loss = lambda*sqrt(T)*vort ;
  }
  else { // SOL_closure == "ESEL_like"
    vort_loss = sigma_vort*vort ;
  }

  return vort_loss ;
}


const Field3D STORM2D::loss_T(const BoutReal &lambda, const Field3D &UNUSED(n), const Field3D &phi, const Field3D &T){
// Parallel loss term for temperature equation

  Field3D T_loss ;
  T_loss.allocate() ;

  if (SOL_closure == "sheath_diss") {
    // Sheath dissipation closure has no parallel conduction because grad_par(T) = 0
    T_loss = (2.0/3.0)*lambda*T*(1.71*V_sheath(phi, T) - 0.71*sqrt(T)) ;
  }
  else if (SOL_closure == "vort_adv"){
    // This is a positive quantity (temperature lost) because curvature of T is -ve at the mid-plane
    T_loss = (2.0/7.0)*(2.0/3.0)*kappa0*lambda*lambda*pow(T,3.5) + (2.0/3.0)*lambda*sqrt(T)*T ;
  }
  else {// (SOL_closure == "ESEL_like")
    T_loss = sigma_T*T ; 
  }
  
  return T_loss ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////Curvature operator for s-alpha geometry in infinite aspect ratio limit
const Field3D STORM2D::curv_op(const Field3D &f){
  return g0*DDZ(f);
}

BOUTMAIN(STORM2D);

