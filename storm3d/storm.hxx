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

#ifndef __FILAMENTS_H__
#define __FILAMENTS_H__

#include <bout/physicsmodel.hxx>
#include <interpolation.hxx>
#include <invert_laplace.hxx>
#include <bout/invert/laplacexy.hxx>
#include <bout_types.hxx>
#include <interpolation.hxx>
#include "../shared/fast_output.hxx"
#include "neutral-model.hxx"

class STORM : public PhysicsModel{
public:
  ~STORM() {}
protected:
  int init(bool restarting);
  int rhs(BoutReal t);
  int timestepMonitor(BoutReal simtime, BoutReal UNUSED(dt)) {
    int ret = 0;
    if(evolving_bcs) compute_bndry_phi(simtime);
    if(fast_output.enable_timestep) {
      ret = fast_output.monitor_method(simtime);
      if (ret) return ret; // return immediately if ret is non-zero (indicating error)
    }
    if( monitor_minmaxmean ) {
      if(minmaxmean_timelast < 0.){
        minmaxmean_timelast = simtime;
      }else if(simtime - minmaxmean_timelast > 5.){
        minmaxmean_timelast = simtime;
        output << "t = " << simtime << endl;
        output.write("\nmin(phi) = %e, max(phi) = %e, mean(phi) = %e\n",
          min(phi,true,"RGN_NOBNDRY"),max(phi,true,"RGN_NOBNDRY"),mean(phi,true,"RGN_NOBNDRY"));
        output.write("min(n) = %e, max(n) = %e, mean(n) = %e\n",
          min(n,true,"RGN_NOBNDRY"),max(n,true,"RGN_NOBNDRY"),mean(n,true,"RGN_NOBNDRY"));
        output.write("min(T) = %e, max(T) = %e, mean(T) = %e\n",
          min(T,true,"RGN_NOBNDRY"),max(T,true,"RGN_NOBNDRY"),mean(T,true,"RGN_NOBNDRY"));
        output.write("min(U) = %e, max(U) = %e, mean(U) = %e\n",
          min(U_aligned,true,"RGN_NOBNDRY"),max(U_aligned,true,"RGN_NOBNDRY"),mean(U_aligned,true,"RGN_NOBNDRY"));
        output.write("min(V) = %e, max(V) = %e, mean(V) = %e\n",
          min(V_aligned,true,"RGN_NOBNDRY"),max(V_aligned,true,"RGN_NOBNDRY"),mean(V_aligned,true,"RGN_NOBNDRY"));
        output.write("min(vort) = %e, max(vort) = %e, mean(vort) = %e\n",
          min(vort,true,"RGN_NOBNDRY"),max(vort,true,"RGN_NOBNDRY"),mean(vort,true,"RGN_NOBNDRY"));
      }
    }
    return ret;
  }
  int outputMonitor(BoutReal simtime, int iteration, int nout);

  int precon(BoutReal t, BoutReal gamma, BoutReal delta);
private:

  Options* globalOptions;
  Options* options;
  Coordinates* coordinates_centre;
  Coordinates* coordinates_stag;

  //////////////////////////////////////////////////////////////
  // Initialisation methods
  //////////////////////////////////////////////////////////////
  void initialise_background();
  void initialise_blob(const char * imp_section) ; 
  void add_noise() ;
  void add_perturbation(Field3D &f, const std::string name);
  void setBoundaryConditionsOptions();
  void check_U_V_x_boundary_conditions();

  void phisolver_1d() ; 
  FieldPerp extrap_sheath_lower(const Field3D &var);
  FieldPerp extrap_sheath_upper(const Field3D &var);
  void set_equilibrium_value(Field3D & f,const char * fname);

  // allocate and set the value of a 1d (y-direction) array
  // do nothing if it has already been allocated
  void set_equilibrium_value_1d_array(BoutReal* &f,const char * fname);

  // Set inner and outer boundary 1d arrays from a Field3D
  void set_equilibrium_value_1d_array(BoutReal* &f_inner, BoutReal* &f_outer, const Field3D &f3d);

  void read_equilibrium_file(BoutReal * data, const char * fname);

  Field3D logn;
  Field3D n;                   // density
  Field3D U;                   // ion parallel velocity
  Field3D V;                   // electron parallel velocity
  Field3D vort;                // vorticity
  Field3D phi;                 // electrostatic potential
  Field2D phi2D;               // axi-symmetric component of phi
  Field3D S;                   // density source
  Field3D logp; 
  Field3D logT; 
  Field3D T;                   // electron temperature
  Field3D qpar_aligned;        // parallel heat flux
  Field3D uE2;                 // square of the absolute value of the perpendicular gradient of phi (used when Boussinesq is NOT used)
  Field3D S_E;                 // energy source

  Field3D psi;                 // parallel component of the magnetic vector potential
  Field3D chiU;                // U+beta/2*psi
  Field3D chiV;                // V-mu*beta/2*psi

  FieldPerp phi_sheath_lower, phi_sheath_upper, n_sheath_lower, n_sheath_upper, T_sheath_lower, T_sheath_upper;
  FieldPerp psi_sheath_lower, psi_sheath_upper;
  FieldPerp V_sheath_lower, V_sheath_upper, U_sheath_upper, U_sheath_lower ;
  Field3D nu_parallel, mu_n, mu_vort,kappa_perp;
  Field3D logn_stag;
  Field3D n_stag ;
  Field3D U_centre;
  Field3D V_centre;
  Field3D logT_stag;
  Field3D T_stag;
  Field3D logn_aligned, logT_aligned;
  Field3D n_aligned, U_aligned, V_aligned, vort_aligned, phi_aligned, T_aligned, psi_aligned;
  Field3D qpar_centre;
  Field3D Tcoef;
  Field3D phi_stag;
  Field3D S_stag;
  Field3D sqrtT; //sqrtT = sqrt(T)
  Field3D p; //p=n*T
  Field3D UmV; //UmV = U-V
  Field3D UmV_centre;
  Field3D psi_centre;

  Field3D Curv_n;   // Curv(n)
  Field3D Curv_T;   // Curv(T)
  Field3D Curv_p;   // Curv(p)
  Field3D Curv_phi; // Curv(phi)

  // Fields and arrays used to initialise save background
  //BoutReal *n_array, *phi_array, *U_array, *V_array ; 
  BoutReal *data, *phi_array_inner, *phi_array_outer;
  Field3D n_eq, U_eq, V_eq, phi_eq, T_eq, q_eq ;
  std::string equilibrium_source;    // specifies which method to use to read in or initialise initial profiles
  std::string equilibrium_file_path; // path to the directory where equilibrium files are stored
  std::string equilibrium_data_file; // netcdf or hdf5 file with initial profiles.
                                // If equilibrium_data_file is an empty string (the default value)
                                // then read 1d initial profiles from .dat binary files

  BoutReal mu_n0 ;              // density diffusion coefficient
  BoutReal mu_vort0 ;           // ion viscosity
  BoutReal mu ;                 // m_i/m_e
  BoutReal nu_parallel0 ;       // electron-ion friction
  BoutReal g0 ; 
  BoutReal kappa0 ;	            // parallel heat conduction
  BoutReal kappa0_perp ;        
  BoutReal diff_perp_U;
  BoutReal diff_perp_V;
  BoutReal phi_wall ;           // potential of target walls. 
  BoutReal beta0 ;              // ratio of plasma to magnetic pressure (with normalization parameters)

  Field3D n_on_sqrt_T ;         // working field for nonuniform diss params
  Field3D Grad_par_phi_stag, Grad_par_T_stag;
  Field3D powT_1_5_stag;        // pow(T_stag,1.5)

  // Parameters used to set normalisations correctly.  
  BoutReal u, e, m_e, epsilon_0, mu_0, Z ;               
  BoutReal B_0, T_e0, T_i0, n_0, m_i, q, R_c ; 
  BoutReal Omega_i, Omega_e, c_s, rho_s ; 
  BoutReal nu_ei, nu_ii ; 
  BoutReal rho_e, rho_i ; 
  BoutReal V_the, V_thi ; 
  BoutReal loglambda ; 
  BoutReal Vsheath_BC_prefactor, Usheath_BC_prefactor ;
  BoutReal Lx, Ly, Lz;

  int ixseps1 ; 
  int ixseps2 ;
  int jyseps1_1, jyseps2_1, jyseps1_2, jyseps2_2, ny_inner;

  bool isothermal ;            // Switch for isothermal simulations
  int boussinesq ;            // Switch for the Boussinesq approximation
  bool split_n0 ;             // Flag to split the n=0 component in the Poisson equation
  bool electromagnetic;        // Switch for evolution of parallel magnetic vector potential, 'electromagnetic effects'
  bool uniform_diss_paras ;    // Switch for global/local values of mu_vort, mu_n
  bool old_phi_wall_value ;    // Switch for setting phi_wall = 0.5*ln(mu/TWOPI) 
  bool S_in_peq ;              // Switch for 0.5*V^2*S/n/mu in pressure equation
  bool run_1d;                 // Run in 1d (parallel only) with phi=ln(n)+const.
  bool hydrodynamic;           // Run with vort=0 and Jpar=n*(U-V)=0. Solve for phi from parallel Ohm's law with no inertia or resistivity
  BoutReal run_1d_T_slowdown;  // When running in 1d, decrease ddt(T) by this factor to allow longer timesteps
  bool symmetry_plane;         // Switch for running in a half-domain, assuming reflection symmetry about the midpoint of the parallel direction
  bool normalise_lengths ;     // Switch to renormalise lengths
  bool add_blob ;              // Switch to add a blob onto the equilibrium
  bool two_blobs ;             // Adds a second blob specified in the section "blob_1" of BOUT.inp
  bool initial_noise ;         // Switch to add random field onto the equilibrium
  bool initial_perturbation ;  // Switch to add perturbations specified in input file onto the equilibrium
  bool verbose ;	       // Add extra output
  bool evolving_bcs;           // Switch to activate averaged Neumann boundary conditions for phi
  bool use_psi_boundary_solver;// Use a Laplace solver to invert for psi on the boundary points, rather than extrapolating
  bool save_aligned_fields;    // Switch to save field-aligned versions of fields to dump files
  bool average_radial_boundaries_core_SOL; // Ad hoc boundary conditions in radial direction
  std::string realistic_geometry; // none for slab, salpha or doublenull
  bool normalise_all;          // flag for normalizing magnetic field, dx, metric coefficients, and curvature operator
  bool monitor_minmaxmean;     // flag to enable the output of the min, max and average of the fields
  bool sources_realistic_geometry;    // flag for sources for realistic geometries
  bool sources_realisticgeometry_background; // flag for background (useful for initial phase)
  bool isshifted;              // chech if shifted parallel transform
  bool increased_dissipation_xbndries; // flag to increase the perpendicular dissipation near the inner and outer boundaries
  bool increased_resistivity_core; // flag to increase resistivity in core region
  bool normalise_sources;      // flag for normalise sources to Ly

  Vector2D bxcv, Curlb_B;      // Vectors for realistic curvature operator
  Field2D B2;                  // Bxy*Bxy
  Field2D G3;                  // Useful to temporary overwrite G3

  // BoutReal *phi_x_boundary ; 
  FILE *file ; 
  FieldGroup comms ;
  BRACKET_METHOD bm;           // Bracket method for advection terms
  Laplacian* phiSolver;  // Laplacian solver for vort -> phi
  Laplacian* psiSolver;  // Laplacian solver for {chiU,chiV} -> psi
  Laplacian* psiBoundarySolver;  // Laplacian solver for n*(U-V) -> psi at y-boundaries, used if use_psi_boundary_solver==true
  LaplaceXY* phiSolverxy; // Laplacian for the axi-symmetric component of phi

  // Neutral gas model
  NeutralModel *neutrals; // Handles evolution of neutral gas
  FieldPerp ionflux_lower, ionflux_upper; // Ion fluxes to the target plates

  //////////////////////////////////////////////////////////////
  // Methods for boundary conditions
  //////////////////////////////////////////////////////////////
  // Electron Velocity Sheath
  void Vsheath_yup_staggered(Field3D &var, const FieldPerp &phisheath, const BoutReal phi_wall, const FieldPerp &T_sheath, const BoutReal Vsheath_BC_prefactor) ;
  void Vsheath_ydown_staggered(Field3D &var, const FieldPerp &phisheath, const BoutReal phi_wall, const FieldPerp &T_sheath, const BoutReal Vsheath_BC_prefactor) ;

  // Ion Velocity Sheath
  void Usheath_yup_staggered(Field3D &U, const FieldPerp &sqrtT, const BoutReal Usheath_BC_prefactor);
  void Usheath_ydown_staggered(Field3D &U, const FieldPerp &sqrtT, const BoutReal Usheath_BC_prefactor);

  //Parallel heat flux
  void qsheath_yup_staggered(Field3D &var, const FieldPerp &T_sheath, const FieldPerp &n_sheath, const FieldPerp &V_sheath, const BoutReal mu) ; 
  void qsheath_ydown_staggered(Field3D &var, const FieldPerp &T_sheath, const FieldPerp &n_sheath, const FieldPerp &V_sheath, const BoutReal mu) ; 

  // Sets x guard cells to values from an array/another field's guard points.  Used for phi perpendicular boundary condition
  void set_xguards(Field3D &f, const BoutReal *g_inner, const BoutReal *g_outer) ;
  void set_xguards(Field3D &f, const Field3D &g) ; 

  // Apply the boundary conditions on psi
  void applyPsiBoundaries();

  // Variables and functions used for averaged Neumann boundary conditions
  void phi_bc_initialise(bool restarting);
  void compute_bndry_phi(const BoutReal time);
  void apply_bndry_phi();
  BoutReal time_last_SOL, time_last_PF, time_last_core;
  BoutReal cstep_SOL, cstep_PF, cstep_core;
  BoutReal time_update_SOL, time_update_PF, time_update_core;
  Field2D phi_bc;
  bool first_step;

  // Internal variable for monitor
  BoutReal minmaxmean_timelast = -1.;

  // Apply neumann boundary condition on the DC (Z) and <->_YZ component
  void average_Z_bndry(Field3D& f, const bool left = true, const bool right = true, const bool periodic = false);
  void average_YZ_bndry(Field3D& f, const bool left = true, const bool right = true, const bool periodic = false);

  //////////////////////////////////////////////////////////////
  // Object to handle fast output
  //////////////////////////////////////////////////////////////
  FastOutput fast_output;

  //////////////////////////////////////////////////////////////
  // STORM operators
  //////////////////////////////////////////////////////////////
  const Field3D Curv(const Field3D& f, const Field3D& f_aligned);
  const Field3D Grad_par_EM(const Field3D &var_aligned, const Field3D &var_outloc,
      const Field3D &Psi, CELL_LOC outloc=CELL_DEFAULT,
      const std::string& method="DEFAULT");
  const Field3D Div_par_EM(const Field3D &var_aligned, const Field3D &var_outloc,
      const Field3D &Psi, CELL_LOC outloc=CELL_DEFAULT,
      const std::string& method="DEFAULT");
  const Field3D Vpar_Grad_par_EM(const Field3D& v_aligned, const Field3D& f_aligned,
      const Field3D& v_outloc, const Field3D& f, const Field3D& Psi,
      CELL_LOC outloc=CELL_DEFAULT, const std::string& method="DEFAULT");
  const Field3D interp_to_fixstag(const Field3D& f, CELL_LOC loc);

  // Grad_perp_dot_Grad_perp implemented without any y-derivatives (eventually
  // including y-derivatives should be optionally included). Avoids
  // intermediate Vector3Ds so is slightly more efficient than
  // Grad_perp(f)*Grad_perp(g) even ignoring the cost of y-derivatives.
  const Field3D Grad_perp_dot_Grad_perp(const Field3D& f, const Field3D& g);
  // Variant to calculate Grad_perp(f).Grad_perp(f)
  const Field3D Grad_perp_dot_Grad_perp(const Field3D& f);

  // Wrapper for interp_to that returns the result in non-field-aligned form
  const Field3D stormInterp_to(const Field3D& f, CELL_LOC outloc) {
    return fromFieldAligned(interp_to(f, outloc, "RGN_NOBNDRY"), "RGN_NOBNDRY");
  }

  //////////////////////////////////////////////////////////////
  // Utility methods
  //////////////////////////////////////////////////////////////

  // Try to get Lx, Ly and Lz from the options. If they are not present, then
  // calculate from dx/dy/dz and metric, checking that all quantities used are
  // constant in order to avoid getting different answers on different
  // processors
  void set_Lx_Ly_Lz();
  Field3D vort_from_phi(const Field3D phi);

  // These two operators can be used to calculate the exact inverse of the
  // multigrid Laplacian operator. Their implementation is copied from BOUT++,
  // tests/integrated/test-multigrid-laplace/test_multigrid_laplace.cxx
  Field3D this_Grad_perp2(const Field3D &f);
  Field3D this_Grad_perp_dot_Grad_perp(const Field3D &f, const Field3D &g);

  // Method to set the ddt guard cells of a variable to 0 when on staggered grid
  void set_lower_ddt_zero(Field3D &var);
  
  // Method to set the sources for realistic geometry simulations
  void set_sources_realistic_geometry();

};

#endif // __ FILAMENTS_H__
