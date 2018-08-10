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
#include <invert_laplace.hxx>
#include <bout_types.hxx>

class STORM : public PhysicsModel{
public:
  virtual ~STORM() {}
protected:
  int init(bool restarting);
  int rhs(BoutReal t);
private:

  Options* globalOptions;
  Options* options;

  void initialise_background(bool restarting);
  void initialise_blob(const char * imp_section);
  FieldPerp extrap_sheath_lower(const Field3D &var);
  FieldPerp extrap_sheath_upper(const Field3D &var);
  void set_equilibrium_value(Field3D & f,const char * fname);

  // allocate and set the value of a 1d (y-direction) array
  // do nothing if it has already been allocated
  void set_equilibrium_value_1d_array(BoutReal* &f,const char * fname);

  // Set inner and outer boundary 1d arrays from a Field3D
  void set_equilibrium_value_1d_array(BoutReal* &f_inner, BoutReal* &f_outer, const Field3D &f3d);

  void read_equilibrium_file(BoutReal * data, const char * fname);

  Field3D n;                   // density
  Field3D U;                   // ion parallel velocity
  Field3D V;                   // electron parallel velocity
  Field3D vort;                // vorticity
  Field3D phi;                 // electrostatic potential
  Field3D S;                   // density source
  Field3D T;                   // electron temperature
  Field3D qpar;                // parallel heat flux
  Field3D S_E;                 // energy source

  FieldPerp phi_sheath_lower, phi_sheath_upper, n_sheath_lower, n_sheath_upper, T_sheath_lower, T_sheath_upper;
  FieldPerp V_sheath_lower, V_sheath_upper, U_sheath_upper, U_sheath_lower;
  Field3D nu_parallel, mu_n, mu_vort,kappa_perp;
  Field3D n_stag;
  Field3D V_centre;
  Field3D T_stag;
  Field3D Tcoef;
  Field3D phi_stag;
  Field3D S_stag;
  Field3D sqrtT; //sqrtT = sqrt(T)
  Field3D p; //p=n*T
  Field3D UmV; //UmV = U-V
  Field3D UmV_centre;

  Field3D Curv_n;               // Curv(n)
  Field3D Curv_T;               // Curv(T)
  Field3D Curv_p;               // Curv(p)
  Field3D Curv_phi;             // Curv(phi)

  // Fields and arrays used to initialise save background
  BoutReal *data, *phi_array_inner, *phi_array_outer;
  Field3D n_eq, U_eq, V_eq, phi_eq, T_eq, q_eq ;
  string equilibrium_source;    // specifies which method to use to read in or initialise initial profiles
  string equilibrium_file_path; // path to the directory where equilibrium files are stored
  string equilibrium_data_file; // netcdf or hdf5 file with initial profiles.
                                // If equilibrium_data_file is an empty string (the default value)
                                // then read 1d initial profiles from .dat binary files

  BoutReal mu_n0;               // density diffusion coefficient
  BoutReal mu_vort0;            // ion viscosity
  BoutReal mu;                  // m_i/m_e
  BoutReal nu_parallel0;        // electron-ion friction
  BoutReal g0;
  BoutReal kappa0;              // parallel heat conduction
  BoutReal kappa0_perp;
  BoutReal phi_wall;            // potential of target walls.

  Field3D n_on_sqrt_T;          // working field for nonuniform diss params
  Field3D Grad_par_phi, Grad_par_T;
  Field3D powT_1_5_stag;        // pow(T_stag,1.5)

  // Parameters used to set normalisations correctly.
  BoutReal u, e, m_e, epsilon_0, mu_0;
  BoutReal B_0, T_e0, T_i0, n_0, m_i, q, R_c;
  BoutReal Omega_i, Omega_e, c_s, rho_s;
  BoutReal nu_ei, nu_ii;
  BoutReal rho_e, rho_i;
  BoutReal V_the, V_thi;
  BoutReal loglambda;
  BoutReal Vsheath_BC_prefactor, Usheath_BC_prefactor;
  BoutReal Lx, Ly, Lz;

  bool isothermal;             // Switch for isothermal simulations
  bool uniform_diss_paras;     // Switch for global/local values of mu_vort, mu_n
  bool run_1d;                 // Run in 1d (parallel only) with phi=ln(n)+const.
  BoutReal run_1d_T_slowdown;  // When running in 1d, decrease ddt(T) by this factor to allow longer timesteps
  bool symmetry_plane;         // Switch for running in a half-domain, assuming reflection symmetry about the midpoint of the parallel direction
  bool add_blob;               // Switch to add a blob onto the equilibrium
  bool verbose;                // Add extra output

  // BoutReal *phi_x_boundary;
  FILE *file;
  FieldGroup comms;
  BRACKET_METHOD bm;           // Bracket method for advection terms
  Laplacian* phiSolver;        // Laplacian solver for vort -> phi

  //////////////////////////////////////////////////////////////
  // Methods for boundary conditions
  //////////////////////////////////////////////////////////////
  // Electron Velocity Sheath
  void Vsheath_yup_staggered(Field3D &var, const FieldPerp &phisheath, const BoutReal phi_wall, const FieldPerp &T_sheath, const BoutReal Vsheath_BC_prefactor);
  void Vsheath_ydown_staggered(Field3D &var, const FieldPerp &phisheath, const BoutReal phi_wall, const FieldPerp &T_sheath, const BoutReal Vsheath_BC_prefactor);

  // Ion Velocity Sheath
  void Usheath_yup_staggered(Field3D &U, const FieldPerp &sqrtT, const BoutReal Usheath_BC_prefactor);
  void Usheath_ydown_staggered(Field3D &U, const FieldPerp &sqrtT, const BoutReal Usheath_BC_prefactor);

  //Parallel heat flux
  void qsheath_yup_staggered(Field3D &var, const FieldPerp &T_sheath, const FieldPerp &n_sheath, const FieldPerp &U_sheath, const FieldPerp &V_sheath, const BoutReal mu);
  void qsheath_ydown_staggered(Field3D &var, const FieldPerp &T_sheath, const FieldPerp &n_sheath, const FieldPerp &U_sheath, const FieldPerp &V_sheath, const BoutReal mu);

  // Sets x guard cells to values from an array/another field's guard points.  Used for phi perpendicular boundary condition
  void set_xguards(Field3D &f, const BoutReal *g_inner, const BoutReal *g_outer);

  //////////////////////////////////////////////////////////////
  // STORM operators
  //////////////////////////////////////////////////////////////
  const Field3D Curv(const Field3D &f);

  //////////////////////////////////////////////////////////////
  // Utility methods
  //////////////////////////////////////////////////////////////

  void phisolver_1d();

  // Try to get Lx, Ly and Lz from the options. If they are not present, then
  // calculate from dx/dy/dz and metric, checking that all quantities used are
  // constant in order to avoid getting different answers on different
  // processors
  void set_Lx_Ly_Lz();
};

#endif // __ FILAMENTS_H__
