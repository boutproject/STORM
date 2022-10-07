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

#include <initialprofiles.hxx>
#include <interpolation.hxx>
#include <field_factory.hxx>
#include <bout/constants.hxx>
#include <dataformat.hxx>
#include <derivs.hxx>

////////////////////////////////////////////////////////////////////////

void STORM::initialise_background() {

  // Always load values, even if restarting, so *_eq fields get set.
  // When restarting, restart values are loaded into fields after init() is done.

  if (equilibrium_source == "1d_profiles") {
    if (mesh->firstX()) {
      set_equilibrium_value_1d_array(phi_array_inner, "phi_eq.dat");
    }
    if (mesh->lastX()) {
      set_equilibrium_value_1d_array(phi_array_outer, "phi_eq.dat");
    }
    
    set_equilibrium_value(n  ,"n_eq.dat");
    set_equilibrium_value(U  ,"U_eq.dat");
    set_equilibrium_value(V  ,"V_eq.dat");
    set_equilibrium_value(phi,"phi_eq.dat");

    vort = 0.0 ;
    if (electromagnetic) psi = 0.;

    if (!isothermal) {
      set_equilibrium_value(T,"T_eq.dat");
    }
  }
  else if (equilibrium_source == "profiles_file") {
    if (equilibrium_data_file == "")
      throw BoutException("Error: no equilibrium_data_file specified. This is required when equilibrium_source=\"profiles_file\"");
    GridFile equilib_profiles(data_format(equilibrium_data_file.c_str()), equilibrium_file_path + "/" + equilibrium_data_file.c_str());
    equilib_profiles.get(mesh, phi, "phi");
    set_equilibrium_value_1d_array(phi_array_inner, phi_array_outer, phi);

    equilib_profiles.get(mesh, n, "n");
    equilib_profiles.get(mesh, U, "U");
    equilib_profiles.get(mesh, V, "V");
    equilib_profiles.get(mesh, vort, "vort");
    if (!isothermal) {
      equilib_profiles.get(mesh, T, "T");
    }
    if (electromagnetic) {
      equilib_profiles.get(mesh, psi, "psi");
    }
  }
  else if (equilibrium_source == "input_file") {
    // This is what BOUT++ does by default, so just don't need to initialise Field3D's here
    // except for phi, since phi is not an evolving field
    if(mesh->get(phi,"phi")) { // Try to read from grid file
      initial_profile("phi",phi); // Get initial profile from options in input file
    }
    // U and V are auxiliary variables in EM version, need to initialize explicitly
    initial_profile("U", U);
    initial_profile("V", V);
    set_equilibrium_value_1d_array(phi_array_inner, phi_array_outer, phi);

    if (electromagnetic) {
      initial_profile("psi", psi);
    }
  }
  else if (equilibrium_source == "single null" || equilibrium_source == "double null") {
    if (realistic_geometry != "singlenull" && realistic_geometry != "doublenull")
      throw BoutException("Need realistic geometry == singlenull or doublenull.");
    n = 1.;
    T = 1.;
    U = 0.;
    // y goes from 0 to 1 in the topological connected region
    BoutReal dy, y;
    for(int ix=mesh->xstart; ix<=mesh->xend; ++ix){
      if (mesh->getGlobalXIndex(ix) < ixseps1) {
        if (mesh->getGlobalYIndexNoBoundaries(mesh->yend) <= jyseps1_1) {
          //PF region in between 0 and jyseps1_1
          for(int iy=mesh->ystart; iy<=mesh->yend; ++iy){
            for(int iz=0; iz<mesh->LocalNz; ++iz){
              dy = 0.5/((double)(jyseps1_1 + 1));
              y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy)) + 0.5)*dy;
              n(ix,iy,iz) = -1.2*y*y + 1.2*y + 0.35;
              T(ix,iy,iz) = -0.23*y*y+0.23*y+0.5;
              U(ix,iy,iz) = 5.*y*y*y - 1.5*5.*y*y + (2.*sqrt(T(ix,iy,iz))+5./2.)*y - sqrt(T(ix,iy,iz));
            }
          }
        }
        else if (mesh->getGlobalYIndexNoBoundaries(mesh->ystart) > jyseps1_1 && mesh->getGlobalYIndexNoBoundaries(mesh->yend) <= jyseps2_1) {
          // Core
          for(int iy=mesh->ystart; iy<=mesh->yend; ++iy){
            for(int iz=0; iz<mesh->LocalNz; ++iz){
              n(ix,iy,iz) = -1.2*0.5*0.5+1.2*0.5+0.35;
              T(ix,iy,iz) = -0.23*0.5*0.5+0.23*0.5+0.5;
            }
          }
        }
        else if (mesh->getGlobalYIndexNoBoundaries(mesh->ystart) > jyseps2_1 && mesh->getGlobalYIndexNoBoundaries(mesh->yend) < ny_inner) {
          // PF region between jyseps2_1 and ny_inner
          for(int iy=mesh->ystart; iy<=mesh->yend; ++iy){
            for(int iz=0; iz<mesh->LocalNz; ++iz){
              dy = 0.5/((double)(ny_inner - jyseps2_1 - 1));
              y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy) - jyseps2_1) - 0.5)*dy + 0.5;
              n(ix,iy,iz) = -1.2*y*y + 1.2*y + 0.35;
              T(ix,iy,iz) = -0.23*y*y+0.23*y+0.5;
              U(ix,iy,iz) = 5.*y*y*y - 1.5*5.*y*y + (2.*sqrt(T(ix,iy,iz))+5./2.)*y - sqrt(T(ix,iy,iz));
            }
          }
        }
        else if (mesh->getGlobalYIndexNoBoundaries(mesh->ystart) >= ny_inner &&  mesh->getGlobalYIndexNoBoundaries(mesh->yend) <= jyseps1_2) {
          // PF region between ny_inner and jyseps1_2
          for(int iy=mesh->ystart; iy<=mesh->yend; ++iy){
            for(int iz=0; iz<mesh->LocalNz; ++iz){
              dy = 0.5/((double)(jyseps1_2 - ny_inner + 1));
              y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy)) - ny_inner + 0.5)*dy;
              n(ix,iy,iz) = -1.2*y*y + 1.2*y + 0.35;
              T(ix,iy,iz) = -0.23*y*y+0.23*y+0.5;
              U(ix,iy,iz) = 5.*y*y*y - 1.5*5.*y*y + (2.*sqrt(T(ix,iy,iz))+5./2.)*y - sqrt(T(ix,iy,iz));
            }
          }
        }
        else if (mesh->getGlobalYIndexNoBoundaries(mesh->ystart) > jyseps1_2-1 && mesh->getGlobalYIndexNoBoundaries(mesh->yend) <= jyseps2_2) {
          // Core
          for(int iy=mesh->ystart; iy<=mesh->yend; ++iy){
            for(int iz=0; iz<mesh->LocalNz; ++iz){
              n(ix,iy,iz) = -1.2*0.5*0.5+1.2*0.5+0.35;
              T(ix,iy,iz) = -0.23*0.5*0.5+0.23*0.5+0.5;
            }
          }
        }
        else if (mesh->getGlobalYIndexNoBoundaries(mesh->ystart) > jyseps2_2) {
          // PF region between jyseps2_2 and GlobalNy
          for(int iy=mesh->ystart; iy<=mesh->yend; ++iy){
            for(int iz=0; iz<mesh->LocalNz; ++iz){
              dy = 0.5/((double)(mesh->GlobalNy - jyseps2_2 - 1));
              y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy) - jyseps2_2) - 0.5)*dy + 0.5;
              n(ix,iy,iz) = -1.2*y*y + 1.2*y + 0.35;
              T(ix,iy,iz) = -0.23*y*y+0.23*y+0.5;
              U(ix,iy,iz) = 5.*y*y*y - 1.5*5.*y*y + (2.*sqrt(T(ix,iy,iz))+5./2.)*y - sqrt(T(ix,iy,iz));
            }
          }
        }
      }
      else if (mesh->getGlobalXIndex(ix) >= ixseps1 && mesh->getGlobalXIndex(ix) < ixseps2) {
        if (mesh->getGlobalYIndexNoBoundaries(mesh->yend) <= jyseps2_1) {
          for(int iy=mesh->ystart; iy<=mesh->yend; ++iy){
            for(int iz=0; iz<mesh->LocalNz; ++iz){
              dy = 0.5/((double)(jyseps2_1 + 1));
              y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy)) + 0.5)*dy;
              n(ix,iy,iz) = -1.2*y*y + 1.2*y + 0.35;
              T(ix,iy,iz) = -0.23*y*y+0.23*y+0.5;
              U(ix,iy,iz) = 5.*y*y*y - 1.5*5.*y*y + (2.*sqrt(T(ix,iy,iz))+5./2.)*y - sqrt(T(ix,iy,iz));
            }
          }
        }
        else if (mesh->getGlobalYIndexNoBoundaries(mesh->ystart) > jyseps2_1 && mesh->getGlobalYIndexNoBoundaries(mesh->yend) < ny_inner) {
          for(int iy=mesh->ystart; iy<=mesh->yend; ++iy){
            for(int iz=0; iz<mesh->LocalNz; ++iz){
              dy = 0.5/((double)(ny_inner - jyseps2_1 - 1));
              y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy) - jyseps2_1) - 0.5)*dy + 0.5;
              n(ix,iy,iz) = -1.2*y*y + 1.2*y + 0.35;
              T(ix,iy,iz) = -0.23*y*y+0.23*y+0.5;
              U(ix,iy,iz) = 5.*y*y*y - 1.5*5.*y*y + (2.*sqrt(T(ix,iy,iz))+5./2.)*y - sqrt(T(ix,iy,iz));
            }
          }
        }
        else if (mesh->getGlobalYIndexNoBoundaries(mesh->ystart) >= ny_inner &&  mesh->getGlobalYIndexNoBoundaries(mesh->yend) <= jyseps1_2) {
          for(int iy=mesh->ystart; iy<=mesh->yend; ++iy){
            for(int iz=0; iz<mesh->LocalNz; ++iz){
              dy = 0.5/((double)(jyseps1_2 - ny_inner + 1));
              y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy)) - ny_inner + 0.5)*dy;
              n(ix,iy,iz) = -1.2*y*y + 1.2*y + 0.35;
              T(ix,iy,iz) = -0.23*y*y+0.23*y+0.5;
              U(ix,iy,iz) = 5.*y*y*y - 1.5*5.*y*y + (2.*sqrt(T(ix,iy,iz))+5./2.)*y - sqrt(T(ix,iy,iz));
            }
          }
        }
        else if (mesh->getGlobalYIndexNoBoundaries(mesh->ystart) > jyseps1_2-1) {
          for(int iy=mesh->ystart; iy<=mesh->yend; ++iy){
            for(int iz=0; iz<mesh->LocalNz; ++iz){
              dy = 0.5/((double)(mesh->GlobalNy - jyseps1_2 - 1));
              y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy) - jyseps1_2) - 0.5)*dy + 0.5;
              n(ix,iy,iz) = -1.2*y*y + 1.2*y + 0.35;
              T(ix,iy,iz) = -0.23*y*y+0.23*y+0.5;
              U(ix,iy,iz) = 5.*y*y*y - 1.5*5.*y*y + (2.*sqrt(T(ix,iy,iz))+5./2.)*y - sqrt(T(ix,iy,iz));
            }
          }
        }
      }
      else if (mesh->getGlobalXIndex(ix) >= ixseps2) {
        if(mesh->getGlobalYIndexNoBoundaries(mesh->yend) < ny_inner){
          for(int iy=mesh->ystart; iy<=mesh->yend; ++iy){
            for(int iz=0; iz<mesh->LocalNz; ++iz){
              dy = 1./((double)(ny_inner));
              y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy)) + 0.5)*dy;
              n(ix,iy,iz) = -1.2*y*y + 1.2*y + 0.35;
              T(ix,iy,iz) = -0.23*y*y+0.23*y+0.5;
              U(ix,iy,iz) = 5.*y*y*y - 1.5*5.*y*y + (2.*sqrt(T(ix,iy,iz))+5./2.)*y - sqrt(T(ix,iy,iz));
            }
          }
        }
        else{
          for(int iy=mesh->ystart; iy<=mesh->yend; ++iy){
            for(int iz=0; iz<mesh->LocalNz; ++iz){
              dy = 1./((double)(mesh->GlobalNy - ny_inner));
              y = ((double)(mesh->getGlobalYIndexNoBoundaries(iy)) - ny_inner + 0.5)*dy;
              n(ix,iy,iz) = -1.2*y*y + 1.2*y + 0.35;
              T(ix,iy,iz) = -0.23*y*y+0.23*y+0.5;
              U(ix,iy,iz) = 5.*y*y*y - 1.5*5.*y*y + (2.*sqrt(T(ix,iy,iz))+5./2.)*y - sqrt(T(ix,iy,iz));
            }
          }
        }
      }
    }
    if (isothermal) T = 1.;
    mesh->communicate(n,U,T);
    average_Z_bndry(n);
    if (!isothermal) {
      average_Z_bndry(T);
    }
    average_Z_bndry(U);
    V = copy(U);
    phi = log(Vsheath_BC_prefactor)*T;
    coordinates_centre->G3 = 0.;
    if (boussinesq > 0) {
      vort = Delp2(phi)/B2 + coordinates_centre->g11*DDX(1./B2)*DDX(phi);
    }
    else {
      vort = n*Delp2(phi)/B2 + coordinates_centre->g11*DDX(n/B2)*DDX(phi);
    }
    coordinates_centre->G3 = G3;
    average_Z_bndry(vort);
    mesh->communicate(vort);
    if (electromagnetic) psi = 0.;
  }
  else {
    throw BoutException("Unrecognized equilibrium_source option '%s'", equilibrium_source.c_str());
  }

  n_eq = n ;
  U_eq = U ;
  V_eq = V ;
  phi_eq = phi ;
  SAVE_ONCE(phi_eq);
  SAVE_ONCE(n_eq);
  SAVE_ONCE(U_eq);
  SAVE_ONCE(V_eq);
  
  if(!isothermal){
    T_eq = T;
    //q_eq = qpar;
    SAVE_ONCE(T_eq);
    //SAVE_ONCE(q_eq);
  }
    
}

////////////////////////////////////////////////////////////////////////

void STORM::set_equilibrium_value(Field3D & f, const char * fname){
  data = new BoutReal[mesh->GlobalNy];
  read_equilibrium_file(data,fname);

  f.allocate();
  for(int i=0; i<mesh->LocalNx; i++){
    if(mesh->getGlobalXIndex(i)>=ixseps1){
      for(int j=0;j<mesh->LocalNy;j++){
        int jglobal = mesh->getGlobalYIndex(j);
        for(int k=0;k<mesh->LocalNz;k++){
          f(i,j,k) = data[jglobal];
        }
      }
    }
  }
}

void STORM::set_equilibrium_value_1d_array(BoutReal* &f, const char * fname){
  // if f has already been set, then do nothing here
  if (f == nullptr) {
    f = new BoutReal[mesh->GlobalNy];
    read_equilibrium_file(f,fname);
  }
}

void STORM::set_equilibrium_value_1d_array(BoutReal* &f_inner, BoutReal* &f_outer, const Field3D &f3d) {
  if (mesh->firstX()) {
    // if f_inner has already been set, then do nothing here
    if (f_inner == nullptr) {
      // Add second set of guard cells to mesh->GlobalNy in case this is used in a
      // double-null simulation. This can be removed when BOUT++ v5 is released, because
      // then GlobalNy will redefined to include all boundary cells; then it will always
      // be true that mesh->getGlobalYIndex(mesh->LocalNy) <= mesh->GlobalNy.
      f_inner = new BoutReal[mesh->GlobalNy + 2*mesh->ystart];
#if CHECK>2
      for (int i = 0; i < mesh->GlobalNy; i++) {
        f_inner[i] = nan("");
      }
#endif
      for(int j = 0; j < mesh->LocalNy; j++){
        int jglobal = mesh->getGlobalYIndex(j);
        f_inner[jglobal] = 0.5*(f3d(mesh->xstart-1, j, 0) + f3d(mesh->xstart, j, 0)) ;
      }
    }
  }

  if (mesh->lastX()) {
    // if f_outer has already been set, then do nothing here
    if (f_outer == nullptr) {
      // Add second set of guard cells to mesh->GlobalNy in case this is used in a
      // double-null simulation. This can be removed when BOUT++ v5 is released, because
      // then GlobalNy will redefined to include all boundary cells; then it will always
      // be true that mesh->getGlobalYIndex(mesh->LocalNy) <= mesh->GlobalNy.
      f_outer = new BoutReal[mesh->GlobalNy + 2*mesh->ystart];
#if CHECK>2
      for (int i = 0; i < mesh->GlobalNy; i++) {
        f_outer[i] = nan("");
      }
#endif
      for(int j = 0; j < mesh->LocalNy; j++){
        int jglobal = mesh->getGlobalYIndex(j);
        f_outer[jglobal] = 0.5*(f3d(mesh->xend, j, 0) + f3d(mesh->xend+1, j, 0)) ;
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////

void STORM::read_equilibrium_file(BoutReal * data, const char * fname){
  std::string datadir;
  int max=512;
  char * fullname = new char[max];
  int size = snprintf(fullname,max,"%s/%s",equilibrium_file_path.c_str(),fname);
  if (size > max){
    delete [] fullname;
    fullname = new char[size+20];
    snprintf(fullname,size+10,"%s/%s",equilibrium_file_path.c_str(),fname);
  }
  
  FILE * file = fopen(fullname,"rb");
  delete [] fullname;
  if (file == NULL){
    throw BoutException("Background file `%s` not found.\n",fname);
  } else {
    int nread = fread(data, sizeof(BoutReal), mesh->GlobalNy, file);
    if (nread != mesh->GlobalNy) {
      // something went wrong
      if (feof(file)) {
        // reached end of file too soon
        throw BoutException("Reached end of file while trying to read %s", fname);
      }
      else if (ferror(file)) {
        // there was some file I/O error
        throw BoutException("Error while reading %s", fname);
      }
      else {
        // something else went wrong, this should not happen
        throw BoutException("Unexpected error while reading %s", fname);
      }
    } else {
      // check we reached the end of the file
      // first try to read another byte from the file
      fgetc(file);
      // now the end-of-file indicator should be set
      if (!feof(file)) {
        throw BoutException("The input data in %s was longer than mesh->GlobalNy", fname);
      }
    }
    fclose(file);
  }
}

////////////////////////////////////////////////////////////////////////

void STORM::initialise_blob(const char * imp_section){
  // Parameters Used for Blob initialisation
  bool boltzmann ;             // Switch for initialising potential of blob to boltzmann response 
  bool conserve_momentum ;     // Switch for adjusting V when seeding blob such that nV is conserved
  BoutReal delta_perp = -1.;    // Perpendicular (x) length scale of density blob
  BoutReal delta_perp_T = -1.;  // Perpendicular (x) length scale of temperature blob
  BoutReal elongation = -1.;    // Elongation of the density blob: ratio between z and x axis
  BoutReal elongation_T = -1.;  // Elongation of the temperature blob: ratio between z and x axis
  BoutReal A = -1.;             // Density blob amplitude
  BoutReal A_T = -1.;           // Temperature blob amplitude
  BoutReal L_b = -1.;           // Parallel extent of density blob
  BoutReal L_b_T = -1.;         // Parallel extent of temperature blob
  BoutReal xoffset = -1.;       // starting z position of density blob (between 0 and 1)
  BoutReal xoffset_T = -1.;     // starting z position of temperature blob (between 0 and 1)
  BoutReal zoffset = -1.;       // starting z position of density blob (between 0 and 1)
  BoutReal zoffset_T = -1.;     // starting z position of temperature blob (between 0 and 1)
  BoutReal angle_blob = -1.;    // tilt of the density blob with respect to vertical axis, in rad
  BoutReal angle_blob_T = -1.;  // tilt of the temperature blob with respect to vertical axis, in rad
  BoutReal delta_front = -1.;   // Parameter used to control gradient of the density blob front. 0 = step function
  BoutReal delta_front_T = -1.; // Parameter used to control gradient of the temperature blob front. 0 = step function
  bool A_relative_to_bg;        // Switch for controlling whether the amplitude of the blob is relative to the midplane density or not.
  Options *blob_options = Options::getRoot()->getSection(imp_section);
  
  OPTION(blob_options, delta_perp,          10) ;
  OPTION(blob_options, elongation,           1) ;
  OPTION(blob_options, angle_blob,           0) ;
  OPTION(blob_options, A,                  2.0) ;
  OPTION(blob_options, L_b,                0.5) ;
  OPTION(blob_options, xoffset,           0.25) ;
  OPTION(blob_options, zoffset,            0.5) ;
  OPTION(blob_options, boltzmann,        false) ;
  OPTION(blob_options, delta_front,        0.3) ;
  OPTION(blob_options, conserve_momentum, true) ;
  OPTION(blob_options, A_relative_to_bg,  true) ;
  if(!isothermal){
    OPTION(blob_options, delta_perp_T,          10) ;
    OPTION(blob_options, elongation_T,           1) ;
    OPTION(blob_options, angle_blob_T,           0) ;
    OPTION(blob_options, A_T,                  2.0) ;
    OPTION(blob_options, L_b_T,                0.5) ;
    OPTION(blob_options, xoffset_T,           0.25) ;
    OPTION(blob_options, zoffset_T,            0.5) ;
    OPTION(blob_options, delta_front_T,        0.3) ;
  }
  
  Field3D x, y, z, x_can, z_can, x_can_T, z_can_T, n_blob, T_blob, nV_eq, nU_eq ;
 
  x = FieldFactory::get()->create3D("x", NULL, mesh, CELL_CENTRE, 0);
  y = FieldFactory::get()->create3D("y", NULL, mesh, CELL_CENTRE, 0)/TWOPI ;
  z = FieldFactory::get()->create3D("z", NULL, mesh, CELL_CENTRE, 0)/TWOPI ;

  if (normalise_lengths){      // delta_perp is in SI
    delta_perp = delta_perp/rho_s ;
    if(!isothermal){
        delta_perp_T /= rho_s;
    }
  }
  
  x_can =  (x-xoffset)*cos(angle_blob) +(z-zoffset)*sin(angle_blob) ;
  z_can = -(x-xoffset)*sin(angle_blob) +(z-zoffset)*cos(angle_blob) ;
  
  BoutReal delta_x = delta_perp/Lx ; 
  BoutReal delta_z = elongation*delta_perp/Lz ; 

  BoutReal delta_x_T = delta_perp_T/Lx ; 
  BoutReal delta_z_T = elongation_T*delta_perp_T/Lz ;
  
  if(!isothermal){
  x_can_T =  (x-xoffset_T)*cos(angle_blob_T) +(z-zoffset_T)*sin(angle_blob_T) ;
  z_can_T = -(x-xoffset_T)*sin(angle_blob_T) +(z-zoffset_T)*cos(angle_blob_T) ;
  }
  
  if (A_relative_to_bg){
    //Get reference temperature and density from midplane equilibrium value
    if (equilibrium_source != "1d_profiles")
      throw BoutException("Cannot use A_relative_to_bg unless equilibrium_source==\"1d_profiles\"");
    BoutReal *n_array,*T_array;
    n_array = new BoutReal[mesh->GlobalNy];
    T_array = new BoutReal[mesh->GlobalNy];
    BoutReal nref = -1.;
    BoutReal Tref = -1.;
    read_equilibrium_file(n_array,"n_eq.dat");
    if(symmetry_plane){
      nref = n_array[mesh->ystart];
    } else {
      nref = n_array[mesh->GlobalNy/2];
    }
    delete [] n_array;

    if (!isothermal) {
      read_equilibrium_file(T_array,"T_eq.dat");
      if(symmetry_plane){
        Tref = T_array[mesh->ystart];
      } else {
        Tref = T_array[mesh->GlobalNy/2];
      }
      delete [] T_array;
    }

    if (symmetry_plane) {
      n_blob = nref*A*exp(-SQ(x_can/delta_x))*exp(-SQ(z_can/delta_z))*0.5*(1.0 - tanh((y-L_b)/delta_front)) ;
      if(!isothermal){
        T_blob = Tref*A_T*exp(-SQ(x_can_T/delta_x_T))*exp(-SQ(z_can_T/delta_z_T))*0.5*(1.0 - tanh((y-L_b_T)/delta_front_T)) ; 
      }
    }else {
      n_blob = nref*A*exp(-SQ(x_can/delta_x))*exp(-SQ(z_can/delta_z))*0.5*(1.0 - tanh((y-0.5*(1.+L_b))/(delta_front/2)))*0.5*(1.0 - tanh((0.5*(1.-L_b)-y)/(delta_front/2))) ;
      if(!isothermal){
        T_blob = Tref*A_T*exp(-SQ(x_can_T/delta_x_T))*exp(-SQ(z_can_T/delta_z_T))*0.5*(1.0 - tanh((y-0.5*(1.+L_b_T))/(delta_front_T/2)))*0.5*(1.0 - tanh((0.5*(1.-L_b_T)-y)/(delta_front_T/2))) ;
      }
    }
  }else{
    if (symmetry_plane) {
      n_blob = A*exp(-SQ(x_can/delta_x))*exp(-SQ(z_can/delta_z))*0.5*(1.0 - tanh((y-L_b)/delta_front)) ;
      if(!isothermal){
         T_blob = A_T*exp(-SQ(x_can_T/delta_x_T))*exp(-SQ(z_can_T/delta_z_T))*0.5*(1.0 - tanh((y-L_b_T)/delta_front_T)) ;
      }
    }else {
      n_blob = A*exp(-SQ(x_can/delta_x))*exp(-SQ(z_can/delta_z))*0.5*(1.0 - tanh((y-0.5*(1.+L_b))/(delta_front/2)))*0.5*(1.0 - tanh((0.5*L_b-y)/(delta_front/2))) ;
      if(!isothermal){
         T_blob = A_T*exp(-SQ(x_can_T/delta_x_T))*exp(-SQ(z_can_T/delta_z_T))*0.5*(1.0 - tanh((y-0.5*(1.+L_b_T))/(delta_front_T/2)))*0.5*(1.0 - tanh((0.5*L_b_T-y)/(delta_front_T/2))) ;
      }
    }
  }

  if(conserve_momentum){
    mesh->communicate(n) ; 
    n_stag = interp_to(n, CELL_YLOW) ;
    nV_eq = n_stag*V ;
    nU_eq = n_stag*U ; 
    n = n + n_blob ;   
    mesh->communicate(n) ; 
    n_stag = interp_to(n, CELL_YLOW) ; 

    V = nV_eq/n_stag ; 
    U = nU_eq/n_stag ;

  }
  else{
    n = n + n_blob ; 
  }
  if(!isothermal){
      T += T_blob;
      mesh->communicate(T);
      T_stag = interp_to(T,CELL_YLOW);
    } 
  if(boltzmann){
    Field3D phi_blob ; 
    phi_blob = log(n) ; 
    vort = Delp2(phi_blob) ; 
  }
}

////////////////////////////////////////////////////////////////////////

void STORM::add_noise(){
  //Add noise to the initial profiles of n, vort, and T
  int MYPE = BoutComm::rank();
  srand (MYPE);
  for (int i = 0; i < mesh->LocalNx ; i++){
    for (int j = 0; j < mesh->LocalNy ; j++){
      for (int k = 0; k < mesh->LocalNz ; k++){
        n(i,j,k)    += 2.*(((double) rand()/(RAND_MAX)) - 0.5)*0.0001;
        vort(i,j,k) += 2.*(((double) rand()/(RAND_MAX)) - 0.5)*0.000001;
        if (!isothermal){
          T(i,j,k)  += 2.*(((double) rand()/(RAND_MAX)) - 0.5)*0.0001;
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////

// Add a function given in the input file to a field
void STORM::add_perturbation(Field3D &f, const std::string name) {
  Options *options = Options::getRoot()->getSection("filaments");
  Field3D perturbation = FieldFactory::get()->create3D(name, options, mesh, f.getLocation(), 0);
  f += perturbation;
}

// For consistency, U and V should have the same radial boundary conditions, as
// V's boundary conditions are applied to UmV_centre
void STORM::check_U_V_x_boundary_conditions() {
  auto options = *globalOptions;
  ASSERT0(options["U"]["bndry_xin"].as<std::string>()
      == options["V"]["bndry_xin"].as<std::string>());
  ASSERT0(options["U"]["bndry_xout"].as<std::string>()
      == options["V"]["bndry_xout"].as<std::string>());
}
