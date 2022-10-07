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
#include <bout_types.hxx>
#include <derivs.hxx>
#include <interpolation.hxx>

const Field3D STORM::Grad_par_EM(const Field3D &var_aligned, const Field3D &var_outloc,
    const Field3D &Psi, CELL_LOC outloc, const std::string& method) {
  AUTO_TRACE();

  ASSERT1(var_aligned.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(var_outloc.getDirectionY() == YDirectionType::Standard);
  ASSERT1(Psi.getDirectionY() == YDirectionType::Standard);

  Field3D result = fromFieldAligned(Grad_par(var_aligned, outloc, method), "RGN_NOBNDRY");

  if (electromagnetic) {
    ASSERT1(Psi.getLocation() == outloc and var_outloc.getLocation() == outloc);

    result -= beta0/2 *bracket(Psi,var_outloc, bm, outloc);
  }

  return result;
}

const Field3D STORM::Div_par_EM(const Field3D &var_aligned, const Field3D &var_outloc,
    const Field3D &Psi, CELL_LOC outloc, const std::string& method) {
  AUTO_TRACE();

  ASSERT1(var_aligned.getDirectionY() == YDirectionType::Aligned);

  if (electromagnetic) {
    ASSERT1(var_outloc.getDirectionY() == YDirectionType::Standard);
    ASSERT1(Psi.getDirectionY() == YDirectionType::Standard);

    auto coords_outloc = mesh->getCoordinates(outloc);
    return coords_outloc->Bxy*Grad_par_EM(var_aligned/var_aligned.getCoordinates()->Bxy,
        var_outloc/coords_outloc->Bxy, Psi, outloc, method);
  } else {
    return fromFieldAligned(Div_par(var_aligned, outloc, method), "RGN_NOBNDRY");
  }
}

const Field3D STORM::Vpar_Grad_par_EM(const Field3D &v_aligned, const Field3D &f_aligned,
    const Field3D &v_outloc, const Field3D &f, const Field3D &Psi, CELL_LOC outloc,
    const std::string& method) {
  AUTO_TRACE();

  ASSERT1(v_aligned.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(f_aligned.getDirectionY() == YDirectionType::Aligned);
  ASSERT1(v_outloc.getDirectionY() == YDirectionType::Standard);
  ASSERT1(f.getDirectionY() == YDirectionType::Standard);
  ASSERT1(Psi.getDirectionY() == YDirectionType::Standard);

  Field3D result = fromFieldAligned(Vpar_Grad_par(v_aligned, f_aligned, outloc, method),
                                    "RGN_NOBNDRY");

  if (electromagnetic) {
    ASSERT1(Psi.getLocation() == outloc and f.getLocation() == outloc);

    result -= beta0/2. * v_outloc * bracket(Psi, f, bm, outloc);
  }

  return result;
}

const Field3D STORM::Grad_perp_dot_Grad_perp(const Field3D &f, const Field3D &g) {

  AUTO_TRACE();

  auto& coords = *f.getCoordinates();

  auto ddx_f = DDX(f);
  auto ddz_f = DDZ(f);
  auto ddx_g = DDX(g);
  auto ddz_g = DDZ(g);
  return  coords.g11*ddx_f*ddx_g + coords.g33*ddz_f*ddz_g
          + coords.g13*(ddx_f*ddz_g + ddz_f*ddx_g);
}

const Field3D STORM::Grad_perp_dot_Grad_perp(const Field3D &f) {

  AUTO_TRACE();

  auto& coords = *f.getCoordinates();

  auto ddx_f = DDX(f);
  auto ddz_f = DDZ(f);
  return  coords.g11*SQ(ddx_f) + coords.g33*SQ(ddz_f)
          + 2.*coords.g13*ddx_f*ddz_f;
}

// Curvature operator
const Field3D STORM::Curv(const Field3D& f, const Field3D& f_aligned){
  AUTO_TRACE();

  ASSERT1(f.getDirectionY() == YDirectionType::Standard);
  ASSERT1(f_aligned.getDirectionY() == YDirectionType::Aligned);

  if (realistic_geometry == "none") {
    return g0*DDZ(f)/sqrt(f.getCoordinates()->g_33);
  } else {
    ASSERT0(f.getLocation() == CELL_CENTRE);
    Vector3D Grad_f;
    Grad_f.covariant = true;
    Grad_f.x = DDX(f);
    Grad_f.y = fromFieldAligned(DDY(f_aligned), "RGN_NOBNDRY");
    Grad_f.z = DDZ(f);
    return Curlb_B*Grad_f;
  }
}

// Interpolate in y, but if double null geometry change the scheme near X-points
const Field3D STORM::interp_to_fixstag(const Field3D& var_fa, CELL_LOC loc){
  if (loc != CELL_YLOW || var_fa.getLocation() != CELL_CENTRE) {
    throw BoutException("interp_to_fixstag is implemented only from CELL_CENTRE to CELL_YLOW!");
  }
  ASSERT0(var_fa.getMesh()->ystart >= 2);

  if(realistic_geometry != "singlenull" && realistic_geometry != "doublenull")
      throw BoutException("interp_to_fixstag implemented only for realistic_geometry==singlenull or realistic_geometry==doublenull"); 

  Field3D result = interp_to(var_fa, loc, "RGN_NOBNDRY");

  // Change the stencil from centered to left/right near the X-points
  int iy;
  if( mesh->getGlobalYIndex(mesh->ystart) == jyseps1_1+1 || 
      mesh->getGlobalYIndex(mesh->ystart) == jyseps1_2+1 ||
      mesh->getGlobalYIndex(mesh->ystart) == jyseps2_1+1 ||
      mesh->getGlobalYIndex(mesh->ystart) == jyseps2_2+1){
    iy = mesh->ystart;
    for(int ix=mesh->xstart; ix<=mesh->xend; ++ix){
      for(int iz=0; iz<mesh->LocalNz; ++iz){
        result(ix,iy,iz) = (15.*var_fa(ix,iy,iz) - 10.*var_fa(ix,iy+1,iz) + 3.*var_fa(ix,iy+2,iz))/8.;
        result(ix,iy+1,iz) = (3.*var_fa(ix,iy,iz) + 6.*var_fa(ix,iy+1,iz) - var_fa(ix,iy+2,iz))/8.;
      }
    }
  }
  if( mesh->getGlobalYIndex(mesh->yend) == jyseps2_1 ||
      mesh->getGlobalYIndex(mesh->yend) == jyseps2_2 ||
      mesh->getGlobalYIndex(mesh->yend) == jyseps1_1 ||
      mesh->getGlobalYIndex(mesh->yend) == jyseps1_2){
    iy = mesh->yend;
    for(int ix=mesh->xstart; ix<=mesh->xend; ++ix){
      for(int iz=0; iz<mesh->LocalNz; ++iz){
        result(ix,iy,iz) = (15.*var_fa(ix,iy,iz) - 10.*var_fa(ix,iy-1,iz) + 3.*var_fa(ix,iy-2,iz))/8.;
        result(ix,iy-1,iz) = (3.*var_fa(ix,iy,iz) + 6.*var_fa(ix,iy-1,iz) - var_fa(ix,iy-2,iz))/8.;
      }
    }
  }

  result = fromFieldAligned(result, "RGN_NOBNDRY");
  return result; 
}

