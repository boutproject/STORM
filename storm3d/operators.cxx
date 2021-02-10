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

// Curvature operator
const Field3D STORM::Curv(const Field3D &f){
  AUTO_TRACE();

  ASSERT1(f.getDirectionY() == YDirectionType::Standard);

  return g0*DDZ(f)/sqrt(f.getCoordinates()->g_33);
}
