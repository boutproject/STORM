/*
 * Mixed diffusive-fluid model, similar to UEDGE model
 * 3D model, diffusive in X-Z and fluid in Y
 */

#ifndef __NEUTRAL_DVPAR_H__
#define __NEUTRAL_DVPAR_H__

#include "neutral-model.hxx"
#include "../shared/BoutEquation/equation.hxx"
#include <invert_laplace.hxx>
#include <interpolation.hxx>

class NeutralDVpar : public NeutralModel {
public:
  NeutralDVpar(Solver *solver, Options &options, Datafile &dump);
  ~NeutralDVpar() {}

  /// Update plasma quantities
  void update(const Field3D& n, const Field3D& n_stag, const Field3D& U, 
              const Field3D& V, const Field3D& T, const Field3D& T_stag,
              BoutReal time);
  
  void setRecycledFlux(const FieldPerp& ionflux_lower, const FieldPerp& ionflux_upper) {
    if (mesh->hasBndryLowerY()) {
        recycled_lower = Rc*ionflux_lower;
    }
    if (mesh->hasBndryUpperY()) {
        recycled_upper = Rc*ionflux_upper;
    }
  }
 
  void precon(BoutReal t, BoutReal gamma, BoutReal delta); 
private:
  Field3D nn, nvn;
  Field3D Tn;
  Field3D lognn;
  Field3D nn_stag, lognn_stag, Tn_stag;
  Field3D vn;
  Field3D lognn_aligned, nn_aligned, nvn_aligned, vn_aligned;
  Field3D Dn;
  Field3D vn_centre;
  FieldPerp recycled_lower, recycled_upper;
  FieldPerp nn_lower, nn_upper;
  FieldPerp vn_lower, vn_upper;
  Field2D profileSn_lower, profileSn_upper, profileSn_lower_stag, profileSn_upper_stag;
  Field2D Snn, Snvn;
  Field2D DDY_Jogyy_oJ;
  Field2D g_12og_22;
  Field2D g_23og_22;
  std::string orderupwindscheme;

  int rhs_counter = 0;
  Equation neutral_density_equation{lognn, "lognn", Options::root()["save_equations"],
                                    dump, rhs_counter};
  Equation neutral_momentum_equation{nvn, "nvn", Options::root()["save_equations"], dump,
                                     rhs_counter};

  BoutReal Tn0;            // Uniform temperature of neutrals
  BoutReal munvn;          // numerical dissipation in neutral momentum equation
  BoutReal lambdamax;      // maximum mean free path for neutrals
  int bndry_neutrals;      // = 1: bndries are applied as dirichlet, = 2: bndries are applied as sources
  BoutReal Rc;             // Recycled fration
  BoutReal dt_update;      // Do not update the collision frequencies every internal timestep, only every dt_update
  bool monitor_minmaxmean; // Flag for output monitor
  BoutReal lambdaSnn;      // Fall off length in y (poloidal plane) for neutral sources
  bool yderivatives;       // Flag to include y derivatives in perpendicular operators

  Mesh* mesh;
  Coordinates* coordinates_centre;
  Coordinates* coordinates_stag;

  bool precon_firststep = true; // flag to detect first step in precon
  Laplacian *inv; // Laplacian inversion used for preconditioning

  bool updaterates;
  BoutReal monitor_timelast = -1.;
  BoutReal minmaxmean_timelast = -1.;
  int timestepmonitor(BoutReal simtime);

  // Need global list of all NeutralDVpar instances, to iterate over in the
  // timestepmonitor_func. Needed because timestepmonitor_func must be passed as a
  // function pointer which can contain no state, so has to get the NeutralDVpar
  // instances from a global variable.
  static std::vector<NeutralDVpar*> neutral_dvpar_instances;
  static int timestepmonitor_func(Solver*, BoutReal simtime, BoutReal);

  void recycleFluxes(BoutReal time);
  void initialiseSource(Options &options);

  const Field3D myInterp_to(const Field3D& f, CELL_LOC outloc) {
    return fromFieldAligned(interp_to(f, outloc, "RGN_NOBNDRY"), "RGN_NOBNDRY");
  }
  const Field3D myGrad_par(const Field3D &var_aligned, CELL_LOC outloc=CELL_DEFAULT) {
    AUTO_TRACE();
    ASSERT1(var_aligned.getDirectionY() == YDirectionType::Aligned);
    return fromFieldAligned( Grad_par(var_aligned, outloc) , "RGN_NOBNDRY");
  }
  const Field3D myDiv_par(const Field3D &var_aligned, CELL_LOC outloc=CELL_DEFAULT) {
    AUTO_TRACE();
    ASSERT1(var_aligned.getDirectionY() == YDirectionType::Aligned);
    return fromFieldAligned( Div_par(var_aligned, outloc) , "RGN_NOBNDRY");
  }

  class OutputMonitor : public Monitor {
  public:
    OutputMonitor(NeutralDVpar* ndvp_in, BoutReal timestep=-1)
      : Monitor(timestep), neutral_dvpar(ndvp_in) {};
    int call(Solver *solver, BoutReal simtime, int iter, int NOUT) override;
  private:
    NeutralDVpar* neutral_dvpar;
  };
  OutputMonitor my_output_monitor;

  // Print minimum, maximum and mean of evolving variables. Used when
  // monitor_minmaxmean=true
  void printMinMaxMean();
};

#endif // __NEUTRAL_DVPAR_H__
