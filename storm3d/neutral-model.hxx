/*
 * Base class for neutral gas model
 */

#ifndef __NEUTRAL_MODEL_H__
#define __NEUTRAL_MODEL_H__

#include "neutral-rates.hxx"
#include <bout/solver.hxx>
#include <bout/mesh.hxx>

class NeutralModel {
public:
  NeutralModel(Options &options) {

    // Better to initialize with NaN
    S  = 0.;
    S_stag = 0.;
    Fe = 0.;
    Fi = 0.;
    Fperp = 0.;
    Rp = 0;
    
    nuiz = 0.0;
    nurc = 0.0;
    nucx = 0.0;
    nuiz_stag = 0.;
    nurc_stag = 0.;
    nuen_stag = 0.;
    nucx_stag = 0.;

    Riz = 0.;
    Rrc = 0.;
    Rcx = 0.;
    Riz_stag = 0.;
    Rrc_stag = 0.;
    Rcx_stag = 0.;

    S_stag.setLocation(CELL_YLOW);
    nuiz_stag.setLocation(CELL_YLOW);
    nurc_stag.setLocation(CELL_YLOW);
    nucx_stag.setLocation(CELL_YLOW);
    nuen_stag.setLocation(CELL_YLOW);
    Riz_stag.setLocation(CELL_YLOW);
    Rrc_stag.setLocation(CELL_YLOW);
    Rcx_stag.setLocation(CELL_YLOW);
    Fe.setLocation(CELL_YLOW);
    Fi.setLocation(CELL_YLOW);

    SAVE_REPEAT3(nuiz,nurc,nucx);
    SAVE_REPEAT4(nuiz_stag,nurc_stag,nuen_stag,nucx_stag);

    OPTION(options, Eionize, 30);      // Energy loss per ionisation [eV]
    hydrogen = HydrogenRates::create(options);
  }
  virtual ~NeutralModel() {}
  
  /*!
   * Creates an instance of NeutralModel, based on given options
   */
  static NeutralModel* create(Solver *solver, Options &options);
  
  /*!
   * Set normalisations for temperature [eV], density [m^-3], 
   * length [m] and frequency [s^-1]
   * 
   */
  void InitialiseNeutrals(BoutReal Te, BoutReal Ne, BoutReal B, BoutReal length, BoutReal freq, BoutReal miome) { 
    T0 = Te; n0 = Ne; B0 = B; L0 = length; F0 = freq; mu = miome;
    output << "==== Using neutral model ====" << endl;
    output << "with T0 = " << T0 << " eV, n0 = " << n0 << "m^-3, B0 = " << B0 << "T, L0 = " << L0 << "m, F0 = " << F0 << "1/s" << endl;
  }

  /*!
   * Update plasma quantities
   */
  virtual void update(const Field3D& n, const Field3D& n_stag, const Field3D& U, 
                      const Field3D& V, const Field3D& T, const Field3D& T_stag,
                      BoutReal time) = 0;
  
  /*!
   *
   */
  virtual void setRecycledFlux(const FieldPerp& UNUSED(ionflux_lower), const FieldPerp& UNUSED(ionflux_upper)) {}
  
  /*!
   * Preconditioning
   */
  virtual void precon(BoutReal UNUSED(t), BoutReal UNUSED(gamma), BoutReal UNUSED(delta)) {}

  Field3D S;     // Plasma particle sink, neutral source
  Field3D S_stag;
  Field3D Fe;    // Electron-neutral friction
  Field3D Fi;    // Ion-neutral friction
  Field3D Fperp; // Ion-neutral friction in vorticity
  Field3D Rp;    // Radiation from the plasma

  Field3D Riz, Rrc, Rcx;
  Field3D Riz_stag, Rrc_stag, Rcx_stag;

protected:
  BoutReal T0, n0, B0, L0, F0; // Normalisations for temperature, density, magnetic field, lengths and frequencies
  BoutReal mu; // ion to electron mass ratio
  
  BoutReal Eionize;   // Energy loss per ionisation [eV]

  HydrogenRates *hydrogen; // Atomic rates

  Field3D nuiz;    // Ionization frequency
  Field3D nurc;    // Recombination frequency
  Field3D nucx;    // Charge-exchange collision frequency
  Field3D nuiz_stag;
  Field3D nurc_stag;
  Field3D nucx_stag;
  Field3D nuen_stag; // Electron-neutral elastic collision frequency
  
  void neutral_rates(const Field3D& n, const Field3D& n_stag, const Field3D& U, const Field3D& V, const Field3D& T, const Field3D& T_stag,
                     const Field3D& nn, const Field3D& nn_stag, const Field3D& vn,
                     bool updaterates = true);

private:
  NeutralModel();
};

#endif // __NEUTRAL_MODEL_H__
