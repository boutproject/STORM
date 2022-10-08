/*!
 *
 */

#include "neutral-model.hxx"
#include "D-vpar.hxx"

NeutralModel *NeutralModel::create(Solver *solver, Options &options, Datafile &dump) {
  // Decide which neutral model to use
  std::string type  = options["type"].withDefault<std::string>("none");

  if (type == "none") {
    // Neutral model which does nothing
    return NULL;
  } else if (type == "d-vpar") {
    // Diffusive in X-Z, fluid in Y
    return new NeutralDVpar(solver, options, dump);
  }
  throw BoutException("Unrecognised neutral model '%s'", type.c_str());
}

/*!
 * Atomic processes
 */
void NeutralModel::neutral_rates(
    const Field3D& n, const Field3D& n_stag, const Field3D& U, const Field3D& V, const Field3D &T, const Field3D& T_stag, // Plasma quantities
    const Field3D& nn, const Field3D& nn_stag, const Field3D& vn, // Neutral gas
    bool updaterates) {

  if (updaterates) {
    // Charge exchange, recombination, ionization and electron-neutral collision frequencies
    BOUT_FOR(i, n.getRegion("RGN_NOBNDRY")) {
      // Plasma quantities in SI units
      BoutReal nplasma = n[i]*n0;
      BoutReal Tplasma = T[i]*T0;
      BoutReal nplasma_stag = n_stag[i]*n0;
      BoutReal Tplasma_stag = T_stag[i]*T0;

      // Charge exchange
      nucx[i] = nplasma*hydrogen->chargeExchange(Tplasma)/F0;
      nucx_stag[i] = nplasma_stag*hydrogen->chargeExchange(Tplasma_stag)/F0;

      // Recombination
      nurc[i] = nplasma*hydrogen->recombination(nplasma, Tplasma)/F0;
      nurc_stag[i] = nplasma_stag*hydrogen->recombination(nplasma_stag, Tplasma_stag)/F0;

      // Ionisation
      nuiz[i] = nplasma*hydrogen->ionisation(Tplasma, nplasma)/F0;
      nuiz_stag[i] = nplasma_stag*hydrogen->ionisation(Tplasma_stag, nplasma_stag)/F0;

      // Electron-neutral elastic collisions
      nuen_stag[i] = nplasma_stag*hydrogen->electronNeutrals(Tplasma_stag)/F0;
    }
  }

  // Reaction rates
  Rcx = nn*nucx;
  Rcx_stag = nn_stag*nucx_stag;
  Rrc = n*nurc;
  Rrc_stag = n_stag*nurc_stag;
  Riz = nn*nuiz;
  Riz_stag = nn_stag*nuiz_stag;

  S = Rrc - Riz;
  S_stag = Rrc_stag - Riz_stag;
  Fi = Rcx_stag*(U - vn) + Rrc_stag*U - Riz_stag*vn;
  Fe = Rrc_stag*V - Riz_stag*(2.*vn - V) + nn_stag*nuen_stag*(V-vn);
  Fperp = Rcx + Rrc;
  Rp = (1.09*T - 13.6/T0)*Rrc + (Eionize/T0)*Riz;
  //Rp = 1.5*T*Rrc + (Eionize/T0)*Riz;
}
