
#ifndef __NEUTRALRATES_H__
#define __NEUTRALRATES_H__

#include <field3d.hxx>
#include <bout_types.hxx>
#include <options.hxx>

class HydrogenRates {
public:
  HydrogenRates() {}
  virtual ~HydrogenRates() {}

  static HydrogenRates* create(Options &options);

  // inputs in SI units, outputs in [m3/s]
  virtual BoutReal ionisation(const BoutReal Te, const BoutReal UNUSED(n)) = 0;
  virtual BoutReal recombination(const BoutReal n, const BoutReal Te) = 0;
  virtual BoutReal chargeExchange(const BoutReal Te) = 0;
  BoutReal electronNeutrals(const BoutReal Te);
private:

};

/// Rates supplied by Eva Havlicova
class HydrogenModel1: public HydrogenRates {
public:
  HydrogenModel1() : HydrogenRates() {}
  ~HydrogenModel1() {}

  // Collision rate coefficient <sigma*v> [m3/s]
  BoutReal ionisation(const BoutReal Te, const BoutReal UNUSED(n));
  
  //<sigma*v> [m3/s]
  BoutReal recombination(const BoutReal n, const BoutReal Te);
  
  // <sigma*v> [m3/s]
  BoutReal chargeExchange(const BoutReal Te);
  
private:
  
};

/*!
 * Hydrogen rates, fitted by Hannah Willett May 2015
 * University of York
 */
class HydrogenModel2: public HydrogenRates {
public:
  HydrogenModel2() : HydrogenRates() {}
  ~HydrogenModel2() {}

  // Ionisation rate coefficient <sigma*v> [m3/s]
  BoutReal ionisation(const BoutReal T, const BoutReal UNUSED(n)); 
  
  // Recombination rate coefficient <sigma*v> [m3/s]
  BoutReal recombination(const BoutReal n, const BoutReal T);
  
  // Charge exchange rate coefficient <sigma*v> [m3/s]
  BoutReal chargeExchange(const BoutReal Te);
  
private:
  
};

// Rates from fitting open-adas database
class HydrogenModel3: public HydrogenRates {
public:
  HydrogenModel3() : HydrogenRates() {}
  ~HydrogenModel3() {}

  // Ionisation rate coefficient <sigma*v> [m3/s]
  BoutReal ionisation(const BoutReal T, const BoutReal n);

  // Recombination rate coefficient <sigma*v> [m3/s]
  BoutReal recombination(const BoutReal n, const BoutReal T);

  // Charge exchange rate coefficient <sigma*v> [m3/s]
  BoutReal chargeExchange(const BoutReal Te);
private:

};


#endif // __NEUTRALRATES_H__
