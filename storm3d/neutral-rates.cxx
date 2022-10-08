
#include "neutral-rates.hxx"

HydrogenRates *HydrogenRates::create(Options &options){
  // Decide which neutral rate model to use
  std::string ratesmodel  = options["ratemodel"].withDefault<std::string>("model1");
  if(ratesmodel == "model1") {
    // model1 : Rates supplied by Eva Havlicova
    return new HydrogenModel1();
  }
  else if(ratesmodel == "model2") {
    // model2 : Fitted by Hannah Willett
    return new HydrogenModel2();
  }
  else if(ratesmodel == "model3") {
    // model3 : Fitted by Fabio Riva
    return new HydrogenModel3();
  }
  else {
    throw BoutException("Unrecognised neutral rate model '%s'", ratesmodel.c_str());
  }
}

// Collision rate coefficient <sigma*v> [m3/s] from fitting Table II, page 40, column Q0 of Atomic and molecular processes in fusion edge plasmas
BoutReal HydrogenRates::electronNeutrals(const BoutReal Te) {
  BoutReal T_log = log10(Te);
  return 1.E-6*pow(10., 8.5189E-6*pow(T_log,10.) + 9.2379E-5*pow(T_log,9.) - 0.00066712*pow(T_log,8.)
                        - 0.00037931*pow(T_log,7.) + 0.005678*pow(T_log,6.) + 0.0015285*pow(T_log,5.)
                        - 0.009523*pow(T_log,4.) - 0.028693*pow(T_log,3.) - 0.18579*SQ(T_log)
                        + 0.12272*T_log - 6.8904);
}

/////////////////////////////////////////////////////////////////////////////

// Collision rate coefficient <sigma*v> [m3/s]
BoutReal HydrogenModel1::ionisation(const BoutReal T, const BoutReal UNUSED(n)) {
  double fION; // collision rate coefficient <sigma*v> [m3/s]
  double TT, X, S;

  TT = T;

  if (TT < 1.0)
    TT = 1.0;
  X = log10(TT);

  if (TT >= 20.0)
    S = -0.5151 * X - 2.563 / X - 5.231;
  else
    S = -3.054 * X - 15.72 * exp(-X) + 1.603 * exp(-X * X);

  fION = pow(10.0, S - 6.0);

  return fION;
}

//<sigma*v> [m3/s]
BoutReal HydrogenModel1::recombination(const BoutReal n, const BoutReal T) {
  double fREC; //<sigma*v> [m3/s]
  double TT, RDNE, RTE, DNE, E, RN, RT, RNJ, RTI, suma;
  int i, j, i1, j1;

  if (n < 1e3) // Log(n) used, so prevent NaNs
    return 0.0;

  double MATA[9][9] = {
      {
          -2.855728479302E+01, -7.664042607917E-01, -4.930424003280E-03,
          -5.386830982777E-03, -1.626039237665E-04, 6.080907650243E-06,
          2.101102051942E-05, -2.770717597683E-06, 1.038235939800E-07,
      },
      {
          3.488563234375E-02, -3.583233366133E-03, -3.620245352252E-03,
          -9.532840484460E-04, 1.888048628708E-04, -1.014890683861E-05,
          2.245676563601E-05, -4.695982369246E-06, 2.523166611507E-07,
      },
      {
          -2.799644392058E-02, -7.452514292790E-03, 6.958711963182E-03,
          4.631753807534E-04, 1.288577690147E-04, -1.145028889459E-04,
          -2.245624273814E-06, 3.250878872873E-06, -2.145390398476E-07,
      },
      {
          1.209545317879E-02, 2.709299760454E-03, -2.139257298118E-03,
          -5.371179699661E-04, -1.634580516353E-05, 5.942193980802E-05,
          -2.944873763540E-06, -9.387290785993E-07, 7.381435237585E-08,
      },
      {
          -2.436630799820E-03, -7.745129766167E-04, 4.603883706734E-04,
          1.543350502150E-04, -9.601036952725E-06, -1.211851723717E-05,
          1.002105099354E-06, 1.392391630459E-07, -1.299713684966E-08,
      },
      {
          2.837893719800E-04, 1.142444698207E-04, -5.991636837395E-05,
          -2.257565836876E-05, 3.425262385387E-06, 1.118965496365E-06,
          -1.291320799814E-07, -1.139093288575E-08, 1.265189576423E-09,
      },
      {
          -1.886511169084E-05, -9.382783518064E-06, 4.729262545726E-06,
          1.730782954588E-06, -4.077019941998E-07, -4.275321573501E-08,
          7.786155463269E-09, 5.178505597480E-10, -6.854203970018E-11,
      },
      {
          6.752155602894E-07, 3.902800099653E-07, -1.993485395689E-07,
          -6.618240780594E-08, 2.042041097083E-08, 3.708616111085E-10,
          -2.441127783437E-10, -9.452402157390E-12, 1.836615031798E-12,
      },
      {
          -1.005893858779E-08, -6.387411585521E-09, 3.352589865190E-09,
          1.013364275013E-09, -3.707977721109E-10, 7.068450112690E-12,
          3.773208484020E-12, -4.672724022059E-14, -1.640492364811E-14,
      },
  };

  RDNE = n;
  RTE = T;

  DNE = RDNE;
  TT = RTE;
  E = DNE * 1.0E-14;

  if (TT < 1.0)
    TT = 1.0;

  RN = log(E);
  RT = log(TT);

  suma = 0.0;
  for (i = 1; i <= 9; i++) {
    i1 = i - 1;
    for (j = 1; j <= 9; j++) {
      j1 = j - 1;
      RNJ = pow(RN, j1);
      if ((RN == 0.0) && (j1 == 0))
        RNJ = 1.0;
      RTI = pow(RT, i1);
      if ((RT == 0.0) && (i1 == 0))
        RTI = 1.0;
      suma = suma + MATA[j - 1][i - 1] * RNJ * RTI;
    }
  }

  fREC = exp(suma) * 1.0E-6;

  return fREC;
}

// <sigma*v> [m3/s]
BoutReal HydrogenModel1::chargeExchange(const BoutReal Te) {
  double fCX; //<sigma*v> [m3/s]
  double TT, S;

  TT = Te;

  if (TT < 1.0)
    TT = 1.0;
  S = -14.0 + log10(TT) / 3.0;

  fCX = pow(10.0, S);

  return fCX;
}

/////////////////////////////////////////////////////////////////////////////

// Collision rate coefficient <sigma*v> [m3/s]
BoutReal HydrogenModel2::ionisation(const BoutReal T, const BoutReal UNUSED(n)) {
  double fION; // Rate coefficient
  double TT;

  TT = T;

  double ioncoeffs[9] = {-3.271397E1,  1.353656E1,   -5.739329,
                         1.563155,     -2.877056E-1, 3.482560e-2,
                         -2.631976E-3, 1.119544E-4,  -2.039150E-6};

  double lograte = 0.0;
  for (int i = 0; i <= 8; i++) {
    lograte = lograte + ioncoeffs[i] * pow(log(TT), i);
  }

  fION = exp(lograte) * 1.0E-6;

  return fION;
}

//<sigma*v> [m3/s]
BoutReal HydrogenModel2::recombination(const BoutReal n, const BoutReal T) {
  double TT, RDNE, RTE, DNE, E, RN, RT, RNJ, RTI, suma, fHAV, fRAD;
  int i, j, i1, j1;

  if (n < 1e3) // Log(n) used, so prevent NaNs
    return 0.0;

  double MATA[9][9] = {
      {
          -2.855728479302E+01, -7.664042607917E-01, -4.930424003280E-03,
          -5.386830982777E-03, -1.626039237665E-04, 6.080907650243E-06,
          2.101102051942E-05, -2.770717597683E-06, 1.038235939800E-07,
      },
      {
          3.488563234375E-02, -3.583233366133E-03, -3.620245352252E-03,
          -9.532840484460E-04, 1.888048628708E-04, -1.014890683861E-05,
          2.245676563601E-05, -4.695982369246E-06, 2.523166611507E-07,
      },
      {
          -2.799644392058E-02, -7.452514292790E-03, 6.958711963182E-03,
          4.631753807534E-04, 1.288577690147E-04, -1.145028889459E-04,
          -2.245624273814E-06, 3.250878872873E-06, -2.145390398476E-07,
      },
      {
          1.209545317879E-02, 2.709299760454E-03, -2.139257298118E-03,
          -5.371179699661E-04, -1.634580516353E-05, 5.942193980802E-05,
          -2.944873763540E-06, -9.387290785993E-07, 7.381435237585E-08,
      },
      {
          -2.436630799820E-03, -7.745129766167E-04, 4.603883706734E-04,
          1.543350502150E-04, -9.601036952725E-06, -1.211851723717E-05,
          1.002105099354E-06, 1.392391630459E-07, -1.299713684966E-08,
      },
      {
          2.837893719800E-04, 1.142444698207E-04, -5.991636837395E-05,
          -2.257565836876E-05, 3.425262385387E-06, 1.118965496365E-06,
          -1.291320799814E-07, -1.139093288575E-08, 1.265189576423E-09,
      },
      {
          -1.886511169084E-05, -9.382783518064E-06, 4.729262545726E-06,
          1.730782954588E-06, -4.077019941998E-07, -4.275321573501E-08,
          7.786155463269E-09, 5.178505597480E-10, -6.854203970018E-11,
      },
      {
          6.752155602894E-07, 3.902800099653E-07, -1.993485395689E-07,
          -6.618240780594E-08, 2.042041097083E-08, 3.708616111085E-10,
          -2.441127783437E-10, -9.452402157390E-12, 1.836615031798E-12,
      },
      {
          -1.005893858779E-08, -6.387411585521E-09, 3.352589865190E-09,
          1.013364275013E-09, -3.707977721109E-10, 7.068450112690E-12,
          3.773208484020E-12, -4.672724022059E-14, -1.640492364811E-14,
      },
  };

  RDNE = n;
  RTE = T;

  DNE = RDNE;
  TT = RTE;
  E = DNE * 1.0E-14;

  RN = log(E);
  RT = log(TT);

  suma = 0.0;
  for (i = 1; i <= 9; i++) {
    i1 = i - 1;
    for (j = 1; j <= 9; j++) {
      j1 = j - 1;
      RNJ = pow(RN, j1);
      if ((RN == 0.0) && (j1 == 0))
        RNJ = 1.0;
      RTI = pow(RT, i1);
      if ((RT == 0.0) && (i1 == 0))
        RTI = 1.0;
      suma = suma + MATA[j - 1][i - 1] * RNJ * RTI;
    }
  }

  fHAV = exp(suma) * 1.0E-6 / (1.0 + 0.125 * TT);

  double A = 3.92E-20;
  double B = 3.0E-124 * pow(1.6E-19, -4.5);
  double Ry = 13.60569;
  double chi = 0.35;

  fRAD = A * pow(Ry, 1.5) * 1.0 / (sqrt(TT) * (Ry + chi * TT));

  double fREC = fHAV + fRAD + B * DNE * pow(TT, -5.0);

  return fREC;
}

// <sigma*v> [m3/s]
BoutReal HydrogenModel2::chargeExchange(const BoutReal T) {
  double fCX;
  double TT;

  TT = T;

  double E = 10.0;

  double cxcoeffs[9][9] = {
      {
          -1.829079582E1, 1.640252721E-1, 3.364564509E-2, 9.530225559E-3,
          -8.519413900E-4, -1.247583861E-3, 3.014307546E-4, -2.499323170E-5,
          6.932627238E-7,
      },
      {
          2.169137616E-1, -1.106722014E-1, -1.382158680E-3, 7.348786287E-3,
          -6.343059502E-4, -1.919569450E-4, 4.075019352E-5, -2.850044983E-6,
          6.966822400E-8,
      },
      {
          4.307131244E-2, 8.948693625E-3, -1.209480567E-2, -3.675019470E-4,
          1.039643391E-3, -1.553840718E-4, 2.670827249E-6, 7.695300598E-7,
          -3.783302282E-8,
      },
      {
          -5.754895093E-4, 6.062141761E-3, 1.075907882E-3, -8.119301728E-4,
          8.911036876E-6, 3.175388950E-5, -4.515123642E-6, 2.187439284E-7,
          -2.911233952E-9,
      },
      {
          -1.552077120E-3, -1.210431588E-3, 8.297212634E-4, 1.361661817E-4,
          -1.008928628E-4, 1.080693990E-5, 5.106059414E-7, -1.299275586E-7,
          5.117133050E-9,
      },
      {
          -1.876800283E-4, -4.052878752E-5, -1.907025663E-4, 1.141663042E-5,
          1.775681984E-5, -3.149286924E-6, 3.105491555E-8, 2.274394089E-8,
          -1.130988251E-9,
      },
      {
          1.125490271E-4, 2.875900436E-5, 1.338839629E-5, -4.340802793E-6,
          -7.003521917E-7, 2.318308730E-7, -6.030983538E-9, -1.755944926E-9,
          1.005189187E-10,
      },
      {
          -1.238982763E-5, -2.616998140E-6, -1.171762874E-7, 3.517971869E-7,
          -4.928692833E-8, 1.756388999E-10, -1.446756796E-10, 7.143183138E-11,
          -3.989884106E-12,
      },
      {
          4.163596197E-7, 7.558092849E-8, -1.328404104E-8, -9.170850254E-9,
          3.208853884E-9, -3.952740759E-10, 2.739558476E-11, -1.693040209E-12,
          6.388219930E-14,
      },
  };

  double lograte = 0.0;
  for (int i = 0; i <= 8; i++) {
    for (int j = 0; j <= 8; j++) {
      lograte = lograte + cxcoeffs[i][j] * pow(log(TT), i) * pow(log(E), j);
    }
  }

  fCX = 1.0E-6 * exp(lograte);

  return fCX;
}

/////////////////////////////////////////////////////////////////////////////

// Collision rate coefficient <sigma*v> [m3/s]
// Obtained by fitting the scd12_h.dat file from the open-adas database with n in [1E17m^-3,1E20m^-3] and Te in [1eV,200eV]
// The fit is 2nd order in logn and 5th order in logT -> fit error smaller than 3%
BoutReal HydrogenModel3::ionisation(const BoutReal Te, const BoutReal n) {
  if(Te < 1.){
    // <sigma*v> < 1e-19 for Te < 1
    return 0.;
  } else {
    BoutReal logT = log(Te);
    BoutReal logn = log(n/5.E17);
    if(n > 1.E20) logn = log(200);

    return 1.E-14*exp(-13.22 + 0.1617*logn + 13.19*logT + 0.01919*SQ(logn) - 0.09323*logn*logT - 5.762*SQ(logT)
                      -0.008934*SQ(logn)*logT + 0.0244*logn*SQ(logT) + 1.463*pow(logT,3.) + 0.002132*SQ(logn*logT)
                      -0.003364*logn*pow(logT,3.) - 0.2004*pow(logT,4.) - 0.0001823*SQ(logn*logT)*logT +
                      0.0002076*logn*pow(logT,4.) + 0.01117*pow(logT,5.));
  }
}

// Collision rate coefficient <sigma*v> [m3/s]
// Obtained by fitting the acd12_h.dat file from the open-adas database with n in [1E17m^-3,1E20m^-3] and Te in [0.2eV,20eV]
// The fit is 5th order in logn and 5th order in logT -> fit error smaller than 5%
BoutReal HydrogenModel3::recombination(const BoutReal n, const BoutReal Te) {
  if(Te > 20.){
    // <sigma*v> < 3e-20 for Te > 20
    return 0.;
  } else {
    BoutReal logT = log(Te);
    if(Te < 0.2) logT = log(0.2);
    BoutReal logn = log(n/5.E17);
    if(n < 1.E17) logn = log(0.2);
    if(n > 1.E20) logn = log(200);

    return 1.E-14*exp(-9.807 + 0.08922*logn - 1.103*logT + 0.01242*SQ(logn) - 0.09195*logn*logT +
                      0.1579*SQ(logT) + 0.002846*pow(logn,3.) - 0.01318*SQ(logn)*logT + 0.03709*logn*SQ(logT)
                      -0.04817*pow(logT,3.) - 0.0001848*pow(logn,4.) - 0.001441*pow(logn,3.)*logT + 0.0039*SQ(logn*logT)
                      -0.004297*logn*pow(logT,3.) + 0.002168*pow(logn,4.) - 0.00001058*pow(logn,5.) +
                      0.0001494*pow(logn,4.)*logT - 0.00001184*logn*SQ(logn*logT) - 0.0001989*logT*SQ(logn*logT)
                      -0.0003826*logn*pow(logT,4.) + 0.0007565*pow(logT,5.));

  }
}

// Collision rate coefficient <sigma*v> [m3/s]
// Obtained by fitting the ccd96_h.dat file from the open-adas database with Te in [0.2eV,200eV]
// The fit is 3rd order in logT -> fit error smaller than 2%
BoutReal HydrogenModel3::chargeExchange(const BoutReal Te) {
  BoutReal logT = log(Te);
  if(Te < 0.2) logT = log(0.2);
  return 1.E-14*exp(0.002098*pow(logT,3.) - 0.02461*SQ(logT) + 0.4602*logT - 0.3722);
}

