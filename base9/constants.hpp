#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "Base9Config.h"

const double EPSILON = 1.e-15;
const double EPS     = 0.000001;
const int MAG_LIST   = 26;   // list length for mag vs. S/N
const int MODEL_LIST = 118;  // age/Mv/V-I entries in WD model array

enum StarStatus { MSRG = 1    // Main-sequence or red giant
                , WD   = 3    // White dwarf
                , NSBH = 4    // Neutron star or black hole
                , BD   = 5    // Brown dwarf/planet
                , DNE  = 9    // Does not exist (zams mass > 100 Msun or < 0.0001)
};

const double R_sun   = 6.96342e10;

const int NPARAMS            = 9;

// Define model sets
//IFMR
const int WEIDEMANN  = 0;
const int WILLIAMS   = 1;
const int SALARISLIN = 2;
const int SALARISPW  = 3;
const int LINEAR     = 4;

//WD Cooling
enum class WdModel {WOOD = 0, MONTGOMERY = 1, ALTHAUS = 2, RENEDO = 3, NEW_MONTGOMERY = 4};
enum class WdAtmosphereModelSet {BERGERON, BERGERON_2019, BERGERON_2020};

//Brown Dwarf/Planet
const int NONE       = 0;
const int BARAFFE    = 1;

const int ZAMS       = 0;
const int NOW        = 1;

const int LOW        = 0;
const int HIGH       = 1;

enum class WdAtmosphere {DA, DB};

//Be careful adding new sample types.  Their order is important.
//There are a few places in the code that test for SAMPLE_TYPE > or <.
const int MASS              = -2;
const int AGE_DURING_WANDER = -1;  // age sampling during AGE_WANDER (defined in samplers.h)

const int AGE               =  0;  // age sampling
const int YYY               =  1;  // helium sampling
const int FEH               =  2;  // metallicity sampling
const int MOD               =  3;  // modulus sampling
const int ABS               =  4;  // absorption sampling;
const int CARBONICITY       =  5;
const int IFMR_INTERCEPT    =  6;
const int IFMR_SLOPE        =  7;
const int IFMR_QUADCOEF     =  8;

// multiPopMCMC makes use of this to keep track of multiple Y values without major changes
const int YYA               = YYY;
const int YYB               = IFMR_INTERCEPT;
const int LAMBDA            = IFMR_SLOPE;

const double DOF    = 6.0;
const double GAMMA6 = -2.0590305444197083417635;   /* GAMMA for DOF=6 */

#endif
