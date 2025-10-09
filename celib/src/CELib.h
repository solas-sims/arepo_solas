#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

/*! \file CELib.h
 * \brief This file has all prototype declarations.
 */

void CELibInit(void);

/*!
 * @enum CELib IMF types
 * IMF types.
 */
enum {
    CELibIMF_Salpeter,
    CELibIMF_DietSalpeter,
    CELibIMF_MillerScalo,
    CELibIMF_Kroupa,
    CELibIMF_Kroupa1993,
    CELibIMF_Kennicutt,
    CELibIMF_Chabrier,
    CELibIMF_Susa,
    CELibIMF_NTypes,
};


/*!
 * @enum CELib solar abundance pattern 
 * Solar abundance pattern.
 */
enum {
    CELibSolarAbundancePattern_A09,  // Asplund et al. (2009)  
    CELibSolarAbundancePattern_GS98, // Grevesse & Sauval (1998)
    CELibSolarAbundancePattern_AG89, // Anders & Grevesse (1989) 
    CELibSolarAbundancePattern_Number,
};


/*!
 * @enum CELib Lifetime types
 * Lifetime types.
 */
enum {
    CELibLifeTime_Original, // Original data with linear interpolation.
    CELibLifeTime_LSF,      // Least square fitting
    CELibLifeTime_Number,
};


/*!
 * @enum CELib SN II yields table IDs
 * SN II yields table IDs.
 */
enum {
    CELibSNIIYieldsTableID_P98, // Portinari et al. 1998
    CELibSNIIYieldsTableID_N13, // Nomoto et al. 2013, default
    CELibSNIIYieldsTableID_Number,
};


/*!
 * @enum CELib SN Ia yields table IDs
 * SN Ia yields IDs.
 */
enum {
    CELibSNIaYieldsTableID_I99,   // Iwamoto et al. 1998
    CELibSNIaYieldsTableID_M10,   // Maeda et al. 2010
    CELibSNIaYieldsTableID_S13,   // Seitenzahl et al. 2013, default
    CELibSNIaYieldsTableID_T04,   // Travaglio et al. 2004
    CELibSNIaYieldsTableID_Number,
};

/*!
 * @enum CELib SN Ia yields table model IDs for Iwamoto et al. 1998
 * SN Ia yields table model IDs.
 */
enum {
    CELibSNIaYieldsTableModelID_I99_W7,   
    CELibSNIaYieldsTableModelID_I99_W70,  
    CELibSNIaYieldsTableModelID_I99_WDD1, 
    CELibSNIaYieldsTableModelID_I99_WDD2, 
    CELibSNIaYieldsTableModelID_I99_WDD3, 
    CELibSNIaYieldsTableModelID_I99_CDD1, 
    CELibSNIaYieldsTableModelID_I99_CDD2, 
    CELibSNIaYieldsTableModelID_I99_Number,
};

/*!
 * @enum CELib SN Ia yields table model IDs for Maeda et al. 2010
 * SN Ia yields table model IDs.
 */
enum {
    CELibSNIaYieldsTableModelID_M10_W7,   
    CELibSNIaYieldsTableModelID_M10_CDEF, 
    CELibSNIaYieldsTableModelID_M10_CDDT, 
    CELibSNIaYieldsTableModelID_M10_ODDT, 
    CELibSNIaYieldsTableModelID_M10_Number,
};

/*!
 * @enum CELib SN Ia yields table model IDs for Seitenzahl et al. 2013
 * SN Ia yields table model IDs.
 */
enum {
    CELibSNIaYieldsTableModelID_S13_N1,       
    CELibSNIaYieldsTableModelID_S13_N3,       
    CELibSNIaYieldsTableModelID_S13_N5,       
    CELibSNIaYieldsTableModelID_S13_N10,      
    CELibSNIaYieldsTableModelID_S13_N20,      
    CELibSNIaYieldsTableModelID_S13_N40,      
    CELibSNIaYieldsTableModelID_S13_N100H,    
    CELibSNIaYieldsTableModelID_S13_N100,     // default
    CELibSNIaYieldsTableModelID_S13_N100L,    
    CELibSNIaYieldsTableModelID_S13_N150,     
    CELibSNIaYieldsTableModelID_S13_N200,     
    CELibSNIaYieldsTableModelID_S13_N300C,    
    CELibSNIaYieldsTableModelID_S13_N1600,    
    CELibSNIaYieldsTableModelID_S13_N1600C,   
    CELibSNIaYieldsTableModelID_S13_N100Z05,  
    CELibSNIaYieldsTableModelID_S13_N100Z01,  
    CELibSNIaYieldsTableModelID_S13_N100Z001, 
    CELibSNIaYieldsTableModelID_S13_N100ZDepend, 
    CELibSNIaYieldsTableModelID_S13_Number,
};

/*!
 * @enum CELib SN Ia yields table model IDs for Travaglio et al. 2004
 * SN Ia yields table model IDs.
 */
enum {
    CELibSNIaYieldsTableModelID_T04_W7,      
    CELibSNIaYieldsTableModelID_T04_c32d512, 
    CELibSNIaYieldsTableModelID_T04_c33d256a,
    CELibSNIaYieldsTableModelID_T04_c33d256b,
    CELibSNIaYieldsTableModelID_T04_b53d256, 
    CELibSNIaYieldsTableModelID_T04_b303d768,
    CELibSNIaYieldsTableModelID_T04_Number,
};


/*!
 * @enum CELib AGB yields IDs
 * AGB yields IDs.
 */
enum {
    CELibAGBYieldsTableID_K10,    // Karakas 2010
    CELibAGBYieldsTableID_K10D14, // Karakas 2010 & Doherty et al. 2014, default
    CELibAGBYieldsTableID_vdHG97, // van den Hoek & Groenewgen 1997
    CELibAGBYieldsTableID_Number,
};

/*!
 * @enum CELib AGB yields IDs for ``extremely low metal''.
 * AGB yields IDs for extremely low metal.
 */
enum {
    CELibAGBZ0YieldsTableID_CL08,    // Campbell & Lattanzio 2008
    CELibAGBZ0YieldsTableID_CL08G13, // Campbell & Lattanzio 2008 & Gil-Pons et al. 2013, default
    CELibAGBZ0YieldsTableID_Number,
};

/*!
 * @enum CELib NSM yields IDs
 * NSM yields IDs.
 */
enum {
    CELibNSMYieldsTableID_W14,    // Wanajo et al. 2014, default
    CELibNSMYieldsTableID_Number,
};


/*!
 * @enum CELib feedback types
 * Feedback types.
 */
enum {
    CELibFeedbackType_SNII,
    CELibFeedbackType_SNIa,
    CELibFeedbackType_AGB,
    CELibFeedbackType_NSM,
    CELibFeedbackType_Number,
};

/*!
 * @enum CELibYields
 * Elements used in CELib.
 */
enum {
    CELibYield_H,      //! Hydrogen
    CELibYield_He,     //! Helium
    CELibYield_C,      //! Carbon
    CELibYield_N,      //! Nitrogen
    CELibYield_O,      //! Oxygen
    CELibYield_Ne,     //! Neon
    CELibYield_Mg,     //! Magnesium
    CELibYield_Si,     //! Silicon
    CELibYield_S,      //! Sulfur
    CELibYield_Ca,     //! Calcium
    CELibYield_Fe,     //! Iron
    CELibYield_Ni,     //! Nickel
    CELibYield_Eu,     //! Europium
    CELibYield_Number, //! Total number of elements.
};


/*!
 * @enum CELib SNIa rate model IDs
 * SNIa rate model IDs.
 */
enum {
    CELibSNIaRateModelID_GR83,     // Greggio & Renzini 1983
    CELibSNIaRateModelID_PowerLaw, // Power law  default
    CELibSNIaRateModelID_PowerLawV13, //  Power law, Vogelsberger et al. 2013
    CELibSNIaRateModelID_Number,
};

/*!
 * @enum CELib AGB rate model IDs
 * SNIa AGB model IDs.
 */
enum {
    CELibAGBRateModelID_LinearBin, // LinearBin
    CELibAGBRateModelID_LogBin,    // LogBin
    CELibAGBRateModelID_Number,
};

/*!
 * @enum CELib NSM rate model IDs
 * NSM rate model IDs.
 */
enum {
    CELibNSMRateModelID_PowerLaw, // Power law, default
    CELibNSMRateModelID_Number,
};

/*! @struct CELibStructRunParameters
 * This structure contains all control parameters of CELib.
 */
struct CELibStructRunParameters{
    bool TestMode;                  //!< Test mode flag.
    int IntegrationSteps;           //!< Number of integration steps
    int IMFType;                    //!< Type of IMF.
    int SolarAbundancePatternType;  //!< Type of solar abundance pattern.
    int LifeTimeType;               //!< Type of stellar life time.

    /* Parameters for Type II SNe */
    int SNIIYieldsTableID;           //!< Yields table ID. 0=Portinari+1998, 1=Nomoto+2013
    int SNIIYieldsModificationP98;   //!< Ad hoc modifications for the Portinari+1998 yields table.  If 1, 0.5*C, 2*Mg, 0.5*Fe
    double SNIIHyperNovaFraction;   //!< Blending ratio of hyper novae. This value should be 0-1. This flag works when the yields table of Nomoto+2013 is used.
    double SNIIEnergy;              //!< Energy of one type II SN.
    double SNIIEnergyPerMass;       //!< Type II SN energy released from a 1Msun SSP particle.
    double SNIINumberPerMass;       //!< Number of type II SNe in a 1Msun SSP particle.
    double SNIIEnergyEfficiency;    //!< The efficiency of type II SNe. Usually set to unity.
    double SNIIUpperMass;           //!<
    double SNIILowerMass;           //!<
    int SNIIExplosionTimeFollowYields; // 0: Use SNIILower-SNIIUpper, 1: Use yields min/max
    double SNIIUnitConvert;         //!< From simulation mass unit to sular mass.

    /* Parameters for Type Ia SNe */
    int SNIaType;         //!< This flag controls the adopted DTD in the Ia model.
    int SNIaYieldsTableID;      //!< 0: Iwamoto+1999, 1: Maeda+2010, 2: Seitenzahl+2013, 3: Travaglio+2004
    int SNIaYieldsTableModelID; //!< For SNIaYieldTableID=0, 0: W7, 1: W70, 2: WDD1, 3: WDD2, 4: WDD3, 5: CDD1, 6: CDD2
                               // For SNIaYieldTableID=1, 0: W7, 1: C-DEF, 2: C-DDT, 3: O-DDT
                               // For SNIaYieldTableID=2, 0: N1, 1: N3, 2: N5, 3: N10, 4: N20, 5: N40, 6: N100H
                               // 7: N100, 8: N100L, 9: N150, 10: N200, 11: N300C, 12: N1600, 13: N1600C,
                               // 14: N100_Z0.5, 15: N100_Z0.1, 16: N100_Z0.01, 17: Metal dependent
                               // For SNIaYieldTableID=3, 0: W7, 1: c32d512, 2: c33d256a, 3: c33d256b 
                               // 4: b53d256, 5: b303d768 
    int SNIaNassociation;  // The size of SNIa association.
    double SNIaUpperMassBinary;    // The maximum total mass of a binary sistem.
    double SNIaLowerMassBinary;    // The minimum total mass of a binary sistem.
    double SNIaUpperMass;
    double SNIaLowerMass;
    double SNIaEnergy;
    double SNIaEnergyPerMass;
    double SNIaEnergyEfficiency;

    /* Parameters for AGB Mass Loss */
    int AGBYieldsTableID;    // 0: Karakas 2010, 1: Karakas 2010 + Doherty+2014
    int AGBBinNumber;       // Number of bin for AGB 
    int AGBBinType;         // If 1, the log bin is adopted. Otherwise, the linear bin is used.
    double AGBBinTimeInterval; 
    double AGBBinUpperAge;
    double AGBBinLowerAge;
    double AGBUpperMass;
    double AGBLowerMass;

    /* Parameters for NSM Mass Loss */
    int NSMYieldsTableID; // 0: Wanajo+2014
    int NSMType; // This flag controls the adopted DTD in the NSM model.
    int NSMNassociation;  // The size of NSM association.
    double NSMNumberPerMass; // Event number per Msun.
    double NSMUpperMass; // Upper mass for NS formation.
    double NSMLowerMass; // Lower mass for NS formation.
    double NSMFraction;  // NS binary formation rate.
    double NSMDTDNormalization;  // Normalization for DTD.
    double NSMDTDPowerLawIndex;  // the power index for DTD.
    double NSMDTDOffsetForPower;  // Time offset for the power law type DTD.
    double NSMEnergy;
    double NSMEnergyPerMass;
    double NSMEnergyEfficiency;


    /* Parameters for PopIII */
    int PopIIIIMF;             //!< PopIII yields flag. If 1, CELib adopts the Susa IMF. Otherwith, CELib adopts the ordinary IMF for extremely low metallicity regime.
    int PopIIISNe;             //!< Pop III type II SNe flag. 
                               //!< This flag works when the Nomoto+2013 yields table is used.
    int PopIIIAGB;             //!< Pop III AGB flag
    int PopIIIAGBYieldsTableID;  //!< Pop III AGB yields table ID
    double PopIIIMetallicity; 
    int PopIIILifeTime; 
};

extern struct CELibStructRunParameters CELibRunParameters;


/*! @struct CELibStructFeedbackOutput
 * This structure contains all quantities evaluated by feedback functions.
 */
struct CELibStructFeedbackOutput{
    double Energy;
    double EjectaMass;
    double RemnantMass;
    double Elements[CELibYield_Number];
};

/*! @struct CELibStructFeedbackInput
 * This structure contains all parameters used in feedback functions.
 */
struct CELibStructFeedbackInput{
    double Mass;                 // SSP particle mass in simulation unit.
    double Metallicity;          // SSP particle Z
    double MassConversionFactor; // A factor to convert Elements[] from the simulation mass unit to Msun.
    double *Elements;            // SSP particle elements in simulation unit.
    int    Count;                // Counter
    int    noPopIII;             // No Pop III flag. If 1, this function returns non Pop III yields 
                                 //   even thought the metallicity is lower than Z_popIII
};


/*! @struct CELibStructNextEventTimeInput
 * This structure contains all parameters used to evaluate next event time.
 */
struct CELibStructNextEventTimeInput{
    double R;                    // A random number in [0,1).
    double InitialMass_in_Msun;  // SSP particle mass
    double Metallicity;          // SSP particle Z
    int    Count;                // Counter
    int    noPopIII;             // No Pop III flag. If 1, this function returns the lifetime of 
                                 //   no Pop III stars even thought the metallicity is lower than Z_popIII
};

/*! @struct CELibStructFeedbackStarbyStarOutput
 * This structure contains all quantities evaluated by the star-by-star feedback functions.
 */
struct CELibStructFeedbackStarbyStarOutput{
    double Energy;             // Energy in erg.
    double EjectaMass;         // Total ejecta mass in Msun.
    double RemnantMass;        // Remnant mass in Msun.
    double Elements[CELibYield_Number]; // Return mass to the ISM (Y, not y) in Msun.
};

/*! @struct CELibStructFeedbackStarbyStarInput
 * This structure contains all parameters used in the star-by-star feedback.
 */
struct CELibStructFeedbackStarbyStarInput{
    double Mass;                 // Mass of the star in simulation unit.
    double Metallicity;          // Metallicity of the star Z
    double MassConversionFactor; // A factor to convert Elements[] from the simulation mass unit to Msun.
    double *Elements;            // Star particle's elements composition in simulation unit.
    int    noPopIII;             // No Pop III flag. If 1, this function returns no Pop III yields 
                                 //   even though the metallicity is lower than Z_popIII
};


/*! @struct CELibStructNextEventTimeStarbyStarInput
 * This structure contains all parameters used to evaluate next event time for
 * the star-by-star case.
 */
struct CELibStructNextEventTimeStarbyStarInput{
    double InitialMass_in_Msun;  // Mass of the star
    double Metallicity;          // Metallicity of the star, Z
    int    noPopIII;             // No Pop III flag. If 1, this function returns the lifetime of 
                                 //   no Pop III stars even thought the metallicity is lower than Z_popIII
};

//_______Functions in RunParameters.c
void CELibSaveRunParameter(void);
void CELibLoadRunParameter(void);
void CELibSetRunParameterTestMode(const bool TestMode);
void CELibSetRunParameterIntegrationSteps(const int IntegrationSteps);
void CELibSetRunParameterIMFType(const int IMFType);
void CELibSetRunParameterSolarAbundancePatternType(const int SolarAbundancePatternType);
void CELibSetRunParameterLifeTimeType(const int LifeTimeType);
void CELibSetRunParameterSNIIYieldsTableID(const int SNIIYieldsTableID);
void CELibSetRunParameterSNIIYieldsModificationP98(const int SNIIYieldsModificationP98);
void CELibSetRunParameterSNIIRange(const double SNIIUpperMass_in_SolarMass, const double SNIILowerMass_in_SolarMass);
void CELibSetRunParameterSNIIHyperNovaFraction(const double SNIIHyperNovaFraction);
void CELibSetRunParameterSNIIExplosionTimeFollowYields(const int SNIIExplosionTimeFollowYields);
void CELibSetRunParameterSNIaType(const int SNIaType);
void CELibSetRunParameterSNIaYieldsTableID(const int SNIaYieldsTableID);
void CELibSetRunParameterSNIaYieldsTableModelID(const int SNIaYieldsTableModelID);
void CELibSetRunParameterSNIaRange(const double SNIaUpperMass_in_SolarMass, const double SNIaLowerMass_in_SolarMass);
void CELibSetRunParameterSNIaUpperMassBinary(const double SNIaUpperMassBinary);
void CELibSetRunParameterSNIaLowerMassBinary(const double SNIaLowerMassBinary);
void CELibSetRunParameterSNIaNassociation(const int Nassociation);
void CELibSetRunParameterAGBYieldsTableID(const int AGBYieldsTableID);
void CELibSetRunParameterAGBRange(const double AGBUpperMass_in_SolarMass, const double AGBLowerMass_in_SolarMass);
void CELibSetRunParameterAGBBinUpperAge(const int AGBBinUpperAge);
void CELibSetRunParameterAGBBinLowerAge(const int AGBBinLowerAge);
void CELibSetRunParameterAGBBinNumber(const int AGBBinNumber);
void CELibSetRunParameterAGBBinType(const int AGBBinType);
void CELibSetRunParameterAGBBinTimeInterval(const double AGBBinTimeInterval);
void CELibSetRunParameterNSMYieldsTableID(const int NSMYieldsTableID);
void CELibSetRunParameterNSMType(const int NSMType);
void CELibSetRunParameterNSMDTDPowerLawIndex(const double NSMDTDPowerLawIndex);
void CELibSetRunParameterNSMDTDOffsetForPower(const double NSMDTDOffsetForPower);
void CELibSetRunParameterNSMNassociation(const int NSMNassociation);
void CELibSetRunParameterPopIIIIMF(const int PopIIIIMF);
void CELibSetRunParameterPopIIISNe(const int PopIIISNe);
void CELibSetRunParameterPopIIIAGB(const int PopIIIAGB);
void CELibSetRunParameterPopIIIAGBYieldsTableID(const int PopIIIAGBYieldsTableID);
void CELibSetRunParameterPopIIIMetallicity(const double PopIIIMetallicity);
void CELibSetRunParameterPopIIILifeTime(const int PopIIILifeTime);

bool CELibGetRunParameterTestMode(void);
int CELibGetRunParameterIntegrationSteps(void);
int CELibGetRunParameterIMFType(void);
void CELibGetRunParameterIMFName(char *IMFName);
int CELibGetRunParameterSolarAbundancePatternType(void);
int CELibGetRunParameterLifeTimeType(void);
int CELibGetRunParameterSNIIYieldsModificationP98(void);
double CELibGetRunParameterSNIIUpperMass(void);
double CELibSetRunParameterSNIILowerMass(void);
int CELibGetRunParameterSNIIYieldsTableID(void);
double CELibGetRunParameterSNIIHyperNovaFraction(void);
int CELibSetRunParameterSNIIExplosionTime(void);
int CELibGetRunParameterSNIaYieldsTableID(void);
int CELibGetRunParameterSNIaYieldsTableModelID(void);
int CELibGetRunParameterSNIaType(void);
double CELibGetRunParameterSNIaUpperMass(void);
double CELibGetRunParameterSNIaLowerMass(void);
double CELibGetRunParameterSNIaUpperMassBinary(void);
double CELibGetRunParameterSNIaLowerMassBinary(void);
int CELibGetRunParameterSNIaNassociation(void);
int CELibGetRunParameterAGBYieldsTableID(void);
double CELibGetRunParameterAGBUpperMass(void);
double CELibGetRunParameterAGBLowerMass(void);
int CELibGetRunParameterAGBBinNumber(void);
int CELibGetRunParameterAGBBinType(void);
double CELibGetRunParameterAGBBinTimeInterval(void);
int CELibGetRunParameterNSMYieldsTableID(void);
int CELibGetRunParameterNSMType(void);
double CELibGetRunParameterNSMDTDPowerLawIndex(void);
double CELibGetRunParameterNSMDTDOffsetForPower(void);
int CELibGetRunParameterNSMNassociation(void);
int CELibGetRunParameterPopIIIIMF(void);
int CELibGetRunParameterPopIIISNe(void);
int CELibGetRunParameterPopIIIAGB(void);
int CELibGetRunParameterPopIIIAGBYieldsTableID(void);
double CELibGetRunParameterPopIIIMetallicity(void);
int CELibGetRunParameterPopIIILifeTime(void);
void CELibShowCurrentRunParameters(void);


//_______Functions in WriteIMF.c
void CELibWriteIMFData(const int Type, const char OutputDir[]);
void CELibWriteIMFCumlative(const char OutputDir[]);


//_______Functions for stellar lifetime in LifeTime.c
void CELibInitLifeTime(void);
extern double (*CELibGetLifeTimeofStar)(const double, const double);
extern double (*CELibGetDyingStellarMass)(const double, const double);
double CELibGetLifeTimeofStarTable(const double Mass, const double Metallicity);
double CELibGetDyingStellarMassTable(const double LifeTime, const double Metallicity);
double CELibGetLifeTimeofStarLSF(const double Mass, const double Metallicity);
double CELibGetDyingStellarMassLSF(const double LifeTime, const double Metallicity);
int CELibGetLifeTimeTableSizeMetallicity(void);


//_______Functions for abundance patterns in InitialMetallicity.c
void CELibSetPrimordialMetallicity(const double Mass, double Elements[restrict]);
void CELibSetMetallicityWithSolarAbundancePattern(const double Mass, double Elements[restrict], const double Metallicity);
void CELibSetSolarMetallicity(const double Mass, double Elements[restrict]);
double CELibGetMetalFractionForSolarChemicalComposision(void);
double CELibGetElementWeight(const int ElementIndex);
double CELibGetElementNumberDensitySolar(const int ElementIndex);
int CELibGetSolarAbundancePatternType(void);
void CELibSetSolarAbundancePatternType(const int ID);


//_______Function for SNII in SNIIRate.c
double CELibGetSNIIExplosionTime(const double Rate, const double Metallicity);


//_______Functions for SNIa in SNIaRate.c
void CELibInitSNIaRate(void);
void CELibWriteSNIaRate(char OutDir[]);
double CELibGetSNIaExplosionTime(const double Rate, const double Metallicity, const double InitialMass, const int Count);
double CELibGetSNIaIntegratedRateGreggioRenzini(const double Age, const double Metallicity);
double CELibGetSNIaIntegratedRatePowerLaw(const double Age);


//_______Functions for NSM in NSMRate.c
void CELibInitNSMRate(void);
void CELibWriteNSMRate(char OutDir[]);
double CELibGetNSMFeedbackTime(const double Rate, const double Metallicity, const double InitialMass, const int Count);


//_______Functions in Info.c
void CELibShowVersion(void);
void CELibShowCurrentStatus(void);


//_______Functions in SNIIYields.c
extern double CELibSetSNIIEmptyElements[CELibYield_Number];
struct CELibStructFeedbackOutput CELibGetSNIIFeedback(struct CELibStructFeedbackInput Input);
void CELibGetSNIIIntegratedYields(const double Metallicity, double Elements[]);
double CELibGetSNIIYieldsMaxExplosionMass(const int IndexMetal);
double CELibGetSNIIYieldsMinExplosionMass(const int IndexMetal);
struct CELibStructFeedbackOutput CELibGetSNIIYieldsIntegratedInGivenMassRange(const int MetalID, const double LowerMass_in_SolarMass, const double UpperMass_in_SolarMass); 
struct CELibStructFeedbackOutput CELibGetSNIIYieldsIntegratedInGivenMassRangeZ(const double Metallicity, const double LowerMass_in_SolarMass, const double UpperMass_in_SolarMass);
double CELibGetSNIIIntegratedEnergy(const double Metallicity);
struct CELibStructFeedbackStarbyStarOutput CELibGetSNIIFeedbackStarbyStar(struct CELibStructFeedbackStarbyStarInput Input);


//_______Functions in SNIaYields.c
//#define CELibSNIaYieldTableModelNumber (7)
struct CELibStructFeedbackOutput CELibGetSNIaFeedback(struct CELibStructFeedbackInput Input);
void CELibGetSNIaCurrentModelName(char *ModelName);
void CELibGetSNIaCurrentYieldsModelName(char *ModelName);


//_______Functions in AGBMassLoss.c
double CELibGetAGBFeedbackTime(const double Rate, const double Metallicity, const int Count);
double CELibGetAGBYieldsReturnMassInGivenTimeRange(const double Time_start_in_year, const double Time_end_in_year, const double Metallicity);
struct CELibStructFeedbackOutput CELibGetAGBFeedback(struct CELibStructFeedbackInput Input);
struct CELibStructFeedbackOutput CELibGetAGBYieldsIntegratedInGivenMassRange(const int MetalID, const double LowerMass_in_SolarMass, const double UpperMass_in_SolarMass);
struct CELibStructFeedbackOutput CELibGetAGBYieldsIntegratedInGivenMassRangeZ(const double Metallicity, const double LowerMass_in_SolarMass, const double UpperMass_in_SolarMass);


//_______Functions in NSMYields.c
struct CELibStructFeedbackOutput CELibGetNSMFeedback(struct CELibStructFeedbackInput Input);


//_______Functions in UnifiedAPIs.c
struct CELibStructFeedbackOutput CELibGetFeedback(struct CELibStructFeedbackInput Input, const int Type);
double CELibGetNextEventTime(struct CELibStructNextEventTimeInput Input, const int Type);
struct CELibStructFeedbackStarbyStarOutput CELibGetFeedbackStarbyStar(struct CELibStructFeedbackStarbyStarInput Input, const int Type);
double CELibGetNextEventTimeStarbyStar(struct CELibStructNextEventTimeStarbyStarInput Input, const int Type);


//_______Functions in YieldsInfo.c
void CELibInitYieldNamesIDs(void);
int CELibGetSNIIYieldElementID(const char *name);
int CELibGetSNIaYieldElementID(const char *name);
int CELibGetAGBMassLossElementID(const char *name);
int CELibGetNSMYieldElementID(const char *name);
int CELibGetYieldElementID(const char *name);
void CELibGetSNIIYieldElementName(const int ID, char *name);
void CELibGetSNIaYieldElementName(const int ID, char *name);
void CELibGetAGBMassLossElementName(const int ID, char *name);
void CELibGetNSMYieldElementName(const int ID, char *name);
void CELibGetYieldElementName(const int ID, char *name);
void CELibShowAllElementName(void);
void CELibShowAllElementID(void);

#ifdef __cplusplus
}
#endif
