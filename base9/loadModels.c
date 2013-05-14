#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "evolve.h"
#include "loadModels.h"
#include "msRgbEvol.h"
#include "gBergMag.h"
#include "wdCooling.h"
#include "gBaraffeMag.h"

// Declared in parent program (mcmc, simCluster, makeCMD)
extern int    verbose;

void loadModels(int needFS, struct cluster *theCluster)

// **************************************************************************************
// last update: 20jul10
// 
// Query the user for model set #.
// 
// Call appropriate function to load that module's models
// 
// Returns the number of the model set.
// **************************************************************************************

{

  //int  modelSet;
  char path[100] = "models/\0";

  //printf("\n Enter path of input models (e.g., models/): ");
  //scanf("%s",path);
  
  printf("Choose a main sequence isochrone set:\n");
  printf("\t0. Girardi\n");
  printf("\t1. Chaboyer-Dotter (w/helium sampling)\n");
  printf("\t2. Yonsei-Yale\n");
  printf("\t3. DSED\n");
  printf("> ");
  scanf("%d",&(theCluster->evoModels.mainSequenceEvol));
  
  if(theCluster->evoModels.mainSequenceEvol < 0 || theCluster->evoModels.mainSequenceEvol > 3) {
      printf("***Error: No models found for main sequence evolution model %d.***\n",theCluster->evoModels.mainSequenceEvol);
      printf("[Exiting...]\n");
      exit(1);
  }
  
  printf("Choose a filter set:\n");
  printf("\t0. Standard (UBVRIJHK)\n");
  printf("\t1. ACS\n");
  printf("\t2. SDSS (ugriz) + 2Mass (JHK)\n");
  printf("> ");
  scanf("%d",&(theCluster->evoModels.filterSet));

  if(theCluster->evoModels.filterSet < 0 || theCluster->evoModels.filterSet > 2) {
      printf("***Error: No models found for filter set %d.***\n",theCluster->evoModels.filterSet);
      printf("[Exiting...]\n");
      exit(1);
  }

  setFilterNames(theCluster->evoModels.filterSet);
  
  printf("Choose a white dwarf filter set:\n");
  printf("\t0. Wood \n");
  printf("\t1. Montgomery \n");
  printf("> ");
  scanf("%d",&(theCluster->evoModels.WDcooling));

  printf("%d\n", theCluster->evoModels.WDcooling);

  if(theCluster->evoModels.WDcooling < 0 || theCluster->evoModels.WDcooling > 1) {
      printf("***Error: No models found for white dwarf filter set %d.***\n",theCluster->evoModels.WDcooling);
      printf("[Exiting...]\n");
      exit(1);
  }

  printf("Choose a white dwarf carbonicity (between 0.0 and 1.0):\n");
  printf("> ");
  scanf("%lf",&(theCluster->carbonicity));

  printf("Choose an initial-final mass relation:\n");
  printf("\t0. Weidemann\n");
  printf("\t1. Williams\n");
  printf("\t2. Salaris Linear\n");
  printf("\t3. Salaris Piecewise Linear\n");
  printf("> ");
  scanf("%d",&(theCluster->evoModels.IFMR));

  printf("Choose a brown dwarf model set:\n");
  printf("\t0. None\n");
  printf("\t1. Baraffe\n");
  printf("> ");
  scanf("%d",&(theCluster->evoModels.brownDwarfEvol));
  //printf("*** %d ***\n",theCluster->evoModels.brownDwarfEvol);
  
  if(theCluster->evoModels.brownDwarfEvol == BARAFFE)
    loadBaraffe(path);

  theCluster->evoModels.WDatm = BERGERON;
  
  printf("\nReading models...\n");

  //setModels(theCluster, modelSet);
  loadMSRgbModels(theCluster, path, needFS);
  loadWDCool(path, theCluster->evoModels.WDcooling);
  loadBergeron(path,  theCluster->evoModels.filterSet);
  
  printf("Models read.\n");

}
