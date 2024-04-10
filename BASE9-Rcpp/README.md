## Local Package Installation
1. Download the latest development archive
2. Unzip the archive
3. `cd` into `base-develop/`
4. Run ./build_local.sh
5. Run R CMD build BASE9
6. Run R CMD check BASE9

## Usage

1. `cd` into `BASE9.Rcheck`
2. Run `R`

Now, from the R prompt:

```R
> install.packages("Rcpp") # If not already installed
> library(BASE9, lib.loc=".")

> initBase("base-models-master/", "PARSEC/PARSEC_2022.model", 1, 1, FALSE) # Note that this needs to be the full path to your models directory, including trailing '/'
> setClusterParameters(8.796, 0.07, 0.0, 0.09, 0.29, 0.5) # These are the reference base9.yaml parameters for the Hyades
> evolve(2.64) # Returns an array containing magnitudes for a 2.64 solar mass star
```

## Interface
### initBase
```C++
void initBase(string modelDir, string msModel, int wdModel, int ifmr, bool distModIsParallax)
```
```R
initBase <- function(modelDir, msModel, wdModel, ifmr, distModIsParallax)
```

Takes a path to the model directory (include trailing ‘/‘), a main sequence model file (relative to the model directory), and the integer equivalents for the WD model and IFMR (see ../conf/base9.yaml).

`distModIsParallax` should be `TRUE` if you intend to use parallax distances. This can only be set by `initBase` and `changeModels` due to limitations in my implementation of the R -> C++ binding.

This function **must** be called before any other function from this package.

### changeModels
```C++
void changeModels(std::string msModel, int wdModel, int ifmr, bool distModIsParallax)
```
```R
changeModels <- function(msModel, wdModel, ifmr, distModIsParallax)
```

Changes the computer models in the same manner as `initBase`. Uses the previously input `modelDir`. This function reloads the models (which can take several seconds) and should be used sparingly.

### setClusterParameters
```C++
void setClusterParameters(double age, double feh, double distMod, double av, double y, double carbonicity);
```
```R
setClusterParameters <- function(age, feh, distMod, av, y, carbonicity)
```

Sets cluster parameters (Θ). This function is more expensive than evolve, so it should be called less often if possible.

This is the most sensible function to call after `initBase`.

### setIFMRParameters
```C++
void setIFMRParameters(double intercept, double slope, double quadCoef)
```
```R
setIFMRParameters <- function(intercept, slope, quadCoef)
```

Sets the IFMR parameters.

### evolve
```C++
std::vector<double> evolve (double mass1, double mass2)
```
```R
evolve <- function(mass1, mass2)
```

Takes a primary and secondary mass and returns magnitudes for all filters available to both the MS and WD models in the order described by `listFilters`.

### listFilters
```C++
std::vector<std::string> listFilters()
```
```R
listFilters <- function()
```

Returns a list of every filter output by `evolve` at its appropriate index.

### getAGBt_zmass
```C++
std::vector<std::string> listFilters()
```
```R
listFilters <- function()
```

Returns the AGB zero-age tip mass for the current isochrone. You must `setClusterParameters` for this to make sense.
