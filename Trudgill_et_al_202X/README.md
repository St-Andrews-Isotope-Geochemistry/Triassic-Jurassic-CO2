# Materials and code associated with Trudgill et al., 202X

## Citation
TBC

## Quick Links
- [Figure X](#figure-x)

## Figures
### Figure X
  ![Figure X][figurex]

Some figures are slightly different to those in the manuscript due to the copy editing process.

## How to download
```
git init
git remote add origin https://github.com/St-Andrews-Isotope-Geochemistry/Publications
git fetch Trudgill_et_al_202X
git checkout Trudgill_et_al_202X
git submodule update --recursive --init
```

## Outline
### Calculate_Age
The calculate age script takes values in [TJ_Age_Calibation](./Data/TJ_Age_Calibration.xlsx) and performs a piecewise linear interpolation to get the age of samples in [TJ_d18O_d13C](./Data/TJ_d18O_d13C.xlsx) and [TJ_d11B_pH](./Data/TJ_d11B_pH.xlsx).

### d18O Temperature
Applies multiple calibrations to the &delta;<sup>18</sup>O values in [TJ_d18O_d13C](./Data/TJ_d18O_d13C.xlsx) to calculate temperature. This is interpolated to the same depths as boron data in [TJ_d11B_pH](./Data/TJ_d11B_pH.xlsx).

### Maximum Initial pH
Calculates the maximum initial pH based on individual values for all inputs. This is done using two carbonate chemistry libraries to ensure parity.

### Minimum pH Change
Calculates pH from the &delta;<sup>11</sup>B values in [TJ_d11B_pH](./Data/TJ_d11B_pH.xlsx). Uncertainties are propagated using a latin hypercube technique to estimate the extrema of input parameters which would act to minimise pH change.

### Minimum pH Change Variation
Calculates pH change while varying each input parameter independently, in order to establish the sensitivity of results to a each parameter.

### TJ_CO2
Calculates pH and atmospheric CO<sub>2</sub> concentration from the &delta;<sup>11</sup>B values in [TJ_d11B_pH](./Data/TJ_d11B_pH.xlsx). Uncertainties are propagated using a latin hypercube technique.

### Run_TJ_CO2
Used to produce additional samples from [TJ_CO2](./Code/Analysis/TJ_CO2.m) and save them to [TJ_CO2_Evolutions](./Data/TJ_CO2_Evolutions.csv).

### Analyse_Ensemble
Parses the results in [TJ_CO2_Evolutions](./Data/TJ_CO2_Evolutions.csv) to produce metrics (such as median and standard deviation).

## Dependencies
This repository uses [git LFS](https://git-lfs.github.com/) to store a large data file and figures. If you want to download these files install git lfs as detailed [here](https://git-lfs.github.com/) and then run: `git lfs fetch`

### Submodules
CO2_Systematics performs some CO2 calculations and is found [here](https://github.com/St-Andrews-Isotope-Geochemistry/CO2_Systematics). For this paper the 'script' branch was used, which is a version of csys which has been updated to use equilibrium coefficients for carbonate chemistry from the [MyAMI model](https://github.com/St-Andrews-Isotope-Geochemistry/MyAMI).  

BuCC performs most of the carbonate chemistry calculations for this paper, and is found [here](https://github.com/St-Andrews-Isotope-Geochemistry/BuCC). It is an object oriented implementation of carbonate chemistry calculations.  

Geochemistry_Helpers provides some useful classes to represent geochemical concepts (such as pX, which represents pH, and delta to represent isotope measurements on the delta scale). Sampling and distribution classes are also found in this repository. It can be found [here](https://github.com/St-Andrews-Isotope-Geochemistry/Geochemistry_Helpers).


[figurex]: ./Figures/xxx.png "xxx"
