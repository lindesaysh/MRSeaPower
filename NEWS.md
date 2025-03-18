# MRSea 1.1-beta (devo)

## Notes

* General:

- Website initiated and under development. Added README for website homepage. 
- removed functions that were duplicated in `MRSea` (`generateNoise`, `rpois.od`).  `MRSea` functions are more up to date so this package is specified as a requirement to run `MRSeaPower`. 
- update to `genChangeData` function.  Allows users to specify the fitted values from the model instead of requiring the model object. 

* Vignettes:

- `UsingMRSeaPower` all code checked. Updates to change to `ggplot` in places. 
- New vignette for imposing changes: `genchangeData`
- New vignette for an introduction to the package: `IntroductiontoMRSeaPower`.  This vignette gives an overview of the power analysis process and possible outputs without the details of coding.

* Other: 

  
## Bug Fixes

* General Code

- fixed error with `powerPlot` function

* Documentation

- fixed documentation error associated with use of percent sign in `genChangeData`, `plot.coverage`, `plot.preds`, `plot.sigdiff`,  `powerSimParallel`

* Other:



