`HMMoce`: an `R` package for improved analysis of marine animal movement
data using hidden Markov models
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

**Camrin D. Braun<sup>1,2</sup>\*, Benjamin Galuardi<sup>3,4</sup>, Paul
Gatti<sup>5</sup>, [Martin C. Arostegui](https://www.whoi.edu/profile/martin.arostegui/)<sup>2</sup>, Simon R. Thorrold<sup>2</sup>**

1.  Massachusetts Institute of Technology-Woods Hole Oceanographic
    Institution Joint Program in Oceanography/Applied Ocean Science and
    Engineering, Cambridge, MA 02139
2.  Biology Department, Woods Hole Oceanographic Institution, Woods
    Hole, MA 02543
3.  School of Marine Science and Technology, University of Massachusetts
    Dartmouth, Fairhaven, MA 02719
4.  Greater Atlantic Regional Fisheries Office, National Marine
    Fisheries Service, National Oceanic and Atmospheric Administration,
    Gloucester, MA 01930
5.  Centre for Fisheries Ecosystems Research, Fisheries and Marine
    Institute of Memorial University of Newfoundland

## Update 2026-11

Version 2.0.0 
Major modernization, accuracy, and performance updates:

- **Modern Spatial Stack (`terra` & `sf`):** Replaced retired spatial packages (`raster`, `rgdal`, `rgeos`, `maptools`) with `terra` (`SpatRaster`) and `sf` across grid setup, resampling, bathymetry/environmental masking, contour generation (`getCtr_sf`), and all likelihood functions.
- **GLORYS Reanalysis & Tag Ingestion:** Added/updated 3D GLORYS temperature-depth profile (`calc.glorys`, `calc.glorys.par`), GLORYS ocean heat content (`calc.ohc.glorys`, `calc.ohc.glorys.par`), GLORYS surface SST (`calc.sst.par.glorys`), combined PDT likelihood helpers (`make.L.pdt`), and Microwave Telemetry archival XLS-to-eTUFF conversion (`microwave_to_etuff`).
- **HMM Filter, Smoother & Likelihood Refinements:** Improved endpoint anchoring, backward smoother normalization, and convolution speed in `hmm.filter` and `hmm.smoother`; updated multiplicative daily likelihood combination (`make.L`, `make.L.mult`) with empty-layer, light, and bathymetry guards; added fast/`RTMB` parameter optimization support in `opt.params`; and fixed graphics device handling in `plotHMM`.
- **Fast Posterior Viterbi (`hmm.posteriorviterbi`):** Vectorized pure-R Viterbi track estimation (~20–30x faster with identical output) using separable 1D max-plus passes and active bounding-box shifts across 1- and 2-state (`'separate'` and `'joint'`) modes.
- **Testing & Package Maintenance:** Standardized `R/*.R` filenames, refreshed `roxygen2` documentation, and added a `testthat` unit test suite covering grids, kernels, likelihoods, `bin_TempTS`, `hmm.filter`/`hmm.smoother`, `opt.params`, and `hmm.posteriorviterbi`.

## Update 2020-05-13

All the tweaks and improvements made over the last 1-2 years have been
folded into v1.2 which is now on the master branch. The next version
will likely be v2.0 and will lack backwards compatibility with previous
versions. The new version is currently under development in the -dev
branch and pull requests are welcome\!

Many thanks to Paul Gatti for catching many small bugs and providing
several (very\!) useful new functions for what will be the next version.
Check out the -dev branch to learn more.

## Package Summary

Electronic tagging of marine animals is common throughout the world
oceans. Many of these studies have deployed archival tags that rely on
light levels and sea-surface temperatures to retrospectively track
movements of tagged animals. However, methodological issues associated
with light-level geolocation have constrained meaningful inference to
species where it is possible to accurately estimate time of sunrise and
sunset. Most studies have largely disregarded the oceanographic profiles
collected by the tag as a potential way of refining light-level
geolocation estimates provided by electronic tags.

Open-source oceanographic measurements and outputs from high-resolution
models are increasingly available and accessible. We integrated
temperature and depth profiles recorded by electronic tags, with
empirical data and model outputs, to construct likelihoods and improve
geolocation estimates for marine animals using an existing, but
modified, state-space hidden Markov model (HMM). Our model (`HMMoce`)
exhibited as much as 6-fold improvement in pointwise error as compared
to traditional light-level geolocation approaches and produced the
lowest mean error in 3 of 4 cases when compared to the state-of-the-art
tag manufacturer’s HMM (GPE3). `HMMoce` contained behavior
state-switching capability not found in other comparable methods. The
use of profile-based likelihood estimates proved useful when we removed
data to emulate data returned from species that yield poor quality light
data. The results demonstrated the general applicability of the `HMMoce`
model to marine animals, particularly those that do not frequent surface
waters during crepuscular periods. Our model is available as an
open-source `R` package, `HMMoce`, that uses a state-space HMM approach
and leverages available tag and oceanographic data to improve position
estimates derived from electronic tags.

## Package Citation

Braun, C. D., Galuardi, B., and Thorrold, S. R. (2018). HMMoce: An R
package for improved geolocation of archival-tagged fishes using a
hidden Markov method. Methods in Ecology and Evolution, 9(5), 1212-1220.

## Package Structure

The package is structured as follows: \* Load the relevant tag data and
establish a study area of interest. \* Get the environmental data to
base the likelihood calculations on. \* Calculate the desired
likelihoods (e.g. depth-temperature profiles, SST, etc) \* Estimate
parameters and run the model. Results are written out along the way. \*
Perform model checking and choose a final model.

## Installation Instructions

`HMMoce` can be installed from CRAN from within `R` using
`install.packages('HMMoce')`. To get the latest developments, get it
from GitHub using `devtools::install_github('camrinbraun/HMMoce',
ref='dev')`

## Examples

For an example use of the package, please see the vignette using
`vignette('HMMoce')`.

[![Travis-CI Build
Status](https://travis-ci.org/camrinbraun/HMMoce.svg?branch=master)](https://travis-ci.org/camrinbraun/HMMoce)
