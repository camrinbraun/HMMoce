`HMMoce`: an `R` package for improved analysis of marine animal movement data using hidden Markov models
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
**Camrin D Braun, Benjamin Galuardi, Simon R. Thorrold**

Summary
-------

1.  While the number of marine animals being tagged and tracked continues to grow, current satellite tracking techniques largely constrain meaningful inference to largescale movements of surface-dwelling species and are inherently prone to significant error. The paradigm in fish tracking is to use light and sea surface temperatures to estimate animal movements, largely disregarding the 3-D data collected by the tag.
2.  Oceanographic data and high-resolution models are constantly improving and are becoming widely available and accessible. This data can be leveraged, along with 3-D tag data, to improve estimates of animal movements. We use a likelihood approach for these observations and convolve them with a theoretical movement model in a hidden Markov framework that accomodates state-switching. This discretized approach also provides efficient handling of grid-based oceanographic data and likelihood surfaces generated within the package.
3.  We present an open-source `R` package, `HMMoce`, that uses a novel state-space HMM approach to improve position estimates derived from electronic tags using three-dimensional oceanographic data.
4.  We demonstrate `HMMoce` with example blue shark (*Prionace glauca*) data that is included in the package. Our findings illustrate how our software leverages all available tag data, along with oceanographic information, to improve position estimates of tagged marine species.

Installation Instructions
-------------------------

`HMMoce` can be installed from CRAN from within `R` using `install.packages('HMMoce')`.

Examples
--------

For an example use of the package, please see the vignette using `vignette('HMMoce')`.
