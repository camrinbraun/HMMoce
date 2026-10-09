#' HMMoce: an R package for improved analysis of marine animal movement data 
#' using hidden Markov models
#' 
#' The HMMoce package provides a workflow for leveraging all available data from
#' tags deployed on marine animals to estimate movements, space use, and 
#' behavior. Marine animals, mostly fish, are notoriously difficult to track 
#' with electronic tags because these devices require occupation of the 
#' surface-air interface to record satellite-based geolocations or occupation of
#' the photic zone to collect light levels to estimate position. It is common 
#' among fishes to avoid the photic zone during daylight hours thus rendering 
#' this geolocation approach useless. In the HMMoce package, we leverage all 
#' tag-based data streams, like depth-temperature profiles, in conjunction with 
#' whatever traditional geolocation data (e.g. light) is available to calculate 
#' the most probable movements of the tagged animal. This is performed in a 
#' hidden Markov framework originally developed by Pedersen et al. 2008.
#' 
#' @references 
#'  \itemize{
#'   \item Pedersen MW, Righton D, Thygesen UH, et al. (2008) Geolocation of
#'   North Sea cod (Gadus morhua) using hidden Markov models and behavioural
#'   switching. Can J Fish Aquat Sci 65:2367-2377.
#'   
#'   \item Pedersen MW, Patterson TA, Thygesen UH, Madsen H (2011) Estimating
#'   animal behavior and residency from movement data. Oikos 120:1281-1290.
#'   
#'   \item Woillez M, Fablet R, Ngo TT, et al. (2016) A 
#'   HMM-based model to geolocate pelagic fish from high-resolution individual 
#'   temperature and depth histories: European sea bass as a case study. Ecol 
#'   Modell 321:10-22.
#'  }
#'   
#' @docType package
#' @name HMMoce

"_PACKAGE"
#> [1] "_PACKAGE"