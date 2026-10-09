#' Estimate parameters
#' 
#' \code{opt.params} is a wrapper for several different parameter optimization routines
#' 
#' @param pars.init is vector of initial parameter set: sigma1 (resident), sigma2 (migrant), p11 switch from behavior 1 to 1 (remain resident), p22 switch from behavior 2 to 2 (remain migrant)
#' @param lower.bounds vector of lower bounds of pars.init
#' @param upper.bounds vector of upper bounds of pars.init
#' @param g grid from \code{\link{setup.grid}}
#' @param L is likelihood array output from \code{make.L}
#' @param alg.opt character indicating which optimization algorithm to use. Options are "optim", "nlminb", "ga", "ga_optim". 
#' @param reformat logical indicating whether to reformat and thus standardize output of the optimization algorithms
#' @param write.results optional. logical indicating whether to write results with \code{saveRDS}.
#' @param hessian optional. logical indicating whether to compute hessian matrix for estimate of parameter standard deviation. Note: usually slows down the computations significantly. Default is FALSE.
#' @param maskL optional. logical indicating whether to apply likelihood mask in \code{hmm.filter}. Default is FALSE.
#' @param optim_method optional. see \code{method} input to \code{optim}. Default is "L-BFGS-B".
#' @param filename optional. if \code{write.results == TRUE}, this should satisfy \code{file} argument to \code{saveRDS}
#' @param max_iter optional. see \code{maxiter} input to \code{GA::ga()}. Default is 100.
#' @param run optional. see \code{run} input to \code{GA::ga()}. Default is 100.
#' @param p_size optional. see \code{popSize} input to \code{GA::ga()}. Default is 100.
#' @param ncores optional. see \code{parallel} input to \code{GA::ga()}. Default is FALSE, therefore not parallelized.
#' @export
#' @author Paul Gatti
#' @author Camrin Braun
#' @examples 
#' \dontrun{
#' pars.init=c(2,.2,.6,.8)
#' lower.bounds=c(0.1,0.001,.1,.1)
#' upper.bounds=c(5,.5,.9,.9)
#' 
#' }

opt.params <- function(pars.init, lower.bounds, upper.bounds, g, L, alg.opt = 'optim', reformat = TRUE, ...){
  
  args <- list(...)
  
  write.results <- ifelse('write.results' %in% names(args), args$write.results, FALSE)
  hessian <- ifelse('hessian' %in% names(args), args$hessian, FALSE)
  maskL <- ifelse('maskL' %in% names(args), args$maskL, FALSE)
  optim_method <- ifelse('optim_method' %in% names(args), args$optim_method, 'L-BFGS-B')
  filename <- ifelse('filename' %in% names(args), args$filename, 'results.RDS')
  #pm <- ifelse('pm' %in% names(args), args$pm, 0.1)
  max_iter <- ifelse('max_iter' %in% names(args), args$max_iter, 100)
  run <- ifelse('run' %in% names(args), args$run, 100)
  p_size <- ifelse('p_size' %in% names(args), args$p_size, 100)
  ncores <- ifelse('ncores' %in% names(args), args$ncores, FALSE)
  
  # optional args include:
    ## hessian = logical. compute hessian matrix for estimate of parameter SD --> usually slow down the computations
    ## maskL logical
    ## optim_method see method in ?optim Default is 'L-BFGS-B'
    ## filename if write.results == TRUE. should satisfy file arg to saveRDS
    ## pm=.1 # probability of mutation # default value of GA
    ## max_iter=100 # maximum number of generations (i.e. iterations)
    ## Run=100 # number of iteration without improvement befor ethe algorithm is stopped
    ## p_size=100 # population size (each individual is a parameter set)
    ## ncores=4 # number of core to use for the parallelization, usually its good to use half of the cores (and let some for the background processes)
  
  #-----------------------------------------------
  # Parameter estimation
  #-----------------------------------------------
  # a. nlminb or optim (gradient based)
  # b. GA genetic algorithm
  #-----------------------------------------------
  
  #-----------------------------------------------
  # a. nlminb or optim
  #-----------------------------------------------
  
  if(alg.opt == 'optim'){
    t0 = Sys.time()
    res <- stats::optim(par = pars.init, 
                    fn = neg.log.lik.fun, 
                    g = g, 
                    L = L, 
                    maskL = maskL, 
                    gr = NULL, 
                    lower = lower.bounds, 
                    upper = upper.bounds, 
                    method = optim_method,
                    hessian = hessian)
    t1 = Sys.time()
    print(t1 - t0) 
    if (write.results) saveRDS(res, file = filename)
    
    ## format optim
    if (reformat) res <- list(par = res$par,
                              value = res$value,
                              convergence = res$convergence,
                              message = res$message)
    
  }
  

  if(alg.opt == 'nlminb'){
    t0 = Sys.time()
    
    res <- tryCatch({stats::nlminb(pars.init, 
                            objective = neg.log.lik.fun,
                            maskL = maskL, 
                            gradient = NULL, 
                            hessian = NULL,
                            g = g,
                            L = L,
                            scale = 1,
                            lower = lower.bounds, 
                            upper = upper.bounds,
                            control = list(iter.max = 100))},
                    error = function(error_message){return(error_message)})
    
    t1 = Sys.time()
    print(t1 - t0) 
    if (write.results) saveRDS(res, file = filename)
    
    ## format nlminb output like optim
    if (reformat) res <- list(par = res$par,
                              value = res$objective,
                              convergence = res$convergence,
                              message = res$message)
    
  }
  
  
  #-----------------------------------------------
  # b. GA genetic algorithm
  #-----------------------------------------------
  # evolutionnary algorithm, completely different approach from gradient based such as optim or nlminb
  # do not need initialization value,
  # although you can still provide a parameter set that might work,
  # and the algorithm will integrate it in its initial population of parameter set
  # IT IS MUCH MORE time/computing consumming BUT also more resilient to local minima
  # has a common feature with SEM, i.e. random sampling (mutations)
  
  if(alg.opt == 'ga' | alg.opt == "ga_optim"){
    
    print('Genetic algorithm with the higher resolution likelihood grid may take a VERY long time to run (several hours).')
    t0 = Sys.time()
    
    res <- GA::ga(type = 'real-valued',
                 fitness = pos.log.lik.fun,
                 g = g,
                 L = L,
                 maskL = maskL,
                 lower = lower.bounds,
                 upper = upper.bounds,
                 popSize = p_size,
                 maxiter = max_iter,
                 run = run,
                 #names = c('sigma1.ncell','sigma2.ncell','pswitch11','pswitch22'),
                 keepBest = T,
                 parallel = ncores
    )
    
    t1 = Sys.time()
    print(t1 - t0)
    if (write.results) saveRDS(res, file = filename)
    
    ## format ga output like optim
    if (reformat) res <- list(par = res@solution,
                              value = res@fitnessValue,
                              convergence = NA,
                              message = NA)

  }
  
  if(alg.opt == "ga_optim"){
    
    print('Refining genetic algorithm optimization results with optim...')
    t0 = Sys.time()
    
    par.inits <- ifelse(reformat, res$par, res@solution)
    
    res <- stats::optim(par = par.inits,
                       fn = neg.log.lik.fun,
                       g = g,
                       L = L,
                       maskL = maskL,
                       gr = NULL, 
                       lower = lower.bounds, 
                       upper = upper.bounds, 
                    method = optim_method,
                    hessian = hessian)
    
    t1 = Sys.time()
    print(t1 - t0) 
    if (write.results) saveRDS(res, file = filename)
    
    ## format optim
    if (reformat) res <- list(par = res$par,
                              value = res$value,
                              convergence = res$convergence,
                              message = res$message)
    
  }
  
  return(res)
}
