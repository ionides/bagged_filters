measles_sim_w_params <- function(U = 40,
                                 N=15*26,
                                 mup=c(alpha=1,
                                       iota=0,
                                       R0=30,
                                       amplitude=0.5,
                                       cohort=0,
                                       gamma=52,
                                       sigma=52,
                                       mu=0.02,
                                       sigmaSE=0.15,
                                       rho=0.5,
                                       psi=0.15,
                                       g=400)){
 measles_sim_U <- U
 measles_uk <- measles(measles_sim_U)
 measles_RPnames <- c("alpha","iota","R0","amplitude","cohort",
                      "gamma","sigma","mu","sigmaSE","rho","psi","g")
 measles_paramnames <- measles_RPnames
 
 measles_params <- rep(NA,length(measles_paramnames))
 names(measles_params) <- measles_paramnames
 measles_unit_params <- mup
 measles_params[measles_RPnames] <- measles_unit_params[measles_RPnames]

 set.seed(34)
 measles_sim <- simulate(measles_uk,params=measles_params)
 measles_sim@data <- measles_sim@data[1:U,1:N]
 time(measles_sim) <- measles_sim@times[1:N]
 
 m_paramnames <- measles_RPnames
 m_params <- measles_params[names(measles_params)%in%m_paramnames]
 coef(measles_sim) <- m_params
 return(measles_sim)
}