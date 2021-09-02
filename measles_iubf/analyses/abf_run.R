library("spatPomp")
library("doRNG")
library("doParallel")

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
 run_level <- 1
} else run_level <- as.numeric(args)

out_files_dir <- paste0("../generated/abf_outs_",run_level,"/")
if(!dir.exists(out_files_dir)) dir.create(out_files_dir)

cores <-  as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset=NA))
if(is.na(cores) & run_level >= 2) cores <- 36
if(is.na(cores)) cores <- detectCores()
registerDoParallel(cores)
registerDoRNG(1598260027L)

source('../analyses/iubf_measles_sim.R')
set.seed(1)
U<-20;N<-2*26
m <- measles_sim_w_params(U=U,N=N)
theta <- coef(m)

if(run_level==1){ # personal machine testing
  abf_nbhd <- function(object, time, unit) {
    nbhd_list <- list()
    if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
    if(time>2) nbhd_list <- c(nbhd_list, list(c(unit, time-2)))
    # if(unit==1 && time > 1) nbhd_list <- c(nbhd_list, list(c(2, time-1)))
    # if(unit==2) nbhd_list <- c(nbhd_list, list(c(1, time)))
    # if(unit==3) nbhd_list <- c(nbhd_list, list(c(1, time)))
    # if(unit==4) nbhd_list <- c(nbhd_list, list(c(3, time)))
    # if(unit==5) nbhd_list <- c(nbhd_list, list(c(4, time)))
    # if(unit==6) nbhd_list <- c(nbhd_list, list(c(1, time)))
    # if(unit==6) nbhd_list <- c(nbhd_list, list(c(4, time)))
    # if(unit==7) nbhd_list <- c(nbhd_list, list(c(1, time)))
    # if(unit==8) nbhd_list <- c(nbhd_list, list(c(1, time)))
    # if(unit==9) nbhd_list <- c(nbhd_list, list(c(1, time)))
    # if(unit==10) nbhd_list <- c(nbhd_list, list(c(5, time)))

    return(nbhd_list)
  }
  abf_nabf <- 1
  abf_nrep <- 35000
  abf_np <- 1
} else if(run_level==2){ # first draft results
  abf_nbhd <- function(object, time, unit) {
    nbhd_list <- list()
    if(unit>1) nbhd_list <- c(nbhd_list, list(c(unit-1, time)))
    if(unit>2) nbhd_list <- c(nbhd_list, list(c(unit-2, time)))
    if(unit>3) nbhd_list <- c(nbhd_list, list(c(unit-3, time)))
    if(unit>4) nbhd_list <- c(nbhd_list, list(c(unit-4, time)))
    if(unit>5) nbhd_list <- c(nbhd_list, list(c(unit-5, time)))
    return(nbhd_list)
  }
  abf_nabf <- 7
  abf_nrep <- 100
  abf_np <- 40
}
eval_params <- theta
eval_params[c("g","R0", "amplitude", "gamma", "sigma", "sigmaSE", "psi")] <- c(400,
                                                                               29.85374,
                                                                               0.5089659,
                                                                               50.80411,
                                                                               47.19193,
                                                                               0.1336978,
                                                                               0.1681658)
measles_abf2 <- list()
for(i in seq_len(abf_nabf)){
  measles_abf2 <- c(measles_abf2,
                    abf(m,
                        params = eval_params,
                        Nrep = abf_nrep,
                        Np=abf_np,
                        nbhd = abf_nbhd))
  print(sapply(measles_abf2, "slot", "loglik"))
}
eval_params <- theta
eval_params[c("g","R0", "amplitude", "gamma", "sigma", "sigmaSE", "psi")] <- c(200,
                                                                               32.08040,
                                                                               0.5182387,
                                                                               46.00447,
                                                                               58.18090,
                                                                               0.1243386,
                                                                               0.1680632)

measles_abf3 <- list()
for(i in seq_len(abf_nabf)){
  measles_abf3 <- c(measles_abf3,
                    abf(m,
                        params = eval_params,
                        Nrep = abf_nrep,
                        Np=abf_np,
                        nbhd = abf_nbhd))
  print(sapply(measles_abf3, "slot", "loglik"))
}

measles_abf4 <- list()
for(i in seq_len(abf_nabf)){
  measles_abf4 <- c(measles_abf4,
                    abf(m,
                        params = theta,
                        Nrep = abf_nrep,
                        Np=abf_np,
                        nbhd = abf_nbhd))
  print(sapply(measles_abf4, "slot", "loglik"))
}

# for U=6, N=39, the truth (with g=400) has abf likelihood around [1] -1307.803
# for U=3, N=52, with Nrep = 1000, Np = 100, c("g","R0", "amplitude", "gamma", "sigma", "sigmaSE", "psi")] <- c(350,33.06345,0.7208187,37.70308,40.08539,0.1738964,0.1482353) has abf ll around [1] -950.1106 -951.9918
# for U=3, N=52, with Nrep = 1000, Np = 100, c("g","R0", "amplitude", "gamma", "sigma", "sigmaSE", "psi")] <- c(350,33.06345,0.7208187,37.70308,40.08539,0.1738964,0.3) has abf ll around [1] -965.3176 -960.9594 -965.8951
# for U=3, N=52, with Nrep = 100, Np = 100, c("g","R0", "amplitude", "gamma", "sigma", "sigmaSE", "psi")] <- c(350,33.06345,0.7208187,37.70308,40.08539,0.1738964,0.1482353) has abf ll around [1] -965.6622 -959.1156 -975.2946 -994.2319 -971.9970 -967.6285 -967.3767
# for U=3, N=52, with Nrep = 100, Np = 100, c("g","R0", "amplitude", "gamma", "sigma", "sigmaSE", "psi")] <- c(350,33.06345,0.7208187,37.70308,40.08539,0.1738964,0.3) has abf ll around [1] -974.9976 -970.4662 -970.9785 -969.8677 -973.6539 -988.6620 -967.7743

# replicate(10, abf(m,
#     params = coef(m),
#     Nrep = 1000,
#     Np=50,
#     nbhd = function(object, time, unit) {
#      nbhd_list <- list()
#      if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
#      if(time>2) nbhd_list <- c(nbhd_list, list(c(unit, time-2)))
#      return(nbhd_list)
#     })
# ) -> ans
# for U=3, N=4*26: [1] -1877.531 -1874.184 -1872.959 -1873.624 -1878.786 -1875.450 -1877.261 -1874.768 -1871.354 -1871.983




