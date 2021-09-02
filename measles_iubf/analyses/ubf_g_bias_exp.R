library("spatPomp")
library("doRNG")
library("doParallel")
library("ggplot2")

# args = commandArgs(trailingOnly=TRUE)
# if(length(args)==0){
#  run_level <- 1
# } else run_level <- as.numeric(args)
# 
# out_files_dir <- paste0("../generated/abf_slice_",run_level,"/")
# if(!dir.exists(out_files_dir)) dir.create(out_files_dir)

cores <-  as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset=NA))
if(is.na(cores) & run_level >= 2) cores <- 36
if(is.na(cores)) cores <- detectCores()
registerDoParallel(cores)
registerDoRNG(1598260027L)

source('../analyses/iubf_measles_sim.R')
set.seed(1)
U<-3;N<-2*26
m <- measles_sim_w_params(U=U,N=N)
theta <- coef(m)
theta_alt <- theta
theta_alt['g'] <- 1100

abf_nbhd <- function(object, time, unit) {
 nbhd_list <- list()
 if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
 if(time>2) nbhd_list <- c(nbhd_list, list(c(unit, time-2)))
 return(nbhd_list)
}
abf_nrep <- 1500
abf_np <- 1
set.seed(1)
a <- replicate(5, abf(m,
                      params = theta,
                      nbhd = abf_nbhd,
                      Nrep = abf_nrep,
                      Np = abf_np)
)
a_alt <- replicate(5, abf(m,
                          params = theta_alt,
                          nbhd = abf_nbhd,
                          Nrep = abf_nrep,
                          Np = abf_np)
)

p <- replicate(5, pfilter(m,
                          params = theta,
                          Np = 20000))
p_alt <- replicate(5, pfilter(m,
                              params = theta_alt,
                              Np = 20000))

g <- replicate(5, girf(m,
                       params = theta,
                       Np = 2000,
                       Nguide = 30,
                       lookahead = 1))
g_alt <- replicate(5, girf(m,
                       params = theta_alt,
                       Np = 2000,
                       Nguide = 30,
                       lookahead = 1))



