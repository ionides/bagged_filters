library("spatPomp")
library("doRNG")
library("doParallel")
library("ggplot2")

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
 run_level <- 1
} else run_level <- as.numeric(args)

out_files_dir <- paste0("../generated/abf_slice_",run_level, "/u5/")
if(!dir.exists(out_files_dir)) dir.create(out_files_dir)

cores <-  as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset=NA))
if(is.na(cores) & run_level >= 2) cores <- 36
if(is.na(cores)) cores <- detectCores()
registerDoParallel(cores)
registerDoRNG(1598260027L)

source('../analyses/iubf_measles_sim.R')
set.seed(1)
U<-5;N<-8*26
m <- measles_sim_w_params(U=U,N=N)
theta <- coef(m)
gs <- seq(from=100, to=1100, by = 100)
slice_params <- as.matrix(replicate(length(gs),theta))
slice_params['g',] <- gs
abf_reps_per_param <- 4
abf_nrep <- 30000
abf_np <- 1
slice_runs <- vector(mode="list", length=length(gs))
abf_nbhd <- function(object, time, unit) {
 nbhd_list <- list()
 if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
 if(time>2) nbhd_list <- c(nbhd_list, list(c(unit, time-2)))
 return(nbhd_list)
}
stew(file=paste0(out_files_dir,"abf_slice_g",".rda"),seed=5981724,{
  for(i in seq_len(ncol(slice_params))){
   slice_runs[[i]] <- vector(mode="list", length = abf_reps_per_param)
   for(j in seq_len(abf_reps_per_param)){
    slice_runs[[i]][[j]] <- abf(m,
                                params = slice_params[,i],
                                Nrep = abf_nrep,
                                Np = abf_np,
                                nbhd = abf_nbhd
                            )
   }
  }
})
slice_logliks <- vector(mode="list", length=length(gs))
for(i in seq_len(ncol(slice_params))){
  for(j in seq_len(abf_reps_per_param)){
    slice_logliks[[i]] <- c(slice_logliks[[i]], slice_runs[[i]][[j]]@loglik)
  }
}
slice_df <- data.frame(g=rep(gs, each = abf_reps_per_param),
                       loglik=unlist(slice_logliks))
ggplot2::ggplot(data = slice_df) +
 ggplot2::geom_point(ggplot2::aes(x = g, y = loglik)) +
  ggplot2::geom_vline(xintercept = coef(m)[["g"]], color = "red")
