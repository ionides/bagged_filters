library("spatPomp")
library("doRNG")
library("doParallel")
library("ggplot2")

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
 run_level <- 1
} else run_level <- as.numeric(args)

print(paste("run level: ", run_level))
out_files_dir <- paste0("../generated/abf_slice_",run_level,"/u20/")
if(!dir.exists(out_files_dir)) dir.create(out_files_dir)

cores <-  as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset=NA))

if(is.na(cores) & run_level >= 2) cores <- 36
if(is.na(cores)) cores <- detectCores()
registerDoParallel(cores)
registerDoRNG(1598260027L)

source('../analyses/iubf_measles_sim.R')
set.seed(1)
U<-20;N<-8*26
m <- measles_sim_w_params(U=U,N=N)
theta <- coef(m)
psis <- seq(from=0.1, to=0.4, by = 0.05)
slice_params <- as.matrix(replicate(length(psis),theta))
slice_params['psi',] <- psis

if(run_level==1){ # personal machine testing
  abf_nbhd <- function(object, time, unit) {
    nbhd_list <- list()
    if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
    if(time>2) nbhd_list <- c(nbhd_list, list(c(unit, time-2)))
    return(nbhd_list)
  }
  job_num <- 0
  abf_reps_per_param <- 3
  abf_nrep <- 40000
  abf_np <- 1
} else if(run_level==2){ # first draft results
  abf_nbhd <- function(object, time, unit) {
    nbhd_list <- list()
    if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
    if(time>2) nbhd_list <- c(nbhd_list, list(c(unit, time-2)))
    return(nbhd_list)
  }
  job_num <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
  abf_reps_per_param <- 3
  abf_nrep <- 10000*job_num
  abf_np <- 1
}

print(paste("job_num: ", job_num))

slice_runs <- vector(mode="list", length=length(psis))
stew(file=paste0(out_files_dir,"abf_slice_psi",job_num,".rda"),seed=5981724,{
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

num_jobs <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_COUNT'))
if(num_jobs > 1){
  for(l in seq_len(num_jobs)){
    load(paste0(out_files_dir,"abf_slice_psi",l,".rda"))
    slice_logliks <- vector(mode="list", length=length(psis))
    for(j in seq_len(ncol(slice_params))){
      for(k in seq_len(abf_reps_per_param)){
        slice_logliks[[j]] <- c(slice_logliks[[j]], slice_runs[[j]][[k]]@loglik)
      }
    }
    slice_df <- data.frame(psi=rep(psis, each = abf_reps_per_param),
                           loglik=unlist(slice_logliks))
    ggplot2::ggplot(data = slice_df) +
      ggplot2::geom_point(ggplot2::aes(x = psi, y = loglik)) +
      ggplot2::geom_vline(xintercept = coef(m)[["psi"]], color = "red")
  }
} else {
  slice_logliks <- vector(mode="list", length=length(psis))
  for(i in seq_len(ncol(slice_params))){
    for(j in seq_len(abf_reps_per_param)){
      slice_logliks[[i]] <- c(slice_logliks[[i]], slice_runs[[i]][[j]]@loglik)
    }
  }
  slice_df <- data.frame(psi=rep(psis, each = abf_reps_per_param),
                         loglik=unlist(slice_logliks))
  ggplot2::ggplot(data = slice_df) +
    ggplot2::geom_point(ggplot2::aes(x = psi, y = loglik)) +
    ggplot2::geom_vline(xintercept = coef(m)[["psi"]], color = "red")
}
