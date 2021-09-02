library("spatPomp")
library("doRNG")
library("doParallel")

dd <- function(){
   rm(list = ls())
}
dd()

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
   run_level <- 1
} else run_level <- as.numeric(args)

cores <-  as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset=NA))
if(is.na(cores)) cores <- detectCores()
registerDoParallel(cores)
registerDoRNG(1598260027L)

if(run_level==1){ # personal machine testing
   U <- 5; N <- 10
   iubf_nbhd <- function(object, time, unit) {
      nbhd_list <- list()
      if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
      if(time>2) nbhd_list <- c(nbhd_list, list(c(unit, time-2)))
      return(nbhd_list)
   }
   iubf_nubf <- 5
   iubf_nrep_per_param <- 5000
   iubf_nparam <- 200
   iubf_prop <- 0.80
   g_lo = 200
   g_hi = 1100
   g_num = 5
   g_nprof = 2
} else if(run_level==2){ # first draft results
   U <- 10; N <- 8*26
   iubf_nbhd <- function(object, time, unit) {
      nbhd_list <- list()
      if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
      if(time>2) nbhd_list <- c(nbhd_list, list(c(unit, time-2)))
      return(nbhd_list)
   }
   iubf_nubf <- 15
   iubf_nrep_per_param <- 35000
   abf_nreps <- 60000
   iubf_nparam <- 250
   iubf_prop <- 0.80
   g_lo = 100
   g_hi = 1100
   g_num = 11
   g_nprof = 3
   all_g = seq(from = g_lo, to = g_hi, length = g_num)
   all_g = append(all_g, 150, after = 1)
   all_g = c(50,all_g)
   abf_nruns = 8
} else if(run_level==3){ # production run
   U <- 20; N <- 8*26
   iubf_nbhd <- function(object, time, unit) {
      nbhd_list <- list()
      if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
      if(time>2) nbhd_list <- c(nbhd_list, list(c(unit, time-2)))
      return(nbhd_list)
   }
   iubf_nubf <- 15
   iubf_nrep_per_param <- 30000
   iubf_nparam <- 250
   iubf_prop <- 0.80
   g_lo = 100
   g_hi = 1100
   g_num = 11
   g_nprof = 3
   all_g = seq(from = g_lo, to = g_hi, length = g_num)
   all_g = c(50,all_g)
   abf_nruns = 8
} else if(run_level==4){ # if need be
   U <- 20; N <- 8*26
   iubf_nbhd <- function(object, time, unit) {
      nbhd_list <- list()
      if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
      if(time>2) nbhd_list <- c(nbhd_list, list(c(unit, time-2)))
      return(nbhd_list)
   }
   iubf_nubf <- 18
   iubf_nrep_per_param <- 10000
   iubf_nparam <- 250
   iubf_prop <- 0.80
   g_lo = 100
   g_hi = 1100
   g_num = 11
   g_nprof = 3
   all_g = seq(from = g_lo, to = g_hi, length = g_num)
   all_g = c(50,all_g)
   abf_nruns = 8
}

source('../analyses/iubf_measles_sim.R')
source('../analyses/iubf_profile_design.R')
out_files_dir <- paste0("../generated/iubf_outs_",run_level,"/")
if(!dir.exists(out_files_dir)) dir.create(out_files_dir)
abf_files_dir <- paste0("../generated/abf_outs_",run_level,"/")
if(!dir.exists(abf_files_dir)) dir.create(abf_files_dir)

num_jobs <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_COUNT'))
job_num <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(run_level==1){
   num_jobs = 1
   job_num = 1
}

for(i in 1:nrow(pd)){
   if((i%%num_jobs)==(job_num%%num_jobs)){
   stew(file=paste0(out_files_dir,"iubf_out_",i,".rda"),seed=5981724,{
         start_params <- unlist(pd[i,])
         print(start_params)
         iubf(object = m,
              Nubf = iubf_nubf,
              Nrep_per_param = iubf_nrep_per_param,
              Nparam = iubf_nparam,
              nbhd = iubf_nbhd,
              params=start_params,
              prop = iubf_prop,
              rw.sd = rw.sd(
                R0=0.02,sigma=0.02,gamma=0.02,sigmaSE=0.02,psi=0.02,amplitude=0.02
              ),
              cooling.type = "geometric",
              cooling.fraction.50 = 0.5,
              verbose=TRUE
         ) -> iubf_out
         final_params <- coef(iubf_out)
         output <- c(U = U,
                     N = N,
                     loglik = iubf_out@loglik,
                     Nabf = iubf_nubf,
                     Nparam = iubf_nparam,
                     Nrep_per_param = iubf_nrep_per_param,
                     prop = iubf_prop,
                     coef(iubf_out)[c("alpha","iota","R0","amplitude","cohort","gamma","sigma","mu","sigmaSE","rho","psi","g")],
                     notes = paste("sim_study",run_level,sep = "_"))
         write.table(x=as.data.frame(t(as.matrix(output))),
                     file="../generated/iubf_out.csv",
                     sep=",",
                     append = TRUE,
                     col.names = FALSE,
                     row.names = FALSE,
                     quote = FALSE)
      })
   }
}

for(i in 1:nrow(pd)){
   if((i%%num_jobs)==(job_num%%num_jobs)){
      stew(file=paste0(abf_files_dir,"abf_loglik_",i,".rda"),seed=5981724,{
         load(paste0(out_files_dir,"iubf_out_",i,".rda"))
         final_params <- coef(iubf_out)
         print(final_params)
         replicate(abf_nruns,
                     bpfilter(object = m,
                              Np = 50000,
                              # Np = 1,
                              # nbhd = iubf_nbhd,
                              params=final_params,
                              block_size=1,
                              verbose=TRUE
                     )) -> abf_runs
         all_ll <- sapply(abf_runs, "slot", "loglik")
         ll <- logmeanexp(all_ll, se = TRUE)
         
         output <- c(U = U,
                     N = N,
                     run_level = run_level,
                     Nabf_runs = abf_nruns,
                     Nrep = 50000,
                     coef(iubf_out)[c("alpha","iota","R0","amplitude","cohort","gamma","sigma","mu","sigmaSE","rho","psi","g")],
                     abf_loglik = ll[1],
                     abf_loglik_se = ll[2],
                     notes = paste("sim_study",run_level,sep = "_"))
         
         write.table(x=as.data.frame(t(as.matrix(output))),
                     file="../generated/abf_out.csv",
                     sep=",",
                     append = TRUE,
                     col.names = FALSE,
                     row.names = FALSE,
                     quote = FALSE)
      })
   }
}
# library(dplyr)
# abf_out <- read.csv('./generated/abf_out.csv', row.names = NULL) %>% dplyr::filter(Nrep==50000)
# plot.profile <- function(mcap1,ylab,xline=2.5,yline=2.5,xlab="",quadratic=FALSE,...){
#    if(missing(ylab)) ylab <- "profile log likelihood"
#    ggplot() + geom_point(aes(x = mcap1$parameter, y = mcap1$lp)) +
#       geom_line(mapping = aes(x = mcap1$fit$parameter, y = mcap1$fit$smoothed),
#                 color = "red") +
#       {if(quadratic) geom_line(mapping = aes(x = mcap1$fit$parameter, y = mcap1$fit$quadratic),
#                                color = "blue",
#                                lwd = 1.25)} +
#       labs(x = xlab, y = ylab) +
#       geom_vline(xintercept = mcap1$ci, color = "red") +
#       geom_hline(yintercept = max(mcap1$fit$smoothed, na.rm = T) - mcap1$delta, color = "red") +
#       theme(panel.border = element_rect(colour = "black", fill=NA))
# }
# 
# g_prof_mcap <- spatPomp::mcap(
#    lp=abf_out[,"abf_loglik"],
#    parameter=abf_out[,"g"],
#    lambda = 0.8
# )
# plot.profile(g_prof_mcap,xlab='g', quadratic = FALSE) +
#    geom_vline(xintercept=400,color="blue",linetype="dashed") +
#    theme(axis.text.x = element_text(size = 12),
#          axis.text.y = element_text(size = 12),
#          axis.title.x = element_text(size = 14),
#          axis.title.y = element_text(size = 14))

   
   
   
   