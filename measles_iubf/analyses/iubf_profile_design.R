m <- measles_sim_w_params(U=U,N=N)
theta <- coef(m)

theta_start_t <- partrans(m,theta,"toEst")
theta_t_hi <- theta_t_lo <- theta_start_t
theta_t_lo[c("sigma",
             "gamma",
             "sigmaSE",
             "psi",
             "amplitude",
             "R0")] <- theta_t_lo[c("sigma",
                                    "gamma",
                                    "sigmaSE",
                                    "psi",
                                    "amplitude",
                                    "R0")] - log(2)
theta_t_hi[c("sigma",
             "gamma",
             "sigmaSE",
             "psi",
             "amplitude",
             "R0")] <- theta_t_hi[c("sigma",
                                    "gamma",
                                    "sigmaSE",
                                    "psi",
                                    "amplitude",
                                    "R0")] + log(2)

if(!file.exists('../generated/iubf_out.csv')) 
        writeLines(text = paste0(
                                c(
                                        c("U","N","loglik","Nabf","Nparam","Nrep_per_param","prop"),
                                        names(theta_start_t),
                                        "notes"
                                ),
                                collapse = ','),
                   con = '../generated/iubf_out.csv')
if(!file.exists('../generated/abf_out.csv')) 
        writeLines(text = paste0(
                c(
                        c("U","N","run_level","Nabf_runs","Nrep"),
                        names(theta_start_t),
                        "abf_loglik",
                        "abf_loglik_se",
                        "notes"
                ),
                collapse = ','),
                con = '../generated/abf_out.csv')


theta_t_lo <- partrans(m,theta_t_lo,"fromEst")
theta_t_lo <- theta_t_lo[!names(theta_t_lo)=='g']
theta_t_hi <- partrans(m,theta_t_hi,"fromEst")
theta_t_hi <- theta_t_hi[!names(theta_t_hi)=='g']

profile_design(
 g=all_g,
 lower=theta_t_lo,
 upper=theta_t_hi,
 nprof=g_nprof
) -> pd

row_order = c()
for(i in seq_len(g_nprof)){
 to_add <- which(seq_len(nrow(pd)) %% g_nprof == i-1)
 row_order <- c(row_order, to_add)
}

pd <- pd[row_order,]
rownames(pd) <- NULL
