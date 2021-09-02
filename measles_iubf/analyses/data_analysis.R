library("spatPomp")
library("doRNG")
library("doParallel")

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
 run_level <- 1
} else run_level <- as.numeric(args)

source('../analyses/iubf_measles_sim.R')
set.seed(1)
U<-10;N<-2*26
m <- measles_sim_w_params(U=U,N=N)
md <- as.data.frame(m)

library(ggplot2)
ggplot(data = md, mapping = aes(x = year, y = cases)) + 
 geom_line() +
 facet_wrap(~city, scales = "free_y")

m2 <- measles_sim_w_params(U=U,N=N,mup=c(alpha=1,iota=0,R0=41.23672,cohort=0,amplitude=0.5354223, gamma=45.61916, sigma=50.27219,mu=0.02,sigmaSE=0.1351309, rho=0.5,psi=0.2071051,g=1100))
md2 <- as.data.frame(m2)
ggplot(data = md2, mapping = aes(x = year, y = cases)) + 
 geom_line() +
 facet_wrap(~city, scales = "free_y")

m3 <- measles_sim_w_params(U=U,N=N,mup=c(alpha=1,iota=0,R0=46.77458,cohort=0,amplitude=0.6040571, gamma=44.81039, sigma=37.18639,mu=0.02,sigmaSE=0.1439365, rho=0.5,psi=0.2084648,g=455.9014))
md3 <- as.data.frame(m3)
ggplot(data = md3, mapping = aes(x = year, y = cases)) + 
 geom_line() +
 facet_wrap(~city, scales = "free_y")
p1 <- coef(m)
p2 <- coef(m)
p2[['g']] <- 10000
# p1[c("g","R0", "amplitude", "gamma", "sigma", "sigmaSE", "psi")] <- c(1100,
#                                                                      41.23672,
#                                                                      0.5354223,
#                                                                      45.61916,
#                                                                      50.27219,
#                                                                      0.1351309,
#                                                                      0.2071051)
# p2[c("g","R0", "amplitude", "gamma", "sigma", "sigmaSE", "psi")] <- c(455.9014,
#                                                                       46.77458,
#                                                                       0.6040571,
#                                                                       44.81039,
#                                                                       37.18639,
#                                                                       0.1439365,
#                                                                       0.2084648)
md <- as.data.frame(m) %>% dplyr::select(year, city, cases, pop, lag_birthrate)
d1 <- simulate(m, nsim = 10, params = p1, format="data.frame")
d2 <- simulate(m, nsim = 10, params = p2, format="data.frame")
d1 <- d1 %>% dplyr::rename(sim_cases = cases,
                           city = unitname)
d2 <- d2 %>% dplyr::rename(sim_cases = cases,
                           city = unitname)
d1_final <- d1 %>% 
 dplyr::inner_join(md,by = c('year','city')) %>%
 dplyr::select(year,.id,city,cases,sim_cases,C,S,E,I,R)
d2_final <- d2 %>% 
 dplyr::inner_join(md,by = c('year','city')) %>%
 dplyr::select(year,.id,city,cases,sim_cases,C,S,E,I,R)

d1 <- d1_final %>% dplyr::mutate(first_term = pnorm(cases+0.5, 
                                              mean = p1[['rho']]*C,
                                              sd = sqrt(p1[['rho']]*C*(1-p1[['rho']]+p1[['psi']]*p1[['psi']]*p1[['rho']]*C)+1)))
d1 <- d1 %>% dplyr::mutate(second_term = pnorm(cases-0.5, 
                                              mean = p1[['rho']]*C,
                                              sd = sqrt(p1[['rho']]*C*(1-p1[['rho']]+p1[['psi']]*p1[['psi']]*p1[['rho']]*C)+1)))
d1 <- d1 %>% dplyr::mutate(meas_density = first_term - second_term)


d2 <- d2_final %>% dplyr::mutate(first_term = pnorm(cases+0.5, 
                                              mean = p2[['rho']]*C,
                                              sd = sqrt(p2[['rho']]*C*(1-p2[['rho']]+p2[['psi']]*p2[['psi']]*p2[['rho']]*C)+1)))
d2 <- d2 %>% dplyr::mutate(second_term = pnorm(cases-0.5, 
                                               mean = p2[['rho']]*C,
                                               sd = sqrt(p2[['rho']]*C*(1-p2[['rho']]+p2[['psi']]*p2[['psi']]*p2[['rho']]*C)+1)))
d2 <- d2 %>% dplyr::mutate(meas_density = first_term - second_term)
d1$meas_density2 <- d2$meas_density
meas_density_diff_by_city <- d1 %>% 
 dplyr::filter(.id==1) %>%
 dplyr::group_by(city) %>%
 dplyr::summarize(prop_true_better = sum(meas_density > meas_density2)/dplyr::n(),
                  prop_true_worse = sum(meas_density < meas_density2)/dplyr::n(),
                  prop_equal = sum(meas_density == meas_density2)/dplyr::n())
meas_density_diff_by_id <- d1 %>% 
 dplyr::group_by(.id) %>%
 dplyr::summarize(prop_true_better = sum(meas_density > meas_density2)/dplyr::n(),
                  prop_true_worse = sum(meas_density < meas_density2)/dplyr::n(),
                  prop_equal = sum(meas_density == meas_density2)/dplyr::n())
bpfilter(m, params = p1, Np = 10000, block_size = 2)@loglik
bpfilter(m, params = p2, Np = 10000, block_size = 2)@loglik

