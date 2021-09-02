library("spatPomp")
library("ggplot2")
read.csv('../generated/iubf_out.csv', row.names = NULL) -> g_prof
g_prof <- g_prof %>% dplyr::filter(U==20)
plot.profile <- function(mcap1,ylab,xline=2.5,yline=2.5,xlab="",quadratic=FALSE,...){
 if(missing(ylab)) ylab <- "profile log likelihood"
 ggplot() + geom_point(aes(x = mcap1$parameter, y = mcap1$lp)) +
  geom_line(mapping = aes(x = mcap1$fit$parameter, y = mcap1$fit$smoothed),
            color = "red") +
  {if(quadratic) geom_line(mapping = aes(x = mcap1$fit$parameter, y = mcap1$fit$quadratic),
                           color = "blue",
                           lwd = 1.25)} +
  labs(x = xlab, y = ylab) +
  geom_vline(xintercept = mcap1$ci, color = "red") +
  geom_hline(yintercept = max(mcap1$fit$smoothed, na.rm = T) - mcap1$delta, color = "red") +
  theme(panel.border = element_rect(colour = "black", fill=NA))
}

g_prof_mcap <- spatPomp::mcap(
 lp=g_prof[,"loglik"],
 parameter=g_prof[,"g"],
 lambda = 0.8
)
plot.profile(g_prof_mcap,xlab='g', quadratic = FALSE) +
 geom_vline(xintercept=400,color="blue",linetype="dashed") +
 theme(axis.text.x = element_text(size = 12),
       axis.text.y = element_text(size = 12),
       axis.title.x = element_text(size = 14),
       axis.title.y = element_text(size = 14))

# ggplot(data = g_prof, mapping = aes(x = g, y = loglik))+
#  geom_point() + 
#   geom_vline(xintercept=400,color="blue",linetype="dashed")
#   

# measles_unit_params <- c(
#   alpha=1,
#   iota=0,  # set to zero for a closed population
#   R0=30,
#   cohort=0,
#   amplitude=0.5, gamma=52, sigma=52,mu=0.02,
#   sigmaSE=0.15, rho=0.5,
#   psi=0.15,
#   g=400
# )

# p1<-ggplot(data = g_prof, mapping = aes(x = g, y = gamma, color = loglik))+
#   geom_point()
# p2<-ggplot(data = g_prof, mapping = aes(x = g, y = sigma, color = loglik))+
#   geom_point()
# p3<-ggplot(data = g_prof, mapping = aes(x = g, y = R0, color = loglik))+
#   geom_point()
# p4<-ggplot(data = g_prof, mapping = aes(x = g, y = amplitude, color = loglik))+
#   geom_point()
# p5<-ggplot(data = g_prof, mapping = aes(x = g, y = sigmaSE, color = loglik))+
#   geom_point()
# p6<-ggplot(data = g_prof, mapping = aes(x = g, y = psi, color = loglik))+
#   geom_point()
# library(grid)
# grid.newpage()
# grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), ggplotGrob(p4), ggplotGrob(p5), ggplotGrob(p6)))
