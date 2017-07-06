getess <-function(filename) {
#library('coda')
#setwd("E:/working/Projects.IC/Projects/isl/mat/Misl/isl-3.3/out/mcmc/ph1n1/20141028");
ps <- read.csv('../../out/mcmc/ph1n1/20141029/mcmc_output_m6_posterior.csv');
ess <- effectiveSize(ps);
return(ess);
}