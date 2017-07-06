

ps <- read.csv('m2_posterior.csv');
pnorm(abs(geweke.diag(mcmc(ps$X0.057884))$z),lower.tail = FALSE)*2
mc.draw <- mcmc(ps)
summary(mh.draws)
http://math.arizona.edu/~piegorsch/675/GewekeDiagnostics.pdf
http://www.people.fas.harvard.edu/~plam/teaching/methods/convergence/convergence_print.pdf
