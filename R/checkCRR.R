library(coda)
P1<-as.mcmc(posterior)

P2<-P1[,1:9];
colnames(P2) <- c("beta","AbB1","AbB2","AbB3","AbB4","PT1","PT2","PT3","PT4")

idx <- sample(1:10000, 400, replace = FALSE, prob = NULL)
pairs(P2[10000:10000+idx,])