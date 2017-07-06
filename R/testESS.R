testESS <- function() {
x <- 1:3
print(x)
library('coda')
ps <- read.csv('temp/posterior.csv')
ess <- effectiveSize(ps);
write.csv(ess, file = "temp/MyData1.csv")
#sink("temp/output_ess.txt")
#print(ess)
#sink()
#fileConn<-file("output_ess.txt")
#writeLines(c(ess), fileConn)
#close(fileConn)
}