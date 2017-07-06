initParameters <- function() {
  # install.packages("rJava")
  #https://darrenjw.wordpress.com/2011/01/01/calling-java-code-from-r/
  #http://www.codophile.com/how-to-integrate-r-with-java-using-rjava/
  Ps <- c(
    maxX = 1,
    maxa = 1,
    maxi = 10,
    maxj = 1,
    maxk = 1
    #stringsAsFactors = FALSE
    )
  
    #matG <- make.G(Ps)
}
make.G.test<-function(ps) {
  print(ps.maxX)
  
}
#X number of strains
#l,m,n are the initial states
#i,j,k are the final states
#G(X,a,l,m,n,i,j,k)  

#maxX = ps.maxX; %number of strains
#maxa = ps.maxa; %number of age groups
#maxl = ps.maxi; %Ab level to strain A
#maxm = ps.maxj; %Ab level to strain B
#maxn = ps.maxk; %Ab level to strain C

#ps.maxX=1
#ps.maxa=1
#ps.maxi=10
#ps.maxj=1
#ps.maxk=1

#for (i in 1:ps.maxi-i) {
#    boosting <- dpois(1:ps.maxi-i, lambda = 3);
#    boosting[length(boosting)] <- 1-sum(boosting[1:length(boosting)-1])
#}

make.G<-function(ps) {
rtn <- array(0, dim=c (ps["maxX"],ps["maxa"],ps["maxi"],ps["maxj"],ps["maxk"],ps["maxi"],ps["maxj"],ps["maxk"]))
#print(ps["maxX"])
for (X in 1:ps["maxX"]) {
     for (a in 1:ps["maxa"]) {
          for (l in 1:ps["maxi"]) {
               for (m in 1:ps["maxj"]) {
                    for (n in 1:ps["maxk"]) {
                        if (identical(as.numeric(X),1)) { # for strain 1
                            boosting <- dpois(1:(ps["maxi"]-l), lambda = 3)
                            boosting[length(boosting)] <- 1-sum(boosting[1:length(boosting)-1])
                            
                            
                            if (l < ps["maxi"]){ 
                                #print(boosting)
                                rtn[X,a,l,m,n,l+1:(ps["maxi"]-l),m,n] = boosting
                            }
                            #
                            #rtn[X,a,l,m,n,l+[1:pa.maxi],m,n] = boosting
                            #boosting <- 0;
                        }
                    }
               }
          }
     }
}
rtn
}


make.F<-function(ps) {
  for (a in 1:ps["maxa"]) {
    for (i in 1:ps["maxi"]) {
      for (j in 1:ps["maxj"]) {
        for (k in 1:ps["maxk"]) {
          rtn(a,i,j,k) <- 1
        }    
      }
    }
  }
  rtn
}

make.H<-function(ps) {
  rtn <- array(0, dim=c (ps.maxX,ps.maxa,ps.maxi,ps.maxj,ps.maxk));
}

make.M<-function(ps) {}

vv <- function(d) {print(2)}
