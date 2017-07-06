function [ g cnt] = make_g( pa )
%make_g_simple Summary of this function goes here
% Function to return the most simple immune state boost
% Written by Sean Yuan (hyuan@imperial.ac.uk) 
maxX = pa.maxX; %number of strains
maxa = pa.maxa; %number of age groups
maxl = pa.maxi; %Ab level to strain A
maxm = pa.maxj; %Ab level to strain B
maxn = pa.maxk; %Ab level to strain C
g = zeros(maxX,maxa,maxl,maxm,maxn,maxl,maxm,maxn);

pa.AbB=3;% temp
if pa.age_flag == 0
  pa.AbB1 = pa.AbB;
  pa.AbB2 = pa.AbB;
  pa.AbB3 = pa.AbB;
  pa.AbB4 = pa.AbB; 
end
%----------------- Based on Poisson --------------------%
for i=1:pa.maxi
    boostpr = ztpoisspdf(0:pa.maxi-i,pa.AbB);
    boostpr1 = ztpoisspdf(0:pa.maxi-i,pa.AbB1);
    boostpr2 = ztpoisspdf(0:pa.maxi-i,pa.AbB2);
    boostpr3 = ztpoisspdf(0:pa.maxi-i,pa.AbB3);
    boostpr4 = ztpoisspdf(0:pa.maxi-i,pa.AbB4);
    boostpr(end) = 1-sum(boostpr(1:end-1));
    boostpr1(end) = 1-sum(boostpr1(1:end-1));
    boostpr2(end) = 1-sum(boostpr2(1:end-1));
    boostpr3(end) = 1-sum(boostpr3(1:end-1));
    boostpr4(end) = 1-sum(boostpr4(1:end-1));
    Boosting(i).abl = boostpr;
    age(1).Boosting(i).abl = boostpr1;
    age(2).Boosting(i).abl = boostpr2;
    age(3).Boosting(i).abl = boostpr3;
    age(4).Boosting(i).abl = boostpr4;
end
cnt = 0;

for X=1:maxX
    for a=1:maxa
        for l=1:maxl
            for m=1:maxm
                for n=1:maxn
                    if X==1
                            g(X,a,l,m,n,l+[0:pa.maxi-l],m,n)=age(a).Boosting(l).abl;
                    end
                    if X==2
                        if(m<=maxm-abl)
                            g(X,a,l,m,n,l,m+abl,n)=1;
                            cnt = cnt + 1;
                        else
                            g(X,a,l,m,n,l,m,n)=1;
                            cnt = cnt + 1;
                        end
                    end
                    if X==3
                        if(n<=maxn-abl)
                            g(X,a,l,m,n,l,m,n+abl)=1;
                            cnt = cnt + 1;
                        else
                            g(X,a,l,m,n,l,m,n)=1;
                            cnt = cnt + 1;
                        end
                    end
                end
            end
        end
    end %--end of age loop
end
end
