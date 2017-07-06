function [ matFOI_byS ] = getMOF_byS( pars,y )
%getMOF_byS Summary of this function goes here
% Return the matrix of Force of Infection by States
% FOI_byS =: beta x M x f
% For Next generation matrix estimates
% Written by Sean Yuan (hyuan@imperial.ac.uk)
% May 21, 2014
maxX = pars.maxX; %number of strains
maxa = pars.maxa; %number of age groups
maxi = pars.maxi; %Ab level to strain A
maxj = pars.maxj; %Ab level to strain B
maxk = pars.maxk; %Ab level to strain C
beta = pars.beta;
arrf = pars.arrf; %infectivity array
arrSlu = pars.arrSlu; %I lookup
arrIlu = pars.arrIlu; %I lookup
matM = pars.matM; 
lambda = 0; %force of infection
matFOI_byS = zeros(maxX, maxa, maxa, maxi, maxj, maxk); %maxFOI(X,a,b,i)
%calculate for of infection
for X=1:maxX
    for a=1:maxa %a represents the age groups in susceptible individuals
        tmpval_transmisibility = 0;
        for b=1:maxa  %b represents the age groups in infected individuals
            %tmpval = 0;
            for i=1:maxi
                for j=1:maxj
                    for k=1:maxk
                        %tmpval = tmpval + arrf(b,i,j,k)*y(arrIlu(X,b,i,j,k));
                        %assuming b infects a
                        %matFOI_byS(X,a,b,i,j,k) = beta*matM(a,b)*arrf(b,i,j,k)*y(arrIlu(X,b,i,j,k));
                        matFOI_byS(X,a,b,i,j,k) = beta*matM(a,b)*arrf(b,i,j,k);
                    end
                end
            end
           %tmpval_transmisibility = tmpval_transmisibility + tmpval*matM(a,b);   %%%<<< Beware to use summation
        end
        %matFOI(X,a) = tmpval_transmisibility*beta; %% something is going wrong
    end
end
return;
end

