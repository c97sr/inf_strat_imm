function [ h ] = make_h( pa )
%make_h Summary of this function goes here
% Takes a scaler set of parameters and returns an array with values for h.
% h(X,a,i,j,k) is the susceptibility of individuals to strain X dependent on their
% age and their antibody levels. Values are referenced to h(X,1,1,1,1)=1, i.e.
% those in the youngest age group with no detectable titres always have a relative
% suceptibility of 1
%   Written by Sean Yuan (hyuan@imperial.ac.uk) 
maxX = pa.maxX; %number of strains
maxa = pa.maxa; %number of age groups
maxi = pa.maxi; %Ab level to strain A
maxj = pa.maxj; %Ab level to strain B
maxk = pa.maxk; %Ab level to strain C
h = zeros(maxX,maxa,maxi,maxj,maxk);
immune_alpha = pa.immune_alpha; %larger immune_alpha more susceptible
immune_beta = pa.immune_beta;
%if pa.immune_flag == 1
  imm_alpha(1) = pa.immune_alpha1;
  imm_alpha(2) = pa.immune_alpha2;
  imm_alpha(3) = pa.immune_alpha3;
  imm_alpha(4) = pa.immune_alpha4;
%elseif pa.immune_flag == 0
%  imm_alpha(1) = pa.immune_alpha;
%  imm_alpha(2) = pa.immune_alpha;
%  imm_alpha(3) = pa.immune_alpha;
%  imm_alpha(4) = pa.immune_alpha; 
%end

for X=1:maxX
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                    if X==1
                        h(X,a,i,j,k) = getSusc(i-1, imm_alpha(a), immune_beta);     %i-1 = titres from 0:8         
                    end
                    if X==2
                        h(X,a,i,j,k)=0;
                    end
                    if X==3
                        h(X,a,i,j,k)=0;
                    end
                end
            end
        end
    end
end
end
