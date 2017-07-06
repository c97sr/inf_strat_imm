function [ f ] = make_f_simple( pa )
%make_f_simple Summary of this function goes here
%   Simple function to return a null infectivity array of the correct dimensons
%   Written by Sean Yuan (hyuan@imperial.ac.uk) 
maxa = pa.maxa; %number of age groups
maxi = pa.maxi; %Ab level to strain A
maxj = pa.maxj; %Ab level to strain B
maxk = pa.maxk; %Ab level to strain C
f = zeros(maxa,maxi,maxj,maxk);
for a=1:maxa
    for i=1:maxi
        for j=1:maxj
            for k=1:maxk
                f(a,i,j,k)=1;
            end    
        end
    end
end

end

