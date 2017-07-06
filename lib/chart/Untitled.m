function [ infecteds] = gen_total_infecteds( y, times, pars, a)
%gen_total_infecteds Summary of this function goes here
%   Detailed explanation goes here
% y: model output
% times: sampling time
% 

    dims_Slu = prod(size(pars.arrSlu)); %[a i j x]
    dims_Ilu = prod(size(pars.arrIlu)); %[X a i j k]
    infecteds = sum(y(:,dims_Slu+1:dims_Slu+dims_Ilu),2); % total infected number, regardless which strains
    % add the number of infecteds by strain

end

