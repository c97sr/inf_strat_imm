function [logmf] = calculateConstant(multi_p, obs)
%calculateConstant Summary of this function goes here
%calculate the constant factor for multinomial distribution likelihood
%called by calculateLogLikelihood()
%   Detailed explanation goes here
% 24 Oct 2014

    instances = histc(obs,[1:length(multi_p)]);
    denom = log(1);
    for i=1:length(multi_p)
        denom = denom + sum(log(1:instances(i)));
    end
    logmf = sum(log(1:(length(obs)))) - denom;

end

