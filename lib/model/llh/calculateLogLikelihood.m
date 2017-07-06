function [ LLH ] = calculateLogLikelihood( multi_p, obs )
%calculateLogLikelihood Summary of this function goes here
%estimate loglikelihood of observed data giving multinomial probability
%Detailed explanation goes here

    LLH = 0;
    for i = 1:length(obs)
         LLH = LLH + log(multi_p(obs(i)));
    end
    %log_multi_factor = calculateConstant(multi_p, obs);
    log_multi_factor = 0;
    LLH = log_multi_factor+LLH;
end


