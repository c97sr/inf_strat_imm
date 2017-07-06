function [Rt0_mean Rt0_lower Rt0_upper] = calculateR0(PosteriorSamples, par)
%Summary of this function goes here
%Calculate Rt and plot it
%Adapted from calR0_byOutput 
%par: from model output
%Without sampling. Calculate using all data.
% Written by Sean Yuan (hyuan@imperial.ac.uk)


lastsamplingday = par.SamplingLastDay; %should be 365d
posterior = table2array(PosteriorSamples);
%retrieve parameters from posterior

    
for idx = 1:length(posterior(:,1))
    vars = PosteriorSamples.Properties.VariableNames;
    for p=1:length(vars)
        if strcmpi('LLH',vars(p))
        else
           [par] = setParameters(par,char(vars(p)),posterior(idx,p));
        end
    end
    Ab = par.Antibody;
    [NGM] = cal_NGM_T0(Ab,par);   
    A = NGM.A;
    [v d] = eig(A);
    Rt0 = max(max(d)); %R0 = 1.1181
    Rt0_list(idx) = Rt0;
    Rt0_mean = mean(Rt0_list);
    if rem(idx,100) == 0
      disp(idx);
    end
end


%Add quantile of Rt
Rt0_lower = quantile(Rt_list,0.025);
Rt0_upper = quantile(Rt_list,0.975);
save(['Rt0_m' num2str(par.model) '.mat'], 'Rt0_mean','Rt0_lower','Rt0_upper');
end

