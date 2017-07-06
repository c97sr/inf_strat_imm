function [Rt0_mean Rt0_lower Rt0_upper] = calculateR0_2(PosteriorSamples, par)
%Summary of this function goes here
%Calculate Rt and plot it
%With sampling
% Written by Sean Yuan (hyuan@imperial.ac.uk)


lastsamplingday = par.SamplingLastDay; %should be 365d
posterior = table2array(PosteriorSamples);
%retrieve parameters from posterior

    samplesize = 20;
    burnIn = 1000;
    total = length(posterior(:,1))-burnIn;
    idx = burnIn + round(rand(1, samplesize) * total);
    
    
for i = 1:length(idx)
    vars = PosteriorSamples.Properties.VariableNames;
    for p=1:length(vars)
        if strcmpi('LLH',vars(p))
        else
           [par] = setParameters(par,char(vars(p)),posterior(idx(i),p));
        end
    end
    Ab = par.Antibody;
    [NGM] = cal_NGM_T0(Ab,par);   
    A = NGM.A;
    [v d] = eig(A);
    Rt0 = max(max(d)) %R0 = 1.1181
    Rt0_list(i) = Rt0;
    if rem(i,10) == 0
      disp(i);
    end
end
Rt0_mean = mean(Rt0_list);
%Add quantile of Rt
Rt0_lower = quantile(Rt0_list,0.025);
Rt0_upper = quantile(Rt0_list,0.975);


end

