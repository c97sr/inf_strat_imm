function [Rt_list beta_list] = calR0_byT(theta_beta,par,mean_posterior)
%calR0_byOutput Summary of this function goes here
%theta_beta: average beta (redundant, think whether I need it or not)
%par: from model output
%mean_posteior: average value of posteriorsamples
%This function won't produce credible interval
%1 Feb 2015
% Written by Sean Yuan (hyuan@imperial.ac.uk)
    %p = path;
    %path(p,'lib/');

    beta_size = length(theta_beta);
    beta_list = zeros(1,beta_size);
    R0_list = zeros(1,beta_size);
    
    [par] = setParameters(par,'beta',mean_posterior(1));
    [par] = setParameters(par,'AbB1',mean_posterior(2));
    [par] = setParameters(par,'AbB2',mean_posterior(3));
    [par] = setParameters(par,'AbB3',mean_posterior(4));
    [par] = setParameters(par,'AbB4',mean_posterior(5));
    [par] = setParameters(par,'immune_alpha',mean_posterior(6));
    [par] = setParameters(par,'immune_beta',mean_posterior(7));
    NGM = cal_NGM_bybetaByT(theta_beta,par.Antibody.K,par);      
   
    for t=1:360
        t
        %beta_list(1,i) = theta_beta(i);
        %[A] = cal_NGM_bybeta2(theta_beta(i),par.Antibody.K,par);   %cal_NGM_bybeta2 involves ODE calculation 
        A = NGM(t).A;
        [v d] = eig(A);
        Rt = max(max(d)); %R0 = 1.1181
        Rt_list(1,t) = Rt;
    end
end

