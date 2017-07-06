function [R0_list beta_list] = calR0_bypar(theta_beta,par)
%calR0_bypar Summary of this function goes here
%Use default parameter set
% Written by Sean Yuan (hyuan@imperial.ac.uk)
    %p = path;
    %path(p,'lib/');

    beta_size = length(theta_beta);
    beta_list = zeros(1,beta_size);
    R0_list = zeros(1,beta_size);
    
    
    for i=1:beta_size
        beta_list(1,i) = theta_beta(i);
        %[A] = cal_NGM_bybeta2(theta_beta(i),par.Antibody.K,par);   %cal_NGM_bybeta2 involves ODE calculation
        [A] = cal_NGM_bybeta_naive(theta_beta(i),par); 
        [v d] = eig(A);
        R0 = max(max(d)); %R0 = 1.1181
        R0_list(1,i) = R0;
    end
end

