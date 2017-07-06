function [Rt0] = calculateR0_fromPars(par)
%Summary of this function goes here
%Calculate Rt and plot it
%Adapted from calR0_byOutput 
%par: from model output
%Without sampling. Calculate using all data.
% Written by Sean Yuan (hyuan@imperial.ac.uk)

    Ab = par.Antibody;
    
    [NGM] = cal_NGM_T0(Ab,par);   
    A = NGM.A;
    [v d] = eig(A);
    Rt0 = max(max(d)); %R0 = 1.1181
    
    %NGM = cal_NGM_bybetaByT(par.Antibody.K,par);   
    %A = NGM.A;
    %[v d] = eig(A);
    %Rt0 = max(max(d)); %R0 = 1.1181
end

