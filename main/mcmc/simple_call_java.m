% Simple script to demonstrate how to call java from matlab to calculate the State changes by time
% Practically, it should be used with matlab ODE solver 


y0 = [0.172479753600000;0.000879998742857143;0.000879998742857143;0.000879998742857143;0.000879998742857143;0;0;0;0;0;0.291059584200000;0.00148499787857143;0.00148499787857143;0.00148499787857143;0.00148499787857143;0;0;0;0;0;0.386119448400000;0.00196999718571429;0.00196999718571429;0.00196999718571429;0.00196999718571429;0;0;0;0;0;0.127679817600000;0.00132999810000000;0.00132999810000000;0.00132999810000000;0.00132999810000000;0;0;0;0;0;2.51428571428572e-07;0;0;0;0;0;0;0;0;0;4.24285714285714e-07;0;0;0;0;0;0;0;0;0;5.62857142857143e-07;0;0;0;0;0;0;0;0;0;1.90000000000000e-07;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
%1. initialize java object
javaaddpath 'e:\workspace\MyJavaProject\bin\matlabjava.jar';
%javaaddpath(par.javapath);
import matlabjava.*
mepar = matlabjava.Parameters;
meser = matlabjava.Serology;
meser.setParameters(mepar);

%2. initialize parameters values
par = InitParameters(); 
s0_imm = [0,0,0,0];
wan = 0.04;
maxi = 10;
arrg = par.arrg; %6D  1x4x10x1x1x10x1x1  (X a l m n l+g m n)
arrh = par.arrh; %1x4x10x1x1
matM = par.matM; %4x4
beta = par.beta; %0.5

% pass parameters values to have object
meser.updateParametersG(arrg);
meser.updateParametersH(arrh);
meser.updateParametersM(matM);
meser.updateParametersBeta(beta);
meser.updateParameters('s0_imm',s0_imm);
meser.updateParameters('wan',wan);
meser.updateParameters('maxi',maxi);

% set status
meser.setStates(y0);

% calculate the status changes
ydot = meser.getStatesSIRDeriv();
