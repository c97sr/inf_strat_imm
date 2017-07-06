function [ ydot ] = odef_islmodjava( t, y, meser )
%odef_islmod Summary of this function goes here
% ODE model to obain disease dynamics
% Written by Sean Yuan (hyuan@imperial.ac.uk) 

% for full titre model, y:160x1 

ydot = zeros(1,length(y));
meser.setStates(y);
%ydot = meser.getStatesDeriv();
ydot = meser.getStatesSIRDeriv();
return;
end

