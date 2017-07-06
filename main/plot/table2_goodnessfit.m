function [p] = table2_goodnessfit(u, v, Nu, Nv)
% Summary of the function table2_goodnessfit
% Compare the consistency of histograms u and v
% Written by Sean Yuan (hyuan@imperial.ac.uk) 
% Plot all age group information.
u = Nu*u;
v = Nv*v;
bin_k = length(u);
delta = u-v;
sigmasq = u+v;
Chi = 0;
for i = 1:bin_k
  Chi = Chi + (delta(i).^2)./sigmasq(i);
end

disp(Chi);
p = (1 - chi2cdf(Chi,10));

end