function [ pd ] = truncate_normpdf( x, mu, sd, lb, ub )
%truncate_normpdf Summary of this function goes here
%   Detailed explanation goes here
     pd = normpdf(x,mu,sd)/(normcdf(ub,mu,sd) - normcdf(lb,mu,sd));
     if x < lb
         pd = 0;
     end
     if x > ub
         pd = 0;
     end
end

