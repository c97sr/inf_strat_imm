function [ g ] = ztpoisspdf( k, lambda )
%ZTPOISS Summary of this function goes here
%Zero-truncated Poisson distribution
%http://en.wikipedia.org/wiki/Zero-truncated_Poisson_distribution
%   Detailed explanation goes here
  g = ((lambda.^k)*(exp(-lambda)))./(factorial(k)*(1-(exp(-lambda))));
  indc = find(k==0);
  if ~isempty(indc)
    g(indc) = 0;
  end
end

