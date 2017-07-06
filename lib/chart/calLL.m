function [total_nll] = calNLL(rep, Xt, Yt, N)
% calculate P(Y|X,rep)
% Use Poisson or binomial
%binopdf(Yt,XtN,P)  
  P = [ones(1,7) geopdf([0:4],0.5) ones(1,32)*0.05];
  P = P'*rep;
  % first 7 weeks 100%
  % 8th week 50%
  % after 9th week 5.2%
  
  lambda = Xt*N.*P;
  total_ll = 0;
  %loglike = log(binopdf(Yt, Xt*N, rep));
  loglike = log(poisspdf(Yt,lambda));
  total_ll = sum(loglike);
  total_nll = -total_ll;
%  for i=1:length(Yt)
%    loglike = log(prob(i)*N*rep)*Yt(i);
%    total_ll = total_ll + loglike;
%  end
%  total_ll = total_ll + loglike;
%  total_nll = -total_ll;
end