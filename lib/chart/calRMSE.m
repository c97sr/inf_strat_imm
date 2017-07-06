function [RMSE_part RMSE_full] = calRMSE(rep, Xt, Yt, N)
% calculate P(Y|X,rep)
% Use Poisson or binomial
%binopdf(Yt,XtN,P)  
  %P = [ones(1,7) geopdf([0:4],0.5) ones(1,32)*0.05];
  P = [ones(1,7) (1+0.052)./2 ones(1,36)*0.052];
  P = P'*rep;
  % first 7 weeks 100%
  % 8th week 50%
  % after 9th week 5.2%
  
  lambda = Xt*N.*P;
  diff = abs(lambda - Yt);
  diff_target = diff(1:23);
  RMSE_part = (diff_target'*diff_target)^0.5;
  
  diff_full = diff(1:end);
  RMSE_full = (diff_full'*diff_full)^0.5;
  
  
end