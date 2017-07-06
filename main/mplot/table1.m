function [S1 T] = table1(PosteriorSamples, burnIn)
% Summary of the function table1
% Calculate the mean and credible intervals for parameters:
% [Beta AbB1 AbB2 AbB3 AbB4 PT1 PT2 PT3 PT4]
% Written by Sean Yuan (hyuan@imperial.ac.uk) 
% Plot all age group information.

%read posterior
posterior_all = table2array(PosteriorSamples);
posterior = posterior_all(burnIn:end,:);
%number of posterior variables
%var_name
%posterior;
post_bar = mean(posterior);
%quantile(posterior, [0.05 0.95])
post_range = quantile(posterior, [0.025 0.975]);
S1 = '';
for i=1:length(post_bar)
    str = [char(PosteriorSamples.Properties.VariableNames(i)) ',' num2str(post_bar(i)) '[' num2str(post_range(1,i)) '-' num2str(post_range(2,i)) ']']; 
    S1 = strvcat(S1, str); 
end

post_conf = {};
for i=1:length(post_bar)-1
    post_conf(i) = {['[' num2str(post_range(1,i)) '-' num2str(post_range(2,i)) ']']}; 
end

RowName = PosteriorSamples.Properties.VariableNames(1:end-1);
PostBar = post_bar(1:end-1);

%Shape = {'Pan';'Round';'Button';'Pan';'Round'};
%Price = [10.0;13.59;10.50;12.00;16.69];
%Stock = [376;502;465;1091;562];
T = table(RowName',PostBar',post_conf','RowNames',RowName)
T.Properties.VariableNames = {'Parameters' 'Mean' 'CI'};



end