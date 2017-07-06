function [FigH Rt_list Rt_lower Rt_upper] = figure5(PosteriorSamples, par, burnIn, samplesize)
%Summary of this function goes here
%Calculate Rt and plot it
%Adapted from calR0_byOutput 
%par: from model output
% Written by Sean Yuan (hyuan@imperial.ac.uk)


lastsamplingday = par.SamplingLastDay; %should be 365d
posterior = table2array(PosteriorSamples);
%retrieve parameters from posterior
if exist('samplesize') == 0
        samplesize = 10;
    end
    if exist('burnIn') == 0
        burnIn = 1000;
    end
    total = length(posterior(:,1))-burnIn;
    idx = burnIn + round(rand(1, samplesize) * total);
    
    %estimate of posterior mean
    mean_posterior = mean(posterior(idx,1:end-1));
    
for i = 1:samplesize
    vars = PosteriorSamples.Properties.VariableNames;
    for p=1:length(vars)
        if strcmpi('LLH',vars(p))
        else
           [par] = setParameters(par,char(vars(p)),posterior(idx(i),p));
        end
    end
  
    NGM = cal_NGM_bybetaByT(par.Antibody.K,par);   
    for t=1:360
        A = NGM(t).A;
        [v d] = eig(A);
        Rt = max(max(d)); %R0 = 1.1181
        Rt_list(i,t) = Rt;
    end
    Rt_mean = mean(Rt_list);
end


%Add quantile of Rt
Rt_lower = quantile(Rt_list,0.025);
Rt_upper = quantile(Rt_list,0.975);

FigH = figure;
plot(Rt_mean);
hold on;
plot(Rt_lower);
plot(Rt_upper);

xlim([1 305]);
set(gca,'XTick',[1:30.4:305]);
%set('YTickLabel',);

Xl = [1 lastsamplingday];
months = [
          'May';
          'Jun';
          'Jul';
          'Aug';
          'Sep';
          'Oct';
          'Nov';
          'Dec';
          'Jan';
          'Feb';
          '   '
        ];
     
ax = axis;     % Current axis limits
axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
Yl = ax(3:4);  % Y-axis limits

%Something is wrong here. Need to change the shift days. 20150627
%%%% MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
%%%%
%%%%
% Place the text labels
Xt = [1:(365-60-1)/10:365-60]+15;
t = text(Xt,Yl(1)*ones(1,length(Xt)),months(1:1:11,:));
set(t,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45);

% Remove the default labels
set(gca,'XTickLabel','')
end

