function [FigH FigL log peak peak_lb peak_ub] = figure2_poster()
% Summary of figur2
% Plot disease and serological dynamics
% For each age group, 5 subfigures are plotted in a row.
% fig2a, heat map of immune dynamics
% fig2b, disease dynamics
% no resampling but onlty mle
% Written by Sean Yuan (hyuan@imperial.ac.uk) 
% Plot all age group information.
% Return [FigH FigL log peak peak_lb peak_ub]
% also return the likelihood ratio test
% H0 Standard threshold model
% H1 Titre model
% Actual cumulative incidence is defined as the average of CI(titre)&CI(threshold)

p = path;
%path(p,'../');
path(p,'lib/');

global proj Antibody;
init_collect = 1;
second_collect = 2;
third_collect = 3;
k = 2;
Ab = Antibody;

%retrieve parameters from posterior
%model:3
if exist('samplesize') == 0
    samplesize = 10;
end
if exist('burnIn') == 0
    burnIn = 1000;
end

%MLE approach
%idx = find(posterior(:,llhidx)==max(posterior(:,llhidx)));
%samplesize = length(idx);






%Mean parameters
FigL = figure;
set(FigL, 'Position', [100, 500, 816, 480]);
set(gca,'xtick',[]);
hold;
dat1 = load('out/p0e05/m1/ph1n1/20151024/mcmc_output_m1_final.mat');
%dat1 = load('out/p0e05/m1.5/ph1n1/20151027/mcmc_output_m1.5_final.mat');
plotI(dat1, 1);

dat2 = load('out/p0e05/m2/ph1n1/20151024/mcmc_output_m2_final.mat');
%dat2 = load('out/p0e05/m2.5/ph1n1/20160226/mcmc_output_m2.5_final.mat');
mark = ':';
plotI(dat2, 1, mark);

function [] = plotI(dat, display, marker)
    
posteriorTable = dat.PosteriorSamples;
pars = dat.par;
post = table2array(posteriorTable);
posterior = repmat(mean(post),3,1);
samplesize = 3;
idx = 1:3;



for i = 1:samplesize
    vars = posteriorTable.Properties.VariableNames;
    for p=1:length(vars)
        if strcmpi('LLH',vars(p))
        else
           [pars] = setParameters(pars,char(vars(p)),posterior(idx(i),p));
        end
    end

    %set parameters
    beta = pars.beta;
    AbB = [pars.AbB1 pars.AbB2 pars.AbB3 pars.AbB4];
    immune_alpha = [pars.immune_alpha1 pars.immune_alpha2 pars.immune_alpha3 pars.immune_alpha4];
    lastsamplingday = pars.SamplingLastDay + 90;% -60 -> +30
    if sum(pars.arrh(1,1:4,2))<0.0000001
      pars.arrh(1,1:4,2);
    end
    
    %setup initial condition
    ab_baseline = Ab.K(init_collect).Abl;
    ab_k = Ab.K(k).Abl;
    if pars.maxi == 2 % only 2 titres
        [yini age_arr s0_imm] = make_ics_naive2titres( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu);
        ab_baseline(find(ab_baseline <3)) = 0;
        ab_baseline(find(ab_baseline >2)) = 1;
        ab_k(find(ab_k <3)) = 0;
        ab_k(find(ab_k >2)) = 1;
    else
        [yini age_arr] = make_ics_naive( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu);
    end
    [yini_k1 age_arr_k1] = make_ics_fromtitres_byage( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, ab_baseline, Ab.K(init_collect).age);
    [yini_k2 age_arr_k2] = make_ics_fromtitres_byage( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, ab_k, Ab.K(k).age);

    %setep simulation time
    T0 = pars.OutbreakStartingDay;
    meanKdays(1) = mean(pars.Antibody.K(1).numdays - T0);
    meanKdays(2) = mean(pars.Antibody.K(2).numdays - T0);
    sample_time_K1 = round(meanKdays(1));
    sample_time_K2 = round(meanKdays(2));
    times = 0:1:lastsamplingday;
    sample_size_K1 = Ab.K(1).samplesize;
    sample_size_K2 = Ab.K(2).samplesize;

    %run simulation
    %initialize objects
    %javaaddpath e:\workspace\MyJavaProject\bin\matlabjava.jar; %use jave working dir 
    javaaddpath(pars.javapath); %set ./java as default dir
    import matlabjava.*
    mepar_2 = matlabjava.Parameters;
    meser_2 = matlabjava.Serology;
    % set parameters
    meser_2.setParameters(mepar_2);
    meser_2.updateParameters('s0_imm',pars.s0_imm);
    meser_2.updateParameters('wan',pars.wan);
    meser_2.updateParameters('maxi', pars.maxi);
    meser_2.updateParametersG(pars.arrg);
    meser_2.updateParametersH(pars.arrh);
    meser_2.updateParametersM(pars.matM);
    meser_2.updateParametersBeta(pars.beta);
    
    
    x0 = yini;  

    %create an array starting from 1 Jan
    yfull = zeros(T0+lastsamplingday,length(x0(1,:))); %make a full y array storing simulated data from day 1 in the year
    [t y] = ode23(@(t,x)odef_islmodjava(t,x, meser_2), times, x0);  
    clear('mepar_2');
    clear('meser_2');
    NDA = 0;
    if isfield(pars,'OutbreakNDA') == 1
      NDA = pars.OutbreakNDA;
    end
    
    yfull(1:T0,:) = repmat(x0,T0,1); %create initial data until T0
    if NDA < 0 
        yfull(1:T0-NDA,:) = repmat(x0,T0-NDA,1); %create initial data until T0-NDA
    end
    
    yfull(T0-NDA:T0-NDA+length(y(:,1))-1,:) = y; %save output into the array 
    T = t;
    a=1:pars.maxa;
    Y_posterior(i,:,:) = yfull; 
    
    %calculate peak date
    z = squeeze(yfull);
    strain = 1;
    [total_infecteds CInfecteds]= gen_total_infecteds(z, times+1, pars, strain, a);
    [total_infecteds1 CInfecteds1]= gen_total_infecteds(z, times+1, pars, strain, 1);
    [total_infecteds2 CInfecteds2]= gen_total_infecteds(z, times+1, pars, strain, 2);
    [total_infecteds3 CInfecteds3]= gen_total_infecteds(z, times+1, pars, strain, 3);
    [total_infecteds4 CInfecteds4]= gen_total_infecteds(z, times+1, pars, strain, 4);
    if ~exist('peak')
        tmp_peak = find(total_infecteds == max(total_infecteds));
        if length(tmp_peak) == 1
            peak(1) = tmp_peak;
        end
    else
        tmp_peak = find(total_infecteds == max(total_infecteds));
        if length(tmp_peak) == 1
            peak(end+1) = tmp_peak;
        end
    end
   

end
    if samplesize>1
        Y_mean = mean(Y_posterior);
        %Y_mean = reshape(Y_mean(1,:,:),[length(Y_mean(1,:,1)) length(Y_mean(1,1,:))]);
        Y_mean = squeeze(Y_mean);
    else
        Y_mean = Y_posterior;
    end
    Y_lb = quantile(Y_posterior,0.025);
    Y_lb = reshape(Y_lb(1,:,:),[length(Y_lb(1,:,1)) length(Y_lb(1,1,:))]);
    Y_ub = quantile(Y_posterior,0.975);
    Y_ub = reshape(Y_ub(1,:,:),[length(Y_ub(1,:,1)) length(Y_ub(1,1,:))]);

    
    
%% subplot figb
%figure;
agegroup = 1:4;
T_rel = times(1:365)+1;  %The relative days after 120d 
    
if exist('marker')
    [total_infecteds sero log CI_T1T2 total_age CI_T2 CI_T1] = plot_immune( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup, display, marker);
    %[total_infecteds sero log CI_T1T2 total_age] = plot_line( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup, marker);
    [total_infecteds sero log CI_T1T2 total_age rmse] = plot_Incidence( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup);
    strain = 1;
    [total_infecteds CInfecteds] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, agegroup); % plot CI
    [total_infecteds1 CInfecteds1] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 1); % plot CI
    [total_infecteds2 CInfecteds2] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 2); % plot CI
    [total_infecteds3 CInfecteds3] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 3); % plot CI
    [total_infecteds4 CInfecteds4] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 4); % plot CI
    plot(CInfecteds/total_age, 'g--');
else
    [total_infecteds sero log CI_T1T2 total_age CI_T2 CI_T1] = plot_immune( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup, display);
    %[total_infecteds sero log CI_T1T2 total_age] = plot_line( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup)
    [total_infecteds sero log CI_T1T2 total_age rmse] = plot_Incidence( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup);

    strain = 1;
    [total_infecteds CInfecteds] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, agegroup); % plot CI
    [total_infecteds1 CInfecteds1] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 1); % plot CI
    [total_infecteds2 CInfecteds2] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 2); % plot CI
    [total_infecteds3 CInfecteds3] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 3); % plot CI
    [total_infecteds4 CInfecteds4] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 4); % plot CI
    plot(CInfecteds/total_age, 'g');
end




new_peak = find(total_infecteds == max(total_infecteds)) + 120;
formatIn = 'mm/dd/yyyy';
ddd = datenum({'01/01/2009'},formatIn);
ddd1 = ddd+new_peak-1;
ddd2 = ddd+mean(peak)-1;

peak_lb = quantile(peak,0.025);
peak_ub = quantile(peak,0.975);
log = [log '\r\n' 'peak occurs at day '  num2str(mean(peak)) '(' num2str(peak_lb) ', '  num2str(peak_ub) ')'];
log = [log '\r\n' 'convert to date: ' datestr(ddd2)];
log = [log '\r\n' 'single_peak from average Incidence '  num2str(mean(new_peak)) '(' datestr(ddd1) ')'];
disp (log);

%[total_infecteds serolb loglb CI_T1T2_lb] = plot_line(  Y_lb(T0:T0+515-60,:), pars, T_rel, 1, agegroup, ':');
%[total_infecteds seroub logub CI_T1T2_ub] = plot_line(  Y_ub(T0:T0+515-60,:), pars, T_rel, 1, agegroup, ':'); 
%log = [log '\r\n' 'CI lb:' loglb];
%log = [log '\r\n' 'CI ub:' logub];

%log = [log '\r\n' 'CIT1T2 mean:' num2str(CI_T1T2)];
%log = [log '\r\n' 'CIT1T2 lb:' num2str(CI_T1T2_lb)];
%log = [log '\r\n' 'CIT1T2 ub:' num2str(CI_T1T2_ub)];

set(gca,'xlim',[1 30.5*10]);
set(gca,'XTick',[1:(365-60-1)/10:365-60]);
line([sample_time_K1 sample_time_K1], [0 1]);
line([sample_time_K2 sample_time_K2], [0 1]);

lastsamplingday = pars.SamplingLastDay + 90;
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
     
      
      %% Set Text labels 
      ax = axis;    % Current axis limits
      axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
      Yl = ax(3:4);  % Y-axis limits
      % Place the text labels
      Xt = [1:(365-60-1)/10:365-60]+15;
      t = text(Xt,Yl(1)*ones(1,length(Xt)),months(1:1:11,:));
      set(t,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45);

      % Remove the default labels
      set(gca,'XTickLabel','')

      % Get the Extent of each text object.  This
      % loop is unavoidable.
      for i = 1:length(t)
        ext(i,:) = get(t(i),'Extent');
      end
      
      % Determine the lowest point.  The X-label will be
      % placed so that the top is aligned with this point.
      LowYPoint = min(ext(:,2));

      % Place the axis label at this point
      XMidPoint = Xl(1)+abs(diff(Xl))/2;
      tl = text(XMidPoint,LowYPoint,'', 'VerticalAlignment','top','HorizontalAlignment','center');

      %set(gca(2),'YTickLabel',[0:0.01:0.02]);
%title(['agegroup' num2str(agegroup)]);
       set(gca, 'FontSize', 11.5);
       
       mTextBox = uicontrol('style','text');
       parentColor = [1 1 1];
       set(mTextBox,'String','Month', 'Position', [20 2   820 20],'FontSize',12,'backgroundcolor',parentColor);
end

end