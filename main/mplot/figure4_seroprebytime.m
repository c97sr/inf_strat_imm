function [FigH FigL log peak peak_lb peak_ub] = figure4_prebytime()
% 1) Plot the seroconversion between 2 models
% 2) calculate rmse
% plot the mean cumulative incidence
% H0 Standard threshold model
% H1 Titre model
% Todo list: 
% Plot a line for threshold model whose R0 is same as the titre model 
% check the line: calculateR0_fromPars(pars);

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



% plot titre model output (seroprevalence and cumulative incidence)
dat1 = load('out/p0e05/m1/ph1n1/20151024/mcmc_output_m1_final.mat');               % default model Titre.Full
resample_post = resamples(dat1, samplesize);
[rmse_1h rmse_1f hs1 hc1 r0_1 sero1 ci1 r0all_1] = plotI(resample_post, dat1, 1, [], 0);

% plot threshold model output (seroprevalence and cumulative incidence)
% with same samples 
dat2 = load('out/p0e05/m2/ph1n1/20151024/mcmc_output_m2_final.mat');               % default threshold model
resample_post2 = resamples(dat2, samplesize);
mark = ':';
[rmse_2h rmse_2f hs2 hc2 r0_2 sero2 ci2 r0all_2] = plotI(resample_post2, dat2, 1, mark, 1);

% plot threshold model fitting the predicted seroprevalence (seroprevalence and cumulative incidence)
mark = '-.';
%r_ratio = r0_1/r0_2;
%Try from 1.02->1.05
%r_ratio = 1.035; % adjust R0
r_ratio = 1.0325;
%r_ratio=1.01:0.005:1.05;
%sero_diff = 0;

%for id = 1:length(r_ratio)
    [rmse_3h rmse_3f hs3 hc3 r0_3 sero3 ci3 r0all_3] = plotI(resample_post2, dat2, 1, mark, 2, r_ratio);
%    sero_diff(id) = sero3 - sero1;
%    ci_diff(id) = ci3 - ci1;
%end
%ci_diff
legend([hs1 hc1 hs2 hc2 hs3],{'Titre model: Seroprevalence','Titre model: Cumulative Incidence','Threshold model: Seroprevalence','Threshold model: Cumulative Incidence','Threshold model(R_0): Seroprevalence'});
sample_time_K1 = 102;
sample_time_K2 = 235;
line([sample_time_K1 sample_time_K1], [0 1]);
line([sample_time_K2 sample_time_K2], [0 1]);
xlim([92 92+30.5*6]);
ylim([0 0.3]);
ylabel('Adjusted seroprevalence & cumulative incidence(%)');

disp(['rmse (until peak) for mean incidence model1:' num2str(rmse_1h)]);
disp(['rmse (until peak) for mean incidence model2:' num2str(rmse_2h)]);
disp(['rmse (full data) for mean incidence model1:' num2str(rmse_1f)]);
disp(['rmse (full data) for mean incidence model2:' num2str(rmse_2f)]);


function [posterior] = resamples(dat, samplesize)
posteriorTable = dat.PosteriorSamples;
pars = dat.par;
post = table2array(posteriorTable);

%resamples
burnIn = 1000;
total = height(posteriorTable(:,1))-burnIn;
idx = burnIn + round(rand(1, samplesize) * total);
posterior = post(idx,1:end-1);
end


function [RMSE_half RMSE_full hs hc R0_m seroprev CI_T2 R0] = plotI(post, dat, display, marker, mode, r_ratio)
    
formatIn = 'dd/mm/yyyy';
SDate = '01/05/2009';
numday_T1 = datenum('11/08/2009',formatIn)-datenum(SDate,formatIn)+1;
numday_T2 = datenum('22/12/2009',formatIn)-datenum(SDate,formatIn)+1;

    
% retrieve posterior from dat
posteriorTable = dat.PosteriorSamples;
pars = dat.par;
%post = table2array(posteriorTable);

% resamples
%samplesize = 4;
%burnIn = 1000;
%total = height(posteriorTable(:,1))-burnIn;
%idx = burnIn + round(rand(1, samplesize) * total);


posterior = post;

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
    
    

    
    
    % recalculate beta to fit R0
    if pars.maxi == 2
    R0 = [];
    for i = 1:samplesize;
    vars = posteriorTable.Properties.VariableNames;
    for p=1:length(vars)
        if strcmpi('LLH',vars(p))
        else
           [pars] = setParameters(pars,char(vars(p)),posterior(i,p));
           %[pars] = setParameters(pars,char(vars(p)),posterior(idx(i),p));
        end
    end
    %if strcmp(marker,'-.')
    R0(i) = calculateR0_fromPars(pars);
    end
    %%R0_ratio = 1.22./mean(R0);  
    %R0_ratio = 1.035;
    %r_ratio = R0_ratio;
    if strcmp(marker,'-.')
        %R0_ratio = 1.22./mean(R0);  
    end
    end
    
    
    %%Test the following line
    posterior_mean = mean(posterior);
    R0 = 0;
for i = 1:samplesize
    vars = posteriorTable.Properties.VariableNames;
    for p=1:length(vars)
        if strcmpi('LLH',vars(p))
        else
           [pars] = setParameters(pars,char(vars(p)),posterior(i,p));
           %[pars] = setParameters(pars,char(vars(p)),posterior_mean(idx(i),p));
           if pars.maxi == 2
           
           if strfind(char(vars(p)),'beta') 
               %if posterior(idx(i),p) > 3  
               if mode == 2  
                [pars] = setParameters(pars,char(vars(p)),posterior(i,p)*r_ratio);
               end
                 %R0_new(i) = calculateR0_fromPars(pars);
                 %[pars] = setParameters(pars,char(vars(p)),0.0494*0.97);
           end
           %if strfind(char(vars(p)),'ContFrac1')
           %    [pars] = setParameters(pars,char(vars(p)),5.0068);
           %end
           end
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
    

    R0(i) = calculateR0_fromPars(pars);

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
    %javaaddpath(pars.javapath); %set ./java as default dir
    javaaddpath('e:\Documents\Github\serodynamics\isltr\java\matlabjava.jar'); %use this one for calculating CI
    import matlabjava.*
    mepar_2 = matlabjava.ParametersSR;
    meser_2 = matlabjava.SerologySR;

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
    x0 = [x0 zeros(1,40) zeros(1,40)];
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
    agegroup = 1:4;
    T_rel = times(1:365)+1; 
    %[total_infecteds sero log CI_T1T2 total_age RMSE(i)] = plot_Incidence( yfull(T0:T0-60+515,:), pars, T_rel, 1, agegroup);
    [total_infecteds CInfecteds]= gen_total_infecteds(z, times+1, pars, strain, a);
    [total_infecteds1 CInfecteds1]= gen_total_infecteds(z, times+1, pars, strain, 1);
    [total_infecteds2 CInfecteds2]= gen_total_infecteds(z, times+1, pars, strain, 2);
    [total_infecteds3 CInfecteds3]= gen_total_infecteds(z, times+1, pars, strain, 3);
    [total_infecteds4 CInfecteds4]= gen_total_infecteds(z, times+1, pars, strain, 4);
    CIT2_tot(i) = CInfecteds(numday_T2+120);
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


    R0_m = mean(R0);
    R0_lb = quantile(R0,0.025);
    R0_ub = quantile(R0,0.975);

 if pars.maxi == 2
 % mean(R0_new)
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
    
if mode >= 1
    [total_infecteds sero log CI_T1T2 total_age CI_T2 CI_T1 CInfecteds Seroprev] = plot_immune( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup, 0, marker);
    %[total_infecteds sero log CI_T1T2 total_age] = plot_line( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup, marker);
    [total_infecteds sero log CI_T1T2 total_age RMSE_half RMSE_full] = plot_Incidence( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup);
    
    Serodiff = Seroprev - [0; Seroprev(1:end-1)];
    maxday = find(Serodiff == max(Serodiff));
    CIdiff = CInfecteds - [0; CInfecteds(1:end-1)];
    maxday_CI = find(CIdiff == max(CIdiff));
    tmp_peak_inc = find(total_infecteds == max(total_infecteds));
    tmp_peak_sero = find(Serodiff == max(Serodiff));
    seroprev = Seroprev(sample_time_K2);
    
    
    
    
    strain = 1;
    %[total_infecteds CInfecteds] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, agegroup); % plot CI
    %[total_infecteds1 CInfecteds1] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 1); % plot CI
    %[total_infecteds2 CInfecteds2] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 2); % plot CI
    %[total_infecteds3 CInfecteds3] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 3); % plot CI
    %[total_infecteds4 CInfecteds4] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 4); % plot CI
    
    hs = '';
    hc = '';
    if display == 1
    if mode == 1
    hs = plot(Seroprev/total_age,'b--');
    hc = plot(CInfecteds/total_age, 'g--');
    plot(maxday,Seroprev(maxday)/total_age,'bx');
    plot(maxday_CI,CInfecteds(maxday_CI)/total_age,'gx');
    elseif mode == 2
        hs = plot(Seroprev/total_age,'b--');
        hc = '';
        plot(maxday,Seroprev(maxday)/total_age,'bx');
    end
    end
        
else % mode = 0, titre model
    Y_sample = Y_posterior(1,:,:); 
    Y_sample = squeeze(Y_sample);
    [total_infecteds sero log CI_T1T2 total_age CI_T2 CI_T1 CInfecteds Seroprev CINaiveInfecteds CIImmInfecteds CINaiveInfectedsFull CIImmInfectedsFull]= plot_immune( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup, 0);
    Seroprev(sample_time_K2)/CInfecteds(sample_time_K2)                                                                                                                                                    
    seroprev = Seroprev(sample_time_K2);
   
    %[total_infecteds_l sero_l log_l CI_T1T2_l total_age_l CI_T2_l CI_T1_l CInfecteds_l]= plot_immune( Y_lb(T0:T0-60+515,:), pars, T_rel, 1, agegroup, 0);
    %[total_infecteds_u sero_u log_u CI_T1T2_u total_age_u CI_T2_u CI_T1_u CInfecteds_u]= plot_immune( Y_ub(T0:T0-60+515,:), pars, T_rel, 1, agegroup, 0); 
     
     
    %[total_infecteds sero log CI_T1T2 total_age] = plot_line( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup)
    [total_infecteds sero log CI_T1T2 total_age RMSE_half RMSE_full] = plot_Incidence( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup);

    
    Serodiff = Seroprev - [0; Seroprev(1:end-1)];
    maxday = find(Serodiff == max(Serodiff));
    CIdiff = CInfecteds - [0; CInfecteds(1:end-1)];
    maxday_CI = find(CIdiff == max(CIdiff));
    tmp_peak_inc = find(total_infecteds == max(total_infecteds));
    tmp_peak_sero = find(Serodiff == max(Serodiff));
    
    strain = 1;
    [total_infecteds CInfecteds] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, agegroup); % plot CI
    [total_infecteds1 CInfecteds1] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 1); % plot CI
    [total_infecteds2 CInfecteds2] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 2); % plot CI
    [total_infecteds3 CInfecteds3] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 3); % plot CI
    [total_infecteds4 CInfecteds4] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 4); % plot CI
    hs = plot(Seroprev/total_age,'b');
    hc = plot(CInfecteds/total_age, 'g');
    plot(maxday,Seroprev(maxday)/total_age,'bo');
    plot(maxday_CI,CInfecteds(maxday_CI)/total_age,'go');
    
    %plot(71:7*1:71+7*1*(length(Seroprev_obs)-1),Seroprev_obs,'o');
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


lastsamplingday = pars.SamplingLastDay + 90;
Xl = [1 lastsamplingday];

months = [
          '   ';
          '   ';
          '   ';
          'Aug';
          'Sep';
          'Oct';
          'Nov';
          'Dec';
          'Jan';
          '   ';
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
       set(mTextBox,'String','Month', 'Position', [20 2   920 20],'FontSize',12,'backgroundcolor',parentColor);
end

end