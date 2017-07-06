function [FigH seroprv_diff] = figureS6_seroprevalence( PosteriorSamples, pars, burnIn, samplesize) 
% Summary of the function figureNew
% The function to produce a new figure for seroprevalence.
% Written by Sean Yuan (hyuan@imperial.ac.uk) 
% Plot all age group information.


global proj Antibody;
init_collect = 1;
second_collect = 2;
third_collect = 3;
k = 2;
Ab = Antibody;
posterior = table2array(PosteriorSamples);

if exist('samplesize') == 0
    samplesize = 3;
end
if exist('burnIn') == 0
    burnIn = 1000;
end
post = mean(posterior(burnIn:end,:));
total = length(posterior(:,1))-burnIn;
idx = burnIn + round(rand(1, samplesize) * total);
seroconv = [];

for i=1:samplesize
%model:3
    vars = PosteriorSamples.Properties.VariableNames;
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
lastsamplingday = pars.SamplingLastDay;

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
  
    mean_sample_time_K1 = round(meanKdays(1));
    mean_sample_time_K2 = round(meanKdays(2));
    sample_time_K1 = pars.Antibody.K(1).numdays - T0;
    sample_time_K2 = pars.Antibody.K(2).numdays - T0;
    sample_time_K1_age = [];
    sample_time_K2_age = [];
    times = 0:1:lastsamplingday;
    sample_size_K1 = Ab.K(1).samplesize;
    for i1=1:pars.maxa
         lage = pars.ages(i1,1);
         uage = pars.ages(i1,2);
         idx_a = find(lage<=Ab.K(1).age & Ab.K(1).age<uage);
         sample_size_K1(end+1) = length(find(lage<=Ab.K(1).age & Ab.K(1).age<uage));
         sample_time_K1_age(i1).time = Ab.K(1).numdays(idx_a)-T0;
    end
    sample_size_K2 = Ab.K(2).samplesize;
    for i2=1:pars.maxa
         lage = pars.ages(i2,1);
         uage = pars.ages(i2,2);
         idx_a = find(lage<=Ab.K(2).age & Ab.K(2).age<uage);
         sample_size_K2(end+1) = length(find(lage<=Ab.K(2).age & Ab.K(2).age<uage));
         sample_time_K2_age(i2).time = Ab.K(2).numdays(idx_a) - T0;
    end

    %run simulation
    %initialize objects
    daysshift = 0;
    daysshift = pars.OutbreakNDA;
    
    %javaaddpath e:\workspace\MyJavaProject\bin\matlabjava.jar; %use jave working dir 
    javaaddpath(pars.javapath); %set ./java as default dir
    mepar_4 = matlabjava.Parameters;
    meser_4 = matlabjava.Serology;
    % set parameters
    meser_4.setParameters(mepar_4);
    meser_4.updateParameters('wan',pars.wan);
    meser_4.updateParameters('s0_imm', pars.s0_imm);
    meser_4.updateParameters('maxi', pars.maxi);
    meser_4.updateParametersG(pars.arrg);
    meser_4.updateParametersH(pars.arrh);
    meser_4.updateParametersM(pars.matM);
    meser_4.updateParametersBeta(pars.beta);  
    x0 = yini;  
    
    [t y] = ode23(@(t,x)odef_islmodjava(t,x, meser_4), times, x0);  
   
    %[t y] = ode23(@(t,x)odef_islmod(t,x,pars), times, x0);
    clear 'mepar_4';
    clear 'meser_4';
    T = t;
      
    %seroconversion
    atotal=1:pars.maxa; %%%
    %total age group
    %%%%%%%
    %%%%%%%
    %% Sampling model output should be consistent to observed data
    %% 20151006
    %Xout_k1 = mean(retrieve_histogram(y, pars, times, sample_time_K1+daysshift, atotal)); % model output
    %Xout_k2 = mean(retrieve_histogram(y, pars, times, sample_time_K2+daysshift, atotal)); % model output
    Xout_k1 = retrieve_histogram(y, pars, times, mean_sample_time_K1+daysshift, atotal); % model output
    Xout_k2 = retrieve_histogram(y, pars, times, mean_sample_time_K2+daysshift, atotal); % model output
    Yout_k0 = retrieve_histogram(yini, pars, times(1), 1, atotal);
    Yout_k1 = retrieve_histogram(yini_k1, pars, times(1), 1, atotal);     % observed data
    Yout_k2 = retrieve_histogram(yini_k2, pars, times(1), 1, atotal);     % observed data
    if pars.maxi > 2
    seroprev_total(i,1) = sum(Yout_k0(4:end)); %Pre-existing antibody
    seroprev_total(i,2) = sum(Yout_k1(4:end)); %T1 Observed
    seroprev_total(i,3) = sum(Yout_k2(4:end)); %T2 Observed
    seroprev_total(i,4) = sum(Xout_k1(4:end)); %T1 Model output
    seroprev_total(i,5) = sum(Xout_k2(4:end)); %T2 Model output
    else
    seroprev_total(i,1) = sum(Yout_k0(2)); %Pre-existing antibody
    seroprev_total(i,2) = sum(Yout_k1(2)); %T1 Observed
    seroprev_total(i,3) = sum(Yout_k2(2)); %T2 Observed
    seroprev_total(i,4) = sum(Xout_k1(2)); %T1 Model output
    seroprev_total(i,5) = sum(Xout_k2(2)); %T2 Model output
    end
    seroconv(i) = seroprev_total(i,5)-seroprev_total(i,4);
    clear Xout_k1;
    clear Xout_k2;
    clear Yout_k0;
    clear Yout_k1;
    clear Yout_k2;
    %mean titre changes
    for a=1:pars.maxa
%        sample_time_K2_age(a).time
    %Xout_k1(a,:) = mean(retrieve_histogram(y, pars, times, sample_time_K1_age(a).time+daysshift, a)); % model output
    %Xout_k2(a,:) = mean(retrieve_histogram(y, pars, times, sample_time_K2_age(a).time+daysshift, a)); % model output
    Xout_k1(a,:) = retrieve_histogram(y, pars, times, mean_sample_time_K1+daysshift, a); % model output
    Xout_k2(a,:) = retrieve_histogram(y, pars, times, mean_sample_time_K2+daysshift, a); % model output
    Yout_k0(a,:) = retrieve_histogram(yini, pars, times(1), 1, a);
    Yout_k1(a,:) = retrieve_histogram(yini_k1, pars, times(1), 1, a); % observed data
    Yout_k2(a,:) = retrieve_histogram(yini_k2, pars, times(1), 1, a); % observed data
    
    agetitres(a).Yout_k1(a,:) = retrieve_histogram(yini_k1, pars, times(1), 1, a); % observed data
    agetitres(a).Yout_k2(a,:) = retrieve_histogram(yini_k2, pars, times(1), 1, a); % observed data
    
    Xout_k1; % observed data
    Xout_k2; % observed data 
    if pars.maxi > 2
        seroprev(a).age(i,1) = sum(Yout_k0(a,4:end)); %T0 
        seroprev(a).age(i,2) = sum(Yout_k1(a,4:end)); %T1 Observed
        seroprev(a).age(i,3) = sum(Yout_k2(a,4:end)); %T2 Observed
        seroprev(a).age(i,4) = sum(Xout_k1(a,4:end)); %T1 Model output
        seroprev(a).age(i,5) = sum(Xout_k2(a,4:end)); %T2 Model output
    else
        seroprev(a).age(i,1) = sum(Yout_k0(a,2)); %T0 
        seroprev(a).age(i,2) = sum(Yout_k1(a,2)); %T1 Observed
        seroprev(a).age(i,3) = sum(Yout_k2(a,2)); %T2 Observed
        seroprev(a).age(i,4) = sum(Xout_k1(a,2)); %T1 Model output
        seroprev(a).age(i,5) = sum(Xout_k2(a,2)); %T2 Model output
    end
    end
    
end

    seroprv_diff.mean = mean(seroconv);
    seroprv_diff.lb = quantile(seroconv,0.025);
    seroprv_diff.ub = quantile(seroconv,0.975);
    %Seroprevalence
    %baseline = seroprev_total(1,1);
    seroprev_mean = mean(seroprev_total(:,:));   % T0, initial, followup, modelT1, modelT2
    seroprev_lb = quantile(seroprev_total,0.025);% lower bound;
    seroprev_ub = quantile(seroprev_total,0.975);% upper bound;

    
    
    pbin = seroprev_mean(1,2); %round1 prevalence
    %find lower bound
    y0 = pbin; %inital value
    sizeage = sample_size_K1(1); % total:523
    options = optimoptions('fsolve','Display','off'); % Option to display output
    lb = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.975','x','pbin','sizeage'),y0,[options],pbin,sizeage);
    seroprev_lb(1,2) = lb;
    %find upper bound
    ub = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.025','x','pbin','sizeage'),y0,[options],pbin,sizeage);
    seroprev_ub(1,2) = ub;
    
    pbin = seroprev_mean(1,3); %round2 prevalence
    %find lower bound
    y0 = pbin; %inital value
    sizeage = sample_size_K2(1);
    lb = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.975','x','pbin','sizeage'),y0,[options],pbin,sizeage);
    seroprev_lb(1,3) = lb;
    %find upper bound
    ub = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.025','x','pbin','sizeage'),y0,[options],pbin,sizeage);
    seroprev_ub(1,3) = ub;
    
 
    for a=1:pars.maxa
        %Seroprevalence
        %baseline_age = seroprev(a).age(1,1);
        seroprev_mean(a+1,1:5) = mean(seroprev(a).age);           % - baseline_age;
        seroprev_lb(a+1,1:5) = quantile(seroprev(a).age,0.025);  % - lower bound
        seroprev_ub(a+1,1:5) = quantile(seroprev(a).age,0.975); % - upper bound;
        
        pbin = seroprev_mean(a+1,2); %round1 prevalence
        %find lower bound
        y0 = pbin; %inital value
        sizeage = sample_size_K1(1+a); % sample size for each age group 
        lb = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.975','x','pbin','sizeage'),y0,[options],pbin,sizeage);
        seroprev_lb(a+1,2) = lb;
        %find upper bound
        ub = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.025','x','pbin','sizeage'),y0,[options],pbin,sizeage);
        seroprev_ub(a+1,2) = ub;
        
        pbin = seroprev_mean(a+1,3); %round2 prevalence
        %find lower bound
        y0 = pbin; %inital value
        sizeage = sample_size_K2(1+a);  % sample size for each age group 
        lb = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.975','x','pbin','sizeage'),y0,[options],pbin,sizeage);
        seroprev_lb(a+1,3) = lb;
        %find upper bound
        ub = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.025','x','pbin','sizeage'),y0,[options],pbin,sizeage);
        seroprev_ub(a+1,3) = ub;
    end
    
 %Plot Seroprevalence
 seroprev_mean(:,1) = [];
 sero_all = seroprev_mean(1,:);
 sero_all(2) = seroprev_mean(1,3);
 sero_all(3) = seroprev_mean(1,2);
 sero_age = seroprev_mean(2:end,:);
 sero_age(:,2) = seroprev_mean(2:end,3);
 sero_age(:,3) = seroprev_mean(2:end,2);
 
 
              %age1 age2 age3 age4
Strain1_Mean= [sero_all(1) sero_age(:,1)']; % data Y baseline
Strain2_Mean= [sero_all(2) sero_age(:,2)']; % data X baseline
Strain3_Mean= [sero_all(3) sero_age(:,3)']; % data Y followup
Strain4_Mean= [sero_all(4) sero_age(:,4)']; % data X followup

seroprev_ub(:,1) = [];
seroub_all = seroprev_ub(1,:);
seroub_all(2) = seroprev_ub(1,3);
seroub_all(3) = seroprev_ub(1,2);
seroub_age = seroprev_ub(2:end,:);
seroub_age(:,2) = seroprev_ub(2:end,3);
seroub_age(:,3) = seroprev_ub(2:end,2);

seroprev_lb(:,1) = [];
serolb_all = seroprev_lb(1,:);
serolb_all(2) = seroprev_lb(1,3);
serolb_all(3) = seroprev_lb(1,2);
serolb_age = seroprev_lb(2:end,:);
serolb_age(:,2) = seroprev_lb(2:end,3);
serolb_age(:,3) = seroprev_lb(2:end,2);

Strain1_lb= Strain1_Mean - [serolb_all(1) serolb_age(:,1)'];  % data Y baseline
Strain2_lb= Strain2_Mean - [serolb_all(2) serolb_age(:,2)'];
Strain3_lb= Strain3_Mean - [serolb_all(3) serolb_age(:,3)'];
Strain4_lb= Strain4_Mean - [serolb_all(4) serolb_age(:,4)'];

Strain1_ub= [seroub_all(1) seroub_age(:,1)'] - Strain1_Mean;
Strain2_ub= [seroub_all(2) seroub_age(:,2)'] - Strain2_Mean;
Strain3_ub= [seroub_all(3) seroub_age(:,3)'] - Strain3_Mean;
Strain4_ub= [seroub_all(4) seroub_age(:,4)'] - Strain4_Mean;

FigH = figure;
%%%Double sides errors
%barwitherr(cat(3,[Strain1_lb' Strain2_lb' Strain3_lb'...
%    Strain4_lb'],[Strain1_ub' Strain2_ub' Strain3_ub'...
%    Strain4_ub']),[1 2 3 4 5],[Strain1_Mean' Strain2_Mean'...
%    Strain3_Mean' Strain4_Mean'],'LineWidth',0.1,...
%    'BarWidth',1)
[var hbar] = barwitherr(cat(3,[Strain1_lb' Strain2_lb' Strain3_lb'...
    Strain4_lb'],[Strain1_ub' Strain2_ub' Strain3_ub'...
    Strain4_ub']),[1 2 3 4 5],[Strain1_Mean' Strain2_Mean'...
    Strain3_Mean' Strain4_Mean'],'LineWidth',0.1,...
    'BarWidth',1);
set(hbar(1),'FaceColor',[0.2 0.3 0.8] ); % Observed Baseline
set(hbar(2),'FaceColor',[0.3 0.3 0.3] ); % Model Baseline
set(hbar(3),'FaceColor',[0.6 0.8 0.9] ); % Observed Follow-up
set(hbar(4),'FaceColor',[0.7 0.7 0.7] ); % Model Follow-up

set(gca,'XTickLabel',{'All', '<20', '20-39', '40-64', '\geq65'});
set(gca,'YTickLabel',0:10:100);
legend('Baseline (observed)','Baseline (model output)','Follow-up (observed)','Follow-up (model output)')
disp 'keep walking';
xlabel('Ages (yrs)', 'FontSize', 12);
ylabel('Seroprevalence (%)', 'FontSize', 12);
%set(get(gca, 'YLabel'), 'String', 'Seroprevalence (%)');
end