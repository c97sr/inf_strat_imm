function [FigH FigL log peak peak_lb peak_ub] = cal_rmse()
% 1) Plot the seroconversion between 2 models
% 2) calculate rmse
% plot the mean cumulative incidence
% H0 Standard threshold model
% H1 Titre model


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
    samplesize = 30;
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
%dat1 = load('out/imm/m4/ph1n1/20160303/mcmc_output_m4.mat');
%dat1 = load('out/p0e05/m1.5/ph1n1/20151027/mcmc_output_m1.5_final.mat');
[rmse_1h rmse_1f] = plotI(dat1, 1);

dat2 = load('out/p0e05/m2/ph1n1/20151024/mcmc_output_m2_final.mat');
%dat2 = load('out/t20/m2/ph1n1/20160313/mcmc_output_m2_final.mat');
%dat2 = load('out/p0e05/m2.2/ph1n1/20151024/mcmc_output_m2.2_final.mat');
%dat2 = load('out/imm/m4/ph1n1/20160303/mcmc_output_m4.mat');
%dat2 = load('out/boost/m3/ph1n1/20160303/mcmc_output_m3_final.mat');
mark = ':';
[rmse_2h rmse_2f] = plotI(dat2, 1, mark);
disp(['rmse (until peak) for mean incidence model1:' num2str(rmse_1h)]);
disp(['rmse (until peak) for mean incidence model2:' num2str(rmse_2h)]);

disp(['rmse (full data) for mean incidence model1:' num2str(rmse_1f)]);
disp(['rmse (full data) for mean incidence model2:' num2str(rmse_2f)]);
function [RMSE_h RMSE_f] = plotI(dat, display, marker)
    
posteriorTable = dat.PosteriorSamples;
pars = dat.par;
post = table2array(posteriorTable);
%posterior = repmat(mean(post),3,1);
%samplesize = 3;
%idx = 1:3;

%resamples
samplesize = 100;
burnIn = 1000;
total = height(posteriorTable(:,1))-burnIn;
idx = burnIn + round(rand(1, samplesize) * total);
posterior = post;



for i = 1:samplesize
    vars = posteriorTable.Properties.VariableNames;
    for p=1:length(vars)
        if strcmpi('LLH',vars(p))
        else
           [pars] = setParameters(pars,char(vars(p)),posterior(idx(i),p));
           if strfind(char(vars(p)),'beta') & length(vars)<5
               %if posterior(idx(i),p) > 3  
               %[pars] = setParameters(pars,char(vars(p)),posterior(idx(i),p)+0.019*posterior(idx(i),p));
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
    %meser_2.updateParameters('wan',0);
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
    agegroup = 1:4;
    T_rel = times(1:365)+1; 
    %[total_infecteds sero log CI_T1T2 total_age RMSE(i)] = plot_Incidence( yfull(T0:T0-60+515,:), pars, T_rel, 1, agegroup);
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

  RMSE_h=[];
  RMSE_f=[];
   for i=1:samplesize
        Y_mean = Y_posterior(i,:,:);
        %Y_mean = reshape(Y_mean(1,:,:),[length(Y_mean(1,:,1)) length(Y_mean(1,1,:))]);
        Y_mean = squeeze(Y_mean);
        [total_infecteds sero log CI_T1T2 total_age CI_T2 CI_T1] = plot_immune( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup, display);
        [total_infecteds sero log CI_T1T2 total_age RMSE_half RMSE_full] = plot_Incidence( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup);
        RMSE_f(i) = RMSE_full;
        RMSE_h(i) = RMSE_half;
   end

end
end

