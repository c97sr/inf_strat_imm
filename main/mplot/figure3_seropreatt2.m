function [FigH FigL log peak peak_lb peak_ub] = figure3_seropreatt2()
% 1) Plot the changes in seroprevalence with 3 diffrent pre-existing
% immunity


% Todo list: 


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
    samplesize = 3;
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
%set(gca,'xtick',[]);
%hold;



% plot titre model output
%dat1 = load('out/p0e05/m1/ph1n1/20151024/mcmc_output_m1_final.mat');               % default model Titre.Full
%resample_post = resamples(dat1, samplesize);
%mode = 0; %0 default immunity, 1 low immunity in elderly
%[seroprev_ch_mean_0 seroprev_ch_lb_0 seroprev_ch_ub_0] = plotI(resample_post, dat1, 1, [], mode);
%mode = 1;
%[seroprev_ch_mean_1 seroprev_ch_lb_1 seroprev_ch_ub_1] = plotI(resample_post, dat1, 1, [], mode);
%mode = 2;
%[seroprev_ch_mean_2 seroprev_ch_lb_2 seroprev_ch_ub_2] = plotI(resample_post, dat1, 1, [], mode);
%intv = 0.15;
%X = [1-intv 2-intv 3-intv 4-intv; 1 2 3 4; 1+intv 2+intv 3+intv 4+intv]';
%C = [seroprev_ch_mean_0; seroprev_ch_mean_1; seroprev_ch_mean_2]';
%L = [seroprev_ch_mean_0-seroprev_ch_lb_0; seroprev_ch_mean_1-seroprev_ch_lb_1; seroprev_ch_mean_2-seroprev_ch_lb_2]';
%U = [seroprev_ch_ub_0-seroprev_ch_mean_0; seroprev_ch_ub_1-seroprev_ch_mean_1; seroprev_ch_ub_2-seroprev_ch_mean_2]';
%errorbar(X,C,L,U,'.');
%xlabel('Age groups');
%ylabel('Changes in seroprevalence after the wave(%)');
%legend('Titre.A','Titre.A(reduced pre-imm in elderly)','Titre.A(increased pre-imm in elderly)');

% plot threshold model output (seroprevalence and cumulative incidence)
% with same samples 
%FigH = figure;
%set(FigH, 'Position', [200, 600, 816, 480]);
%dat2 = load('out/fix/m5/ph1n1/20160304/mcmc_output_m5_final.mat');               % default threshold model
%resample_post2 = resamples(dat2, samplesize);
%mode = 0; %0 default immunity, 1 low immunity in elderly
%[seroprev_ch_mean_0 seroprev_ch_lb_0 seroprev_ch_ub_0] = plotI(resample_post2, dat2, 1, [], mode);
%mode = 1;
%[seroprev_ch_mean_1 seroprev_ch_lb_1 seroprev_ch_ub_1] = plotI(resample_post2, dat2, 1, [], mode);
%mode = 2;
%[seroprev_ch_mean_2 seroprev_ch_lb_2 seroprev_ch_ub_2] = plotI(resample_post2, dat2, 1, [], mode);
%intv = 0.15;
%X = [1-intv 2-intv 3-intv 4-intv; 1 2 3 4; 1+intv 2+intv 3+intv 4+intv]';
%C = [seroprev_ch_mean_0; seroprev_ch_mean_1; seroprev_ch_mean_2]';
%L = [seroprev_ch_mean_0-seroprev_ch_lb_0; seroprev_ch_mean_1-seroprev_ch_lb_1; seroprev_ch_mean_2-seroprev_ch_lb_2]';
%U = [seroprev_ch_ub_0-seroprev_ch_mean_0; seroprev_ch_ub_1-seroprev_ch_mean_1; seroprev_ch_ub_2-seroprev_ch_mean_2]';
%errorbar(X,C,L,U,'.');
xlabel('Age groups');
ylabel('Changes in seroprevalence after the wave(%)');
legend('Titre.D','Titre.D(reduced pre-imm in elderly)','Titre.D(increased pre-imm in elderly)');


% plot titre model output
dat1 = load('out/boost/m3/ph1n1/20160304/mcmc_output_m3_final.mat');               % default model Titre.Full
resample_post = resamples(dat1, samplesize);
mode = 0; %0 default immunity, 1 low immunity in elderly
[seroprev_ch_mean_0 seroprev_ch_lb_0 seroprev_ch_ub_0] = plotI(resample_post, dat1, 1, [], mode);
mode = 1;
[seroprev_ch_mean_1 seroprev_ch_lb_1 seroprev_ch_ub_1] = plotI(resample_post, dat1, 1, [], mode);
mode = 2;
[seroprev_ch_mean_2 seroprev_ch_lb_2 seroprev_ch_ub_2] = plotI(resample_post, dat1, 1, [], mode);
intv = 0.15;
X = [1-intv 2-intv 3-intv 4-intv; 1 2 3 4; 1+intv 2+intv 3+intv 4+intv]';
C = [seroprev_ch_mean_0; seroprev_ch_mean_1; seroprev_ch_mean_2]';
L = [seroprev_ch_mean_0-seroprev_ch_lb_0; seroprev_ch_mean_1-seroprev_ch_lb_1; seroprev_ch_mean_2-seroprev_ch_lb_2]';
U = [seroprev_ch_ub_0-seroprev_ch_mean_0; seroprev_ch_ub_1-seroprev_ch_mean_1; seroprev_ch_ub_2-seroprev_ch_mean_2]';
errorbar(X,C,L,U,'.');
xlabel('Age groups');
ylabel('Changes in seroprevalence after the wave(%)');
legend('Titre.B','Titre.B(reduced pre-imm in elderly)','Titre.B(increased pre-imm in elderly)');


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


function [seroprev_a_change_mean seroprev_a_change_lb seroprev_a_change_ub] = plotI(post, dat, display, marker, mode)
    
formatIn = 'dd/mm/yyyy';
SDate = '01/05/2009';
numday_T1 = datenum('11/08/2009',formatIn)-datenum(SDate,formatIn)+1;
numday_T2 = datenum('22/12/2009',formatIn)-datenum(SDate,formatIn)+1;

    
% retrieve posterior from dat
posteriorTable = dat.PosteriorSamples;
pars = dat.par;

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
        %[yini age_arr] = make_ics_naive( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu);
    end
    
    
    %setup initial immunity
    pre_sero = 0;
    pre_sero4 = 0;

    if mode == 2
        pars.inittitres_flag = 3 % with other initial immunity setting
        %obs=0.4:   40% seroconversion
        %obs=0.2:   20%  seroconversion
        seroconvert_obs = 0.08;
        prev = [0.5 0.5 0.5 0.5];
        naive = 1 - sum(prev)*seroconvert_obs;
        pars.init_prev = [naive prev*seroconvert_obs zeros(1,pars.maxi-length(prev)-1)];
        pre_sero = 0.02*0.867+seroconvert_obs*0.133;
        pre_sero_a = [0.02 0.02 0.02 seroconvert_obs];
    end
    if mode == 1
        pars.inittitres_flag = 3 % with other initial immunity setting
        %obs=0.4:   40% seroconversion
        %obs=0.02:  2%  seroconversion
        seroconvert_obs = 0.02;
        prev = [0.5 0.5 0.5 0.5];
        naive = 1 - sum(prev)*seroconvert_obs;
        pars.init_prev = [naive prev*seroconvert_obs zeros(1,pars.maxi-length(prev)-1)];
        pre_sero = 0.02*0.867+seroconvert_obs*0.133;
        pre_sero_a = [0.02 0.02 0.02 seroconvert_obs];
    end
    if mode == 0
        %don't need to do anythign
        pre_sero = 0.0227;
        pre_sero_a = [0.02 0.02 0.02 0.04];
    end
    [yini age_arr] = make_ics_naive( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu);



    
    
    %%Test the following line
    posterior_mean = mean(posterior);
    R0 = 0;
for i = 1:samplesize
    vars = posteriorTable.Properties.VariableNames;
    for p=1:length(vars)
        if strcmpi('LLH',vars(p))
        else
           [pars] = setParameters(pars,char(vars(p)),posterior(i,p));
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
    if pars.model == 5
     sample_time_K2 = round(meanKdays(2)) + 60;
    end
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
    agegroup = 1:4;
    T_rel = times(1:365)+1; 


end


    if samplesize>1
        Y_mean = mean(Y_posterior);
        Y_lb = quantile(Y_posterior,0.025);
        Y_ub = quantile(Y_posterior,0.975);   
        Y_mean = squeeze(Y_mean);
        Y_lb = squeeze(Y_lb);
        Y_ub = squeeze(Y_ub);
    else
        Y_mean = Y_posterior;
        Y_lb = Y_mean;
        Y_ub = Y_mean;
    end
    %adjust the starting date
    Y_meanT0 = Y_mean(T0:T0-60+515,:);
    Y_lbT0 = Y_lb(T0:T0-60+515,:);
    Y_ubT0 = Y_ub(T0:T0-60+515,:);
    %Y_lb = reshape(Y_lb(1,:,:),[length(Y_lb(1,:,1)) length(Y_lb(1,1,:))]);
    %Y_ub = reshape(Y_ub(1,:,:),[length(Y_ub(1,:,1)) length(Y_ub(1,1,:))]);

    
%% subplot figb
%figure;
agegroup = 1:4;
T_rel = times(1:365)+1;  %The relative days after 120d 
    
    %average seroprevalence
    [total_infecteds sero log CI_T1T2 total_age CI_T2 CI_T1 CInfecteds Seroprev CINaiveInfecteds CIImmInfecteds CINaiveInfectedsFull CIImmInfectedsFull CITreInfectedsFull CInfecteds_age]= plot_immune( Y_meanT0, pars, T_rel, 1, agegroup, 0);
    sero1 = squeeze(sero(1,sample_time_K2,:));
    sero2 = squeeze(sero(2,sample_time_K2,:));
    sero3 = squeeze(sero(3,sample_time_K2,:));
    sero4 = squeeze(sero(4,sample_time_K2,:));
    seroprev1_change = sum(sero1(4:end))/sum(sero1)- pre_sero_a(1); 
    seroprev2_change = sum(sero2(4:end))/sum(sero2)- pre_sero_a(2);
    seroprev3_change = sum(sero3(4:end))/sum(sero3)- pre_sero_a(3); 
    seroprev4_change = sum(sero4(4:end))/sum(sero4)- pre_sero_a(4);
    seroprev_a_change_mean = [seroprev1_change seroprev2_change seroprev3_change seroprev4_change];
    
    %seroprevalence lb bound
    [total_infecteds sero log CI_T1T2 total_age CI_T2 CI_T1 CInfecteds Seroprev CINaiveInfecteds CIImmInfecteds CINaiveInfectedsFull CIImmInfectedsFull CITreInfectedsFull CInfecteds_age]= plot_immune( Y_lbT0, pars, T_rel, 1, agegroup, 0);
    sero1 = squeeze(sero(1,sample_time_K2,:));
    sero2 = squeeze(sero(2,sample_time_K2,:));
    sero3 = squeeze(sero(3,sample_time_K2,:));
    sero4 = squeeze(sero(4,sample_time_K2,:));
    seroprev1_change_lb = sum(sero1(4:end))/sum(sero1)- pre_sero_a(1); 
    seroprev2_change_lb = sum(sero2(4:end))/sum(sero2)- pre_sero_a(2);
    seroprev3_change_lb = sum(sero3(4:end))/sum(sero3)- pre_sero_a(3); 
    seroprev4_change_lb = sum(sero4(4:end))/sum(sero4)- pre_sero_a(4);
    seroprev_a_change_lb = [seroprev1_change_lb seroprev2_change_lb seroprev3_change_lb seroprev4_change_lb];
    %seroprevalence ub bound
    [total_infecteds sero log CI_T1T2 total_age CI_T2 CI_T1 CInfecteds Seroprev CINaiveInfecteds CIImmInfecteds CINaiveInfectedsFull CIImmInfectedsFull CITreInfectedsFull CInfecteds_age]= plot_immune( Y_ubT0, pars, T_rel, 1, agegroup, 0);
    sero1 = squeeze(sero(1,sample_time_K2,:));
    sero2 = squeeze(sero(2,sample_time_K2,:));
    sero3 = squeeze(sero(3,sample_time_K2,:));
    sero4 = squeeze(sero(4,sample_time_K2,:));
    seroprev1_change_ub = sum(sero1(4:end))/sum(sero1)- pre_sero_a(1); 
    seroprev2_change_ub = sum(sero2(4:end))/sum(sero2)- pre_sero_a(2);
    seroprev3_change_ub = sum(sero3(4:end))/sum(sero3)- pre_sero_a(3); 
    seroprev4_change_ub = sum(sero4(4:end))/sum(sero4)- pre_sero_a(4);
    seroprev_a_change_ub = [seroprev1_change_ub seroprev2_change_ub seroprev3_change_ub seroprev4_change_ub];

    
end

end