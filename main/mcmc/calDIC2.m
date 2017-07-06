function [DIC Dbar PD] = calDIC2(filename)
%calDIC Summary of this function goes here
%This is the correct version to be used to generate the DIC in the manuscript
%05/07/2016
% Written by Sean Yuan (hyuan@imperial.ac.uk)
    %load(filename);
    
    filename = 'out/p0e05/m1/ph1n1/20151024/mcmc_output_m1_final.mat';
    %filename =
    %'out/p0e05/m1.12/ph1n1/20151024/mcmc_output_m1.12_final.mat'; %seed=3
    %filename = 'out/p0e05/m1.13/ph1n1/20151024/mcmc_output_m1.13_final.mat'; %seed=30
    %filename = 'out/p0e05/m1.14/ph1n1/20151024/mcmc_output_m1.14_final.mat'; %seed=100
    %filename = 'out/boost/m3/ph1n1/20160304/mcmc_output_m3_final.mat';
    %filename = 'out/imm/m4/ph1n1/20160304/mcmc_output_m4_final.mat';
    %filename = 'out/fix/m5/ph1n1/20160304/mcmc_output_m5_final.mat';
    
    %filename = 'out/p0e05/m2/ph1n1/20151024/mcmc_output_m2_final.mat';
    %filename = 'out/p0e05/m2.12/ph1n1/20151024/mcmc_output_m2.12_final.mat'; %seed=3
    %filename = 'out/p0e05/m2.13/ph1n1/20151024/mcmc_output_m2.13_final.mat'; %seed=30
    %filename = 'out/p0e05/m2.14/ph1n1/20151024/mcmc_output_m2.14_final.mat'; %seed=100
    %filename = 'out/classic/m6/ph1n1/20160307/mcmc_output_m6_final.mat';
    %filename = 'out/t20/m2/ph1n1/20160313/mcmc_output_m2_final.mat';
    
    %load('out/p0e05/m1.14/ph1n1/20151024/mcmc_output_m1.14_final.mat');
    %load('out/p0e05/m1.15/ph1n1/20151021/mcmc_output_m1.15_final.mat');
    %filename = 'out/p0e05/m2/ph1n1/20151024/mcmc_output_m2_final.mat';
    %load('out/imm/m4/ph1n1/20160303/mcmc_output_m4.mat');
    %filename = 'out/imm/m4/ph1n1/20160303/mcmc_output_m4.mat';
    %load('out/boost/m3/ph1n1/20160303/mcmc_output_m3.mat');
    load(filename);
    sid = regexp(filename,'/');
    outdir = filename(1:sid(end));
    Xbar = mean(table2array(PosteriorSamples));
    thetabar = Xbar(1:end-1);
    LLHbar = Xbar(end);
    Dbar = -2*LLHbar;
    Antibody = par.Antibody;
    SampleSize = Antibody.samplesize;
    Abl = Antibody.Abl;
    age = Antibody.age;
    Abl(find(Abl>9)) = 9; %substitute Ab level >9

    
    %Only 0 and 1 immune status
    %par.maxt = 1;
    AbCut = 3; %AbCut=3->1:40
    Abl(find(Abl<AbCut)) = 0;
    Abl(find(Abl>=AbCut)) = 1;

    Antibody = par.Antibody;
    OutbreakStartingDay = par.OutbreakStartingDay; %default 120 days 
    corrected_numdays = Antibody.numdays-OutbreakStartingDay;
    
%Transform observed titres Abl into 2D-array [time x titres]
for a=1:par.maxa
    obs_titres = zeros(par.SamplingLastDay+1, 2);
for i=1:length(Abl)
    titres = Abl(i);
    time = corrected_numdays(i);
    ind_age = age(i);
    if ind_age>=par.ages(a,1) & ind_age<par.ages(a,2)
        if time<par.SamplingLastDay+1
            obs_titres(time,titres+1) = obs_titres(time,titres+1)+1;
        end
    end
end
    observe(a).obs_titres = obs_titres;
    observe(a).obs_titres_numdays = OutbreakStartingDay+(0:par.SamplingLastDay)'; % from 1 Jan 2009
end
obs_titres = observe(1).obs_titres + observe(2).obs_titres + observe(3).obs_titres + observe(4).obs_titres;
numdays = Antibody.numdays;
    par.Antibody.Abl = Abl;  
    
    PriorMeta = Par_stat.ode.mcmc.Prior;
    x_name = {PriorMeta.name};
    %% Pre-existing Titres 
    %setup initial conditions
    if par.maxi == 2 % only 2 titres
        [y0 age_arr s0_imm] = make_ics_naive2titres( par, par.arrSlu, par.arrIlu, par.arrCIlu, Ab.age);
    else             % full titres
        [y0 age_arr s0_imm] = make_ics_naive( par, par.arrSlu, par.arrIlu, par.arrCIlu, Ab.age);
    end
    par = setParameters(par,'age_arr',age_arr);       
    par = setParameters(par,'s0_imm',s0_imm);

    
    %if length(thetabar) == 10
        %LogLikelihood = inline('-getNegLLHAgejava(pp, par, meser, y0)','pp','par','meser','y0');
    %    LogLikelihood = inline('-getNegLLHbyArrayjava_DIC(pp, x_name, par, meser, y0, observe)','pp','x_name','par','meser','y0', 'observe');
    %end
    if length(thetabar) >= 1
        LogLikelihood = inline('-getNegLLHbyArrayjava_DIC(pp, x_name, par, meser, y0, observe)','pp','x_name','par','meser','y0', 'observe');
        %LogLikelihood = inline('-getNegLLHjava(pp, par, meser, y0)','pp','par','meser','y0');
    end
    
    javaaddpath(par.javapath);
    import matlabjava.*
    mepar = matlabjava.Parameters;
    meser = matlabjava.Serology;
    meser.setParameters(mepar);
    meser.updateParameters('s0_imm',par.s0_imm);
    meser.updateParameters('wan',par.wan);
    meser.updateParameters('maxi', par.maxi); %set the model to be threshold
    

    samplesize = 80;
    burnIn = 1000;
    total = height(PosteriorSamples(:,1))-burnIn;
    idx = burnIn + round(rand(1, samplesize) * total);

 for i = 1:samplesize
    %    vars = PosteriorSamples.Properties.VariableNames;
    %    for p=1:length(vars)
    %        if strcmpi('LLH',vars(p))
    %        else
    %        [pars] = setParameters(par,char(vars(p)),PosteriorSamples(idx(i),p));
    %    end
    %    end
    pp = table2array(PosteriorSamples(idx(i),1:end-1));
    D_LLH_i = LogLikelihood(pp, x_name, par, meser, y0, observe);
    %[y0 age_arr] = make_ics_fromtitres_byage( par, par.arrSlu, par.arrIlu, par.arrCIlu, Ab(1).Abl, Ab(1).age); %use round1 data as initial condition
    %LLH_thetabar = LogLikelihood(thetabar, par, meser, y0); 
    D_i(i) = -2*D_LLH_i;
    %PD = Dbar - D_thetabar;
    %DIC = PD + Dbar;
 end
    theta_mean = mean(table2array(PosteriorSamples(:,1:end-1)));
    D_thetabar = -2*LogLikelihood(theta_mean, x_name, par, meser, y0, observe);
    PD1 = mean(D_i) - D_thetabar;
    PD2 = 0.5*var(D_i);
    DIC = D_thetabar + 2*PD1; 
    
    save([outdir '/dicvalue.mat'], 'DIC', 'D_i', 'D_thetabar');
end

