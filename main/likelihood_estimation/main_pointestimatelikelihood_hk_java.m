function [ nLLH ] = main_pointestimatelikelihood_hk_java( par, x_opt )
% main_pointestimatelikelihood_hk_java
% To estimate likelihood for obserbed antibody titres given parameters par
% Detailed explanation
% Run setISL to setup the running environment
% Call java jre model from matlab
% x_opt is the table of the parameters
% Hsiang-Yu Yuan
% Check output [proj]/herdimmunity_hk_ph1n1/likelihood_2params.mat
% 2 Dec 2014

%initialize global variables
setISL;
Ab = Antibody; 
Pr = proj;

%setup initial conditions
abl_ini= Antibody(1).Abl(find(Antibody(1).age>=par.ages(1,1) & Antibody(1).age<par.ages(par.maxa,2))); %all age groups
[yini age_arr s0_imm] = make_ics_naive( par, par.arrSlu, par.arrIlu, par.arrCIlu, Antibody.age);

%create observerd data object
OutbreakStartingDay = par.OutbreakStartingDay; %default 120 days 
Antibody = par.Antibody;
SampleSize = Antibody.samplesize;
Abl = Antibody.Abl;
age = Antibody.age;
numdays = Antibody.numdays-OutbreakStartingDay;
Abl(find(Abl>par.maxt)) = par.maxt; %substitute Ab level >maxt

%transform observed titres into 2D-array [time x titres]
for a=1:par.maxa
    obs_titres = zeros(par.SamplingLastDay+1, par.maxi);
for i=1:length(Abl)
    titres = Abl(i);
    time = numdays(i);
    ind_age = age(i);
    if ind_age>=par.ages(a,1) & ind_age<par.ages(a,2)
        if time<par.SamplingLastDay+1
            obs_titres(time,titres+1) = obs_titres(time,titres+1)+1;
        end
    end
end
    observe(a).obs_titres = obs_titres;
end


% estimate likelihood for each sampling time

%initialize java objects
javaaddpath e:\workspace\MyJavaProject\bin\matlabjava.jar;
mepar = matlabjava.Parameters
meser = matlabjava.Serology
meser.setParameters(mepar);

%theta = [0.0583,5.7553,4.8019,3.6589,3.1168,1.7794,30]
%theta_name =  {'beta','AbB1','AbB2','AbB3','AbB4','immune_alpha','immune_beta'};
theta = table2array(x_opt);
theta(end) = [];
theta_name = x_opt.Properties.VariableNames;

nLLH = getNegLLHbyArrayjava(theta,theta_name,par,meser,yini,observe)





