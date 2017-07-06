function [ x_opt negLL ] = main_findmin_sera_age_java( model )
%MAIN_FINDMIN Summary of this function goes here
%   Detailed explanation goes here
%no age subgroup
%http://www.mathworks.co.uk/help/optim/ug/fmincon.html?s_tid=doc_12b
%http://www.mathworks.co.uk/help/matlab/ref/function_handle.html
%fminbnd(@sqr,-10,10)
% 19 August 2014

%%
%Assign global variables to local variables; When I use java file, global
%variables will be lost after java is running.
setISL;
Ab = Antibody; 
Pr = proj;

%%
par = InitParameters(); 
par.Antibody = Ab;
%par = setParameters(par,'age_flag',1);
%par = setParameters(par,'model',model);
%par = setParameters(par,'frac_flag',1);
%par = setParameters(par,'ContFrac1', 4); %contfrac1 = 4 still won't produce good immune_alpha
%par = setParameters(par,'PUAb', 0.25); %default pa.PUAb = 0.06; 

    [par Par_stat] = initialModel(model, par);
    A = Par_stat.ode.opt.A;
    b = Par_stat.ode.opt.b;
    x0 = Par_stat.ode.opt.x0;
    x_name = Par_stat.ode.opt.x_name;
    lb = Par_stat.ode.opt.lb;
    ub = Par_stat.ode.opt.ub;

%setup initial conditions from titres among all age groups 
abl_ini= Antibody(1).Abl(find(Antibody(1).age>=par.ages(1,1) & Antibody(1).age<par.ages(par.maxa,2))); %all age groups
%[yini age_arr] = make_ics_fromtitres_byage( par, par.arrSlu, par.arrIlu, par.arrCIlu, Antibody(1).Abl, Antibody(1).age); %use round1 data as initial condition
[yini age_arr s0_imm initS] = make_ics_naive( par, par.arrSlu, par.arrIlu, par.arrCIlu, Antibody.age);
par = setParameters(par,'initS',initS);
par = setParameters(par,'age_arr',age_arr);
par = setParameters(par,'s0_imm',s0_imm);
%A = 1;
%b = 5;
%x0 = 0.3;
options = optimset('fmincon');
options = optimset(options,'MaxFunEvals', 200000, 'Algorithm', 'active-set'); %Defualt = 100 * number of variables
options.MaxIter = 10^(13);
options.ObjectiveLimit = -1e25;

%lb = [0.2; 2; -0.01];
%ub = [1; 8; 0.999];

%initialize objects
javaaddpath e:\workspace\MyJavaProject\bin\matlabjava.jar;
mepar = matlabjava.Parameters
meser = matlabjava.Serology
meser.setParameters(mepar);
meser.updateParameters('wan',par.wan);
meser.updateParameters('s0_imm',par.s0_imm);


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


%optimization
%[x_opt negLL] = fmincon(@(x)getNegLLHAgejava(x,x_name,par,meser,yini),x0,A,b,[],[],lb,ub,[],options);
[x_opt negLL] = fmincon(@(x)getNegLLHbyArrayjava(x,x_name,par,meser,yini,observe),x0,A,b,[],[],lb,ub,[],options);

% Update parameters
for i=1:length(x_opt)
    par = setParameters(par,x_name{i},x_opt(i));
end



%Beta = x_opt(1);
%[Rt_list beta_list] = calR0_bypar(Beta,par); %R0 at what time point?
%Rt = Rt_list(1);
%par.Rt = Rt;

%save output
mainoutdir = 'out/opt';
mainproj = 'ph1n1';
outfile = ['opt_maxllh_m' num2str(model) '.mat'];
[ out_dir ] = set_projectoutput( mainoutdir, mainproj)
save([out_dir outfile] ,'x_opt','negLL','par','Par_stat');

setISL
plot_dynamics_bypar(par)
print('-dpng','-r0',[out_dir 'dynamics_m' num2str(par.model) '.png'])
plot_infected_distribution(par)
print('-dpng','-r0',[out_dir 'secondinfection_m' num2str(par.model) '.png'])


%save('mle_10pars.mat','x_opt');

end

