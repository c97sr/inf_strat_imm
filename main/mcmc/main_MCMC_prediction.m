function [ output_args ] = main_MCMC_prediction( m, outdir, nsteps, sampleno, burnIn, savefigflag, targetdate)
%main_MCMC The main script parameter estimation of the models using MCMC
% MCMC metropolis-hasting algorithm to estimate parameters
% The algorithm accepts or rejects the proposed state based on the density
% of the target distribution
% main_MCMC(m, outdir, nsteps, sampleno)
% input parameters: model m, output directory outdir, number of MCMC steps nsteps, number of the samples sampleno   
% example1: run model 2.0 by steps=1000, sampleno=10 
%   main_MCMC(2.0, 'out/m2_0', 1000, 10)
% 24 Aug, 2015
% Hsiang-Yu Yuan

%% INITIALIZE THE MODEL AND ITS PARAMETERS
%ENVIRONMENTAL VARIABLES
%Assign global variables to local variables; Because global variables are not accessible when java is running.
setISL;
Ab = Antibody; 
Pr = proj;

%Define model parameters
par = InitParameters(); 
par.Antibody = Ab;

%Initialize Prior density function for the Model
par = setParameters(par,'model',m);
if exist('nsteps')
    [par Par_stat] = initialModel(m, par, nsteps);
else
    [par Par_stat] = initialModel(m, par);
end
disp(['Total number of iterations: ' num2str(nsteps)]);
disp(['Initial Parameters: ']);
disp(par);
%Create observerd data object Abl
OutbreakStartingDay = par.OutbreakStartingDay; %default 120 days 
Antibody = par.Antibody;
SampleSize = Antibody.samplesize;
Abl = Antibody.Abl;
age = Antibody.age;
corrected_numdays = Antibody.numdays-OutbreakStartingDay;
formatIn = 'dd/mm/yyyy';
sera_sample_end = datenum({'31/09/2009'},formatIn) - datenum({'01/01/2009'},formatIn)  - OutbreakStartingDay;
par = setParameters(par,'sera_sample_end',sera_sample_end);

Abl(find(Abl>9)) = 9; %substitute Ab level >9

%Only 0 and 1 immune status
if par.maxi == 2
    Abl(find(Abl<3)) = par.mint;
    Abl(find(Abl>=3)) = par.maxt;
end

%Transform observed titres Abl into 2D-array [time x titres]
for a=1:par.maxa
    obs_titres = zeros(par.SamplingLastDay+1, par.maxi);
for i=1:length(Abl)
    titres = Abl(i);
    time = corrected_numdays(i);
    ind_age = age(i);
    if ind_age>=par.ages(a,1) & ind_age<par.ages(a,2)
        %if time<par.SamplingLastDay+1
        if time<sera_sample_end+1 % only extract sera samples before sera_sample_end+1
            obs_titres(time,titres+1) = obs_titres(time,titres+1)+1;
        end
    end
end
    observe(a).obs_titres = obs_titres;
end
observe(a).obs_titres_numdays = OutbreakStartingDay+(0:par.SamplingLastDay)'; % from 1 Jan 2009
obs_titres = observe(1).obs_titres + observe(2).obs_titres + observe(3).obs_titres + observe(4).obs_titres;


%% Pre-existing Titres 
%setup initial conditions
if par.maxi == 2 % only 2 titres
    [y0 age_arr s0_imm] = make_ics_naive2titres( par, par.arrSlu, par.arrIlu, par.arrCIlu, Ab.age);
else             % full titres
    [y0 age_arr s0_imm] = make_ics_naive( par, par.arrSlu, par.arrIlu, par.arrCIlu, Ab.age);
end
par = setParameters(par,'age_arr',age_arr);       
par = setParameters(par,'s0_imm',s0_imm); 

%% CREATE MCMC OUTPUT
%Setup the folder directory
if exist('outdir')
    mainoutdir = [pwd '/' outdir];
else
    mainoutdir = 'out/mcmc';
end
mainproj = 'ph1n1';
outfile = ['mcmc_output_m' num2str(m)];
modelparfile =  ['parameters_m' num2str(m) '.mat'];
if ~exist('targetdate', 'var')
    [out_dir ] = set_projectoutput( mainoutdir, mainproj);
else
    [out_dir ] = set_projectoutput( mainoutdir, mainproj, targetdate);
end

%% CREATE JAVA MODEL
%javaaddpath 'e:\workspace\MyJavaProject\bin\matlabjava.jar';
javaaddpath(par.javapath);
import matlabjava.*
mepar = matlabjava.Parameters;
meser = matlabjava.Serology;
meser.setParameters(mepar);
meser.updateParameters('s0_imm',par.s0_imm);
meser.updateParameters('wan',par.wan);
meser.updateParameters('maxi', par.maxi);
% S0_imm and PUAB, are they same?

%% MCMC simulation
tic;
Par_stat.ode.mcmc.nSteps = nsteps;
Par_stat.ode.mcmc.burnIn = burnIn;
[PosteriorSamples PriorPoint] = runMCMC(meser, par, y0, observe, Par_stat);
elapsed = toc;
clear('mepar');
clear('meser');

%% Save the Parameters and Posterior into m file
sys_par.nSamples = sampleno;
sys_par.burnIn = burnIn;
sys_par.PriorPoint = PriorPoint;
sys_par.PriorMeta = Par_stat.ode.mcmc.Prior;
Par_stat.maxlikelihood = max(PosteriorSamples.LLH); %Add maximum likelihood

outfilename = [out_dir outfile '.mat'];
if exist(outfilename, 'file') == 2
   newfileid = 2;
   for i=2:10
       outfilename = [out_dir outfile '(' num2str(i) ').mat'];
       if exist(outfilename) == 2
           newfileid = i+1;
       else
           break;
       end
   end
   outfilename = [out_dir outfile '(' num2str(newfileid) ').mat'];
end
outfilename
if nsteps < 10
    disp('not save output files');
else
    save([outfilename] ,'PosteriorSamples','Ab','par','sys_par','Par_stat','elapsed');
    disp('save output files');
end

%% Get figures output
if ~exist('savefigflag', 'var')
  savefigflag = 0;
end
if savefigflag == 1
    PosteriorSamples = getOutput( out_dir, outfile, burnIn);
    save_figures(PosteriorSamples,par,burnIn,sampleno, out_dir, outfile);
    save([out_dir outfile '_final.mat'] ,'PosteriorSamples','Ab','par','sys_par','Par_stat','elapsed');
end

end



