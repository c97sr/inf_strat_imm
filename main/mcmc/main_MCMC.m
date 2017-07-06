function [ output_args ] = main_MCMC( m, outdir, nsteps, sampleno, burnIn, savefigflag, targetdate)
%main_MCMC The main script parameter estimation of the models using MCMC
% MCMC metropolis-hasting algorithm to estimate parameters
% The algorithm accepts or rejects the proposed state based on the density
% of the target distribution
% main_MCMC(m, outdir, nsteps, sampleno)
% input parameters: model m, output directory outdir, number of MCMC steps nsteps, number of the samples sampleno   
% example1: run model 2.0 by steps=1000, sampleno=10 
%   main_MCMC(2.0, 'out/m2_0', 1000, 10)
%   main_MCMC(1.0, 'out/m1_0', 1000, 50, 100, 1, '')
% 6 Aug, 2014
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
par.strain = params.strain;

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

%pa.startdate = '01/01/2001' was defined during antibody extraction
%main_producetitres_bytime_mrc.m
no_season = length(Ab.K);
prompt = ['virus ' par.strain ' is selected. \n' ... 
          'Enter which season[2-' num2str(no_season) '] the serological profile you want to view: \n'];          
          %'Samples were recrutied from years ' num2str(str2num(start_year)) ' to ' num2str(str2num(start_year)+no_season-1) '\n'...
          

in_x = input(prompt);
second_collect = in_x; %strain we want to extract
init_collect = second_collect - 1;

par.initK = init_collect;
par.targetK = second_collect;
        
token = regexp(params.startdate, '/');
startyear = params.startdate(token(end)+1:end); %retrieve to start year from Ab data
par.startdate = params.startdate; %set the startdate
%when par.initK == 2, par.firstcase = '01/09/2002';
%define that first case starts on 44 weeks
%https://www.gov.uk/government/uploads/system/uploads/attachment_data/file/526405/Flu_Annual_Report_2015_2016.pdf
par.firstcase = ['01/10/' num2str(str2num(startyear) + par.initK-1)];
par.OutbreakStartingDay = datenum(par.firstcase,'dd/mm/yyyy') - datenum(par.startdate,'dd/mm/yyyy');

OutbreakStartingDay = par.OutbreakStartingDay;
Antibody = par.Antibody;
SampleSize = Antibody.samplesize;
Abl = Antibody.Abl;
age = Antibody.age;
corrected_numdays = Antibody.numdays-OutbreakStartingDay;

%% Re-formatting titres
Abl(find(Abl>9)) = 9; %substitute Ab level >9
%Only 0 and 1 immune status
AbCut = 3; %3 for 1:40, 2 for 1:20
if par.maxi == 2
    Abl(find(Abl<AbCut)) = par.mint; % 
    Abl(find(Abl>=AbCut)) = par.maxt;
end

%% Pre-existing titres 
%setup initial conditions
%take the titres from the initial season or from corrected_numdays -200 to 15 to generate pre-existing titres
if par.maxi == 2 % only 2 titres
    [y0 age_arr] = make_ics_naive2titres( par, par.arrSlu, par.arrIlu, par.arrCIlu, Ab.age);
else             % full titres
    [y0 age_arr] = make_ics( par);
end
par = setParameters(par,'age_arr',age_arr);       

%% Observed titres
%transform observed titres Abl into 2D-array [time x titres]
for a=1:par.maxa
    obs_titres = zeros(par.SamplingLastDay+1, par.maxi);
for i=1:length(Abl)
    titres = Abl(i);
    time = corrected_numdays(i);
    if time < 90 % only extract serosurveillance in next season
      continue;
    end
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

%% CREATE MCMC OUTPUT
%Setup the folder directory
if exist('outdir')
    mainoutdir = [pwd '/' outdir];
else
    mainoutdir = 'out/mcmc';
end
mainproj = params.proj;
outfile = ['mcmc_output_m' num2str(m)];
modelparfile =  ['parameters_m' num2str(m) '.mat'];
if ~exist('targetdate', 'var')
    out_dir_full = set_projectoutput( mainoutdir, [mainproj '/m' num2str(m)] );
    out_dir_rel = set_projectoutput( outdir, [mainproj '/m' num2str(m)]); 
else
    out_dir_full = set_projectoutput( mainoutdir, [mainproj '/m' num2str(m)], targetdate);
    out_dir_rel = set_projectoutput( outdir, [mainproj '/m' num2str(m)], targetdate); 
end

out_dir_rel

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

outfilename = [out_dir_full outfile '.mat'];
nofile = 0;
if exist(outfilename, 'file') == 2
   nofile = nofile + 1;
   newfileid = 2;
   for i=2:10
       outfilename = [out_dir_full outfile '(' num2str(i) ').mat'];
       if exist(outfilename) == 2
           newfileid = i+1;
       else
           break;
       end
   end
   outfilename = [out_dir_full outfile '(' num2str(newfileid) ').mat'];
end
outfilename

par.sample_time_K = round(mean(par.Antibody.K(par.targetK).numdays));
par.out_dir = out_dir_full;
par.out_dir_rel = out_dir_rel;

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
    PosteriorSamples = getOutput( out_dir_full, outfile, burnIn); % also remove MCMC results before burnIn
    save([out_dir_full outfile '_final.mat'] ,'PosteriorSamples','Ab','par','sys_par','Par_stat','elapsed');
    save_figures(PosteriorSamples,par,burnIn,sampleno, out_dir_full, outfile);
end
end





