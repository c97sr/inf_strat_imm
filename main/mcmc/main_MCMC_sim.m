function [ output_args ] = main_MCMC_sim( m, outdir, nsteps, sampleno, burnIn)
%main_MCMC The main script parameter estimation of the models using MCMC
% MCMC metropolis-hasting algorithm to estimate parameters
% The algorithm accepts or rejects the proposed state based on the density
% of the target distribution
% main_MCMC(m, outdir, nsteps, sampleno)
% input parameters: model m, output directory outdir, number of MCMC steps nsteps, number of the samples sampleno   
% example1: run model 2.0 by steps=1000, sampleno=10 
%   main_MCMC(2.0, 'out/m2_0', 1000, 10)
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
        if time<par.SamplingLastDay+1
            obs_titres(time,titres+1) = obs_titres(time,titres+1)+1;
        end
    end
end
    observe(a).obs_titres = obs_titres;
    observe(a).obs_titres_numdays = OutbreakStartingDay+(0:par.SamplingLastDay)'; % from 1 Jan 2009
end
obs_titres = observe(1).obs_titres + observe(2).obs_titres + observe(3).obs_titres + observe(4).obs_titres;
load('sim_sample_set4106.mat');

%% Pre-existing Titres 
%setup initial conditions
if par.maxi == 2 % only 2 titres
    [y0 age_arr s0_imm] = make_ics_naive2titres( par, par.arrSlu, par.arrIlu, par.arrCIlu, Ab.age);
else             % full titres
    [y0 age_arr s0_imm] = make_ics_naive( par, par.arrSlu, par.arrIlu, par.arrCIlu, Ab.age);
end
par = setParameters(par,'age_arr',age_arr);       
par = setParameters(par,'s0_imm',s0_imm); 

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
[PosteriorSamples PriorPoint] = runMCMC(meser, par, y0, sim_samples, Par_stat); % change observe to sim_samples
elapsed = toc;
clear('mepar');
clear('meser');

%% CREATE MCMC OUTPUT
%Setup the folder directory
if exist('outdir')
    mainoutdir = [pwd '/' outdir];
else
    mainoutdir = 'out/mcmc';
end
mainproj = 'ph1n1';
outfile = ['mcmc_output_m' num2str(m) '.mat'];
[out_dir ] = set_projectoutput( mainoutdir, mainproj)

%% Save the Parameters and Posterior into m file
disp('save output files');
sys_par.nSamples = sampleno;
sys_par.burnIn = burnIn;
sys_par.PriorPoint = PriorPoint;
sys_par.PriorMeta = Par_stat.ode.mcmc.Prior;
Par_stat.maxlikelihood = max(PosteriorSamples.LLH); %Add maximum likelihood
save([out_dir outfile] ,'PosteriorSamples','Ab','par','sys_par','Par_stat','elapsed');
csvwrite([out_dir 'posterior.csv'],table2array(PosteriorSamples));

%% Save the figure and close the figure windows
setISL
disp('save posterior');
[FigH] = mcmc_posterior_hist(PosteriorSamples);
print('-dpng','-r0',[out_dir 'mcmc_posterior_histogram_m' num2str(par.model) '.png'])
savefig(FigH,[out_dir 'mcmc_posterior_histogram_m' num2str(par.model) '.fig']);
close(FigH);
if (exist('sampleno'))
    SampleNo = sampleno;
else
    SampleNo = 20;    
end
BurnIn = sys_par.burnIn;
setISL;
disp('save figure2');
%[FigH FigL log] = figure2(PosteriorSamples,par,BurnIn,SampleNo); %use 400 for final
[FigH FigL log peak peak_lb peak_ub] = figure2(PosteriorSamples,par,BurnIn,SampleNo); %use 400 for final
savefig(FigH,[out_dir 'fig2a.fig']);
close(FigH);
savefig(FigL,[out_dir 'fig2b.fig']);
close(FigL);

output.out_dir = out_dir;
output.peak = peak;
output.peak_lb = peak_lb;
output.peak_ub = peak_ub;

save([out_dir outfile] ,'PosteriorSamples','Ab','par','sys_par','Par_stat','elapsed','output');

fileID = fopen([out_dir 'mylogs.txt'],'w');
str = sprintf(log)
fprintf(fileID, '%s', str);
fclose(fileID);

if par.maxi > 2
setISL;
disp('save figure3');
FigH = figure3(PosteriorSamples,par,BurnIn,SampleNo);
savefig(FigH,[out_dir 'fig3.fig']);
close(FigH);


setISL;
FigH = figure4( PosteriorSamples,par,BurnIn,SampleNo);
savefig(FigH,[out_dir 'fig4.fig']);
close(FigH);

setISL;
[S1 T1]= table1(PosteriorSamples, BurnIn);
writetable(T1,[out_dir 'Table1.csv'],'WriteRowNames',false);

setISL;
[Table2Left Table2Right p] = table2( PosteriorSamples, par,BurnIn,SampleNo);
writetable(Table2Left,[out_dir 'Table2a.csv']);
writetable(Table2Right,[out_dir 'Table2b.csv']);
end

setISL;
disp('save figure5');
[FigH Rt] = figure5(PosteriorSamples,par,BurnIn,SampleNo);
savefig(FigH,[out_dir 'fig5.fig']);
close(FigH);
save([out_dir 'Rt_m' num2str(par.model) '.mat'], 'Rt');
setISL; 

csvwrite(['temp/posterior.csv'],table2array(PosteriorSamples));
Tess = calEffectiveSampleSize;
writetable(Tess,[out_dir 'TableESS.csv'],'WriteRowNames',false);
end



