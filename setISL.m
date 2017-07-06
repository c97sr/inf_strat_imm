restoredefaultpath;

% Set path
p = path;
path(p,[pwd '/lib/']);
p = path;
path(p,[pwd '/lib/model']);
p = path;
path(p,[pwd '/lib/model/llh']);
p = path;
path(p,[pwd '/lib/model/susc']);
p = path;
path(p,[pwd '/lib/model/rt']);
p = path;
path(p,[pwd '/lib/chart']);
p = path;
path(p,[pwd '/lib/optim']);
p = path;
path(p,[pwd '/lib/sys']);

p = path;
path(p,[pwd '/main/likelihood_estimation']);
p = path;
path(p,[pwd '/main/calculate_R0']);
p = path;
path(p,[pwd '/main/extract_antibody_titres']);
p = path;
path(p,[pwd '/main/mcmc']);
p = path;
path(p,[pwd '/main/mcmc/mode2']);
p = path;
path(p,[pwd '/main/mcmc/mode4']);
p = path;
path(p,[pwd '/main/plot']);
p = path;
path(p,[pwd '/main/mplot']);
p = path;
path(p,[pwd '/main/mplot/back']);
p = path;
path(p,[pwd '/main/produce_samples']);

% Declare global variables
global params proj dat Antibody
%proj = 'dat/20141206/hk_ph1n1/'; %exlude vaccinated individuals; paired sera
%load([proj 'h1n1_titres.mat']);
%proj = 'dat/20161114/world_h3n2/';
proj = 'dat/20161117/world_h3n2/';
load([proj 'h3n2_titres_2.mat']);  %Syn1997
%load([proj 'h3n2_titres_3.mat']);   %Fuj2002
%proj = 'dat/20161011/world_h1n1/';
%load([proj 'h1n1_titres.mat']);


