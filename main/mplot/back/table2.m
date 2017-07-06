function [ ] = plot_dynamics_bypar( pars )
% Summary of the function table2
% Plot disease and immune dynamics
% For each age group, 5 subfigures are plotted in a row. 
% subplot a) disease dynamics, b) heat map of immune dynamics, c) histogram of antibody levels
% d) HK sera(final) after outbreak and e) d) HK sera(final) after outbreak
% Written by Sean Yuan (hyuan@imperial.ac.uk) 
% Plot all age group information.

p = path;
%path(p,'../');
path(p,'lib/');
setISL;

global proj Antibody;
init_collect = 1;
second_collect = 2;
third_collect = 3;
k = 2;

%retrieve parameters
beta = pars.beta;
AbB = [pars.AbB1 pars.AbB2 pars.AbB3 pars.AbB4];
immune_alpha = [pars.immune_alpha1 pars.immune_alpha2 pars.immune_alpha3 pars.immune_alpha4];
lastsamplingday = pars.SamplingLastDay;

%setup initial condition
[yini age_arr] = make_ics_naive( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu);
[yini_k1 age_arr_k1] = make_ics_fromtitres_byage( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, Antibody.K(init_collect).Abl, Antibody.K(init_collect).age);
[yini_k2 age_arr_k2] = make_ics_fromtitres_byage( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, Antibody.K(k).Abl, Antibody.K(k).age);
[yini_k3 age_arr_k3] = make_ics_fromtitres_byage( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, Antibody.K(k).Abl, Antibody.K(k).age);


%setep simulation time
T0 = pars.OutbreakStartingDay;
meanKdays(1) = mean(pars.Antibody.K(1).numdays - T0);
meanKdays(2) = mean(pars.Antibody.K(2).numdays - T0);
meanKdays(3) = mean(pars.Antibody.K(3).numdays - T0);
sample_time_K1 = round(meanKdays(1));
sample_time_K2 = round(meanKdays(2));
times = 0:1:lastsamplingday;

%run simulation
%initialize objects
javaaddpath e:\workspace\MyJavaProject\bin\matlabjava.jar;
mepar = matlabjava.Parameters
meser = matlabjava.Serology
meser.setParameters(mepar);
% set parameters
meser.updateParametersG(pars.arrg);
meser.updateParametersH(pars.arrh);
meser.updateParametersM(pars.matM);
meser.updateParametersBeta(pars.beta);  
meser.updateParameters('wan',pars.wan);
x0 = yini;  
[t y] = ode23(@(t,x)odef_islmodjava(t,x, meser), times, x0);  
%[t y] = ode23(@(t,x)odef_islmod(t,x,pars), times, x0);

 
%plot disease dynamics
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1500, 360]);

for a=1:4
strain = 1;


[ titremat] = gen_strain_titres( yini, 1, a, pars );

% Plot 3rd column
% susceptible population dynamics from model output
Xout_k1 = retrieve_histogram(y, pars, times, sample_time_K1, a); % model output for K1
Xout_k2 = retrieve_histogram(y, pars, times, sample_time_K2, a); % model output for K2

% susceptible population dynamics from HK data
Yout_k1 = retrieve_histogram(yini_k1, pars, 1, 1, a);
Yout_k2 = retrieve_histogram(yini_k2, pars, 1, 1, a);

Yout3 = [Yout_k1' Yout_k2' Xout_k2'];
Yout4 = [Yout_k1' Xout_k1' Yout_k2' Xout_k2'];


titres = 0:1:pars.maxi; % x axis
count_titres = Yout4(1:pars.maxi,:); % y axis
count_titres(1,:) = count_titres(1,:)/5; % normalize titres level 0
for i=pars.maxi:-1:2
  count_titres(i+1,:) = count_titres(i,:); 
end
count_titres(2,:) = 0;

[hax,hline,hbar] = plotyy([1],[0],titres,count_titres,'plot','bar');
set(hax(1),'YLim',[0 1.001]);
set(hax(1),'YTick',[0:0.5:1]);
set(hax(1),'YTickLabel',[0:0.5:1]);
set(hax(2),'YLim',[0 0.2]);
set(hax(1),'XLim',[-0.5 10.5]);
set(hax(1),'XTickLabel','');
set(hax(2),'XLim',[-0.5 10.5]);
set(hax(2),'XTick',[0:1:10]);
set(hax(2),'XTickLabel',['0';' ';'1';'2';'3';'4';'5';'6';'7';'8';'9']);
title(['HK sera(T' num2str(k) ')']);
end %--end of age loop

end