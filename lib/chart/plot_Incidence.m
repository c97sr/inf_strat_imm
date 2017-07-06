function [ total_infecteds serodist_mat log CInfectedsT1T2 total_age_pop rmse_p rmse_f] = plot_Incidence( y, pars, times, strain, agegroup, marker )
%plot_line Summary of this function goes here
% plot the total susceptible individuals along the time
% Written by Sean Yuan (hyuan@imperial.ac.uk)
% Modify the code on Jun 16, 2014
% To return an infecteds vector including each age groups

arrSlu = pars.arrSlu;
arrIlu = pars.arrIlu;

%Setup some standard auxilliary functions
tmin = min(times);
tmax = max(times);

% Antibody level
% Susceptible = Antibody <1:40
% Recovered = Antibody >=1:40
if pars.maxi == 10
    abl = 4; %1:40
end
if pars.maxi == 2
    abl = 2;
end

% Define age groups
arrIlu = pars.arrIlu; %I lookup
arrSlu = pars.arrSlu; %S lookup
serodist = zeros(length(times),pars.maxi); % [times x titres(9)] 
serodist_mat = []; % [age x times x titre]
for a = 1:pars.maxa
    serodist = gen_strain_titres(y, times, a, pars); % [times x titre]
    serodist_mat(a,:,:) = reshape(serodist, [1, size(serodist)]);  % [age x times x titre]
end

% Total
if ~exist('agegroup')
    agegroup = 1:pars.maxa;
end
serodist_tot3d = sum(serodist_mat(agegroup,:,:),1);  % [1 x times x titre]
serodist_tot = reshape(serodist_tot3d, size(serodist_tot3d,2), []); % [times x titre]
susceptible_tot = sum(serodist_tot(:,1:abl-1),2); % [times x 1], sum of individuals for titres <1:40; check abl=4
recovered_tot = sum(serodist_tot(:,abl:pars.maxi),2); % [times x 1], sum of individuals for titres >=1:40
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []); % No Xtick label

strain = 1;
[total_infecteds CInfecteds]= gen_total_infecteds(y, times, pars, strain, agegroup);
total_age_pop = susceptible_tot(1) + recovered_tot(1);
log = ['CI(group' num2str(agegroup) '):' num2str(CInfecteds(end)/total_age_pop)];

formatIn = 'dd/mm/yyyy';
SDate = '01/05/2009';
numday_T1 = datenum('11/08/2009',formatIn)-datenum(SDate,formatIn)+1;
numday_T2 = datenum('22/12/2009',formatIn)-datenum(SDate,formatIn)+1;
CInfectedsT1T2 = CInfecteds(numday_T2) -  CInfecteds(numday_T1);

T = 1:length(times);
It = getI();
Tobs = 1:7:7*(length(It)-1)+1;
Xt = total_infecteds(Tobs);
Yt = It;
rep = 1; %reporting rate
x = rep;
N = pars.N;
[rmse_p rmse_f] = calRMSE(rep, Xt, Yt, pars.N);
%x_hat = fminsearch(@(x)calLL(x, Xt, Yt));
options = optimset('Display','iter');
%x_hat = fminsearch(@calLL,x,options,Xt,Yt,N);
%x_hat
%fun = @(x) ;

    %[AX,H1,H2] = plotyy([ 0],[ 0],T',total_infecteds, 'plot'); 
    %set(H1(1),'Color','b');
    %set(H2,'Color','r');
    %set(AX(2),'ylim',[0,0.02]);%
    %set(AX(2),'YTick',[0:0.004:0.02]);
    %set(AX(2),'YTickLabel',[0:0.4:2]);
    %set(get(AX(1),'Ylabel'),'String','Susceptible, immune protected (%)', 'FontSize', 12)  
    %set(get(AX(2),'Ylabel'),'String','Infected (%)', 'FontSize', 12)  
        
    %set(AX(1),'YTick',[0:0.2:1]);
    %set(AX(1),'YTickLabel',[0:20:100]);
    %set(AX(1),'xlim',[0,30.5*10]);
    %set(AX(2),'xlim',[0,30.5*10]);
    %set(AX(1),'XTickLabel','');
    %set(AX(2),'XTickLabel','');
    %set(gca,'xtick',[]);

end

function i = getI()
i = [
2
0
1
9
47
23
72
275
384
239
343
558
930
1261
1344
1082
1630
1760
1659
2400
2914
3326
1660
1013
469
257
234
123
115
142
296
319
260
249
237
258
210
214
202
150
113
98
81
64
];
end




