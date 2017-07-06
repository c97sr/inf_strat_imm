function [ total_infecteds serodist_mat log CInfectedsT1T2 total_age_pop CInfectedsT2 CInfectedsT1 CInfecteds recovered_tot] = plot_immune( y, pars, times, strain, agegroup, display, marker )
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
recovered_tot = recovered_tot - recovered_tot(1); % extract the baseline

set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []); % No Xtick label
%xlabel('Days');

%plot(recovered_tot,'g');
%legend('recovered');
strain = 1;
[total_infecteds CInfecteds]= gen_total_infecteds(y, times, pars, strain, agegroup);
%plot(total_infecteds,'r');
%total_age_pop = susceptible_tot(1) + recovered_tot(1); % don't use it
%because recovered_tot already extract baseline
total_age_pop = 1; % 
log = ['CI(group' num2str(agegroup) '):' num2str(CInfecteds(end)/total_age_pop)];

formatIn = 'dd/mm/yyyy';
SDate = '01/05/2009';
numday_T1 = datenum('11/08/2009',formatIn)-datenum(SDate,formatIn)+1;
numday_T2 = datenum('22/12/2009',formatIn)-datenum(SDate,formatIn)+1;
CInfectedsT2 = CInfecteds(numday_T2);
CInfectedsT1 = CInfecteds(numday_T1);
CInfectedsT1T2 = CInfectedsT2  -  CInfectedsT1;
T = 1:length(times);


  I = getI();
  days = 1:7:43*7+1;
  x = days;
  y1 = I;
  yfull = zeros(1,length(T));
  for i1=1:T(end)
    idx = find(days==i1);
    if idx > 0 
       yfull(i1) = y1(idx);
    else
        if i1>1
           yfull(i1) = yfull(i1-1); 
        end
    end
  end
yfull = yfull./(max(yfull)*125);  
yfull = zeros(1,length(T'))-100;
if display == 0
else
%[AX,H1,H2] = plotyy(T',[recovered_tot/total_age_pop],T',yfull', 'plot','bar');
%[AX,H1,H2] = plotyy(T',[recovered_tot/total_age_pop],-1,-1, 'plot','bar');

%set(AX,{'ycolor'},{'b';'r'})  % Left color green, right color red...

    ylim('manual');
    %set(AX(1),'ylim',[0:0.1:1]);
if exist('marker')
    %set(H1,'Color','b');
    %set(H2, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', [0.8 0.8 0.8]);
    %uistack(AX(1));
    %set(AX(1), 'Color', 'none');
    %set(AX(2), 'Color', 'w');
    %set(AX(2),'ylim',[0,0.02]);%
    %set(AX(2),'YTick',[0:0.004:0.02]);
    %set(AX(2),'YTickLabel',[0:0.4:2]);
    %set(AX(2),'ylim',[0,0.02]);
    
    [AX,H1,H2] = plotyy(T',[recovered_tot/1],-1,-1, 'plot','bar');
    set(AX,{'ycolor'},{'b';'b'})  % Left color green, right color red...
    set(H1,'Color','b');
    set(H2, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', [0.8 0.8 0.8]);
    set(AX(1),'ylim',[0,0.4]);%
    set(AX(1),'ylim',[0,0.4]);%
    set(AX(1),'YTick',[0:0.1:0.4]);
    set(AX(1),'YTickLabel',[0:0.1:0.4]);
    set(AX(2),'YTickLabel','');
    set(H1,'Color','b');
    set(H1, 'linestyle', '--');
    set(H1,'linewidth',0.5) 
    set(H2, 'linestyle', '--');
    
else
    [AX,H1,H2] = plotyy(T',[recovered_tot/1],-1,-1, 'plot','bar');
    set(AX,{'ycolor'},{'b';'b'})  % Left color green, right color red...
    set(H1,'Color','b');
    set(H2, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', [0.8 0.8 0.8]);
    uistack(AX(1));
    set(AX(1), 'Color', 'none');
    set(AX(2), 'Color', 'w');
    set(AX(1),'ylim',[0,0.4]);%
    set(AX(2),'ylim',[0,0.4]);%
    %set(AX(2),'ylim',[0,0.02]);%
    %set(AX(2),'YTick',[0:0.004:0.02]);
    set(AX(1),'YTickLabel',[0:0.1:0.4]);
    set(AX(2),'YTickLabel','');
    set(get(AX(1),'Ylabel'),'String','Seroconversion, cumulative incidence (%)', 'FontSize', 12)  
    %set(get(AX(2),'Ylabel'),'String','Infected (%)', 'FontSize', 12) 
end
    set(AX(1),'YTick',[0:0.1:1]);
    set(AX(1),'YTickLabel',[0:10:100]);

set(AX(1),'xlim',[0,30.5*10]);
set(AX(2),'xlim',[0,30.5*10]);
set(AX(1),'XTickLabel','');
set(AX(2),'XTickLabel','');
set(gca,'xtick',[]);
disp(['Seroprevalence:' num2str(recovered_tot(end)-recovered_tot(1))]);
disp done;
end
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

function [ infecteds CInfecteds] = gen_total_infecteds__( y, times, pars, strain, agegroup)
    Slu = pars.arrSlu(agegroup,:,:,:);
    Ilu = pars.arrIlu(strain,agegroup,:,:,:);
    CIlu = pars.arrCIlu(strain,agegroup,:,:,:);
    dims_Slu = prod(size(pars.arrSlu(strain,:,:,:))); %[a i j x]
    dims_Ilu = prod(size(pars.arrIlu(strain,:,:,:,:))); %[X a i j k]
    
    %infecteds = sum(y(:,dims_Slu+1:dims_Slu+dims_Ilu),2); % total infected number, regardless which strains
    infecteds = sum(y(:,Ilu),2); % total infected number, regardless which strains
    CInfecteds = sum(y(:,CIlu),2);
end

