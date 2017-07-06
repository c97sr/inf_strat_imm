function [ total_infecteds serodist_mat log CInfectedsT1T2] = plot_line( y, pars, times, strain, agegroup, marker )
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

%calculate cumulated incidence
%CI = sum(y(end,pars.maxi*pars.maxa*3+1:pars.maxi*pars.maxa*4));


%plot(times', susceptible_tot);
%legend('susceptible(<=1:40)');
%xlim([0 466]);
%hold;

set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []); % No Xtick label
%xlabel('Days');

%plot(recovered_tot,'g');
%legend('recovered');
strain = 1;
[total_infecteds CInfecteds]= gen_total_infecteds(y, times, pars, strain, agegroup);
%plot(total_infecteds,'r');
total_age_pop = susceptible_tot(1) + recovered_tot(1);
log = ['CI(group' num2str(agegroup) '):' num2str(CInfecteds(end)/total_age_pop)];

formatIn = 'dd/mm/yyyy';
SDate = '01/05/2009';
numday_T1 = datenum('11/08/2009',formatIn)-datenum(SDate,formatIn)+1;
numday_T2 = datenum('22/12/2009',formatIn)-datenum(SDate,formatIn)+1;
CInfectedsT1T2 = CInfecteds(numday_T2) -  CInfecteds(numday_T1);

%[AX,H1,H2] = plotyy(times',[susceptible_tot/total_age_pop recovered_tot/total_age_pop CInfecteds/total_age_pop],times',total_infecteds, 'plot');
T = 1:length(times);

%[AX,H1,H2] = plotyy(T',[susceptible_tot/total_age_pop CInfecteds/total_age_pop],T',total_infecteds, 'plot');
%use Cumulated Infecteds insteads of recovered individuals >= 1:40
[AX,H1,H2] = plotyy(T',[susceptible_tot/total_age_pop recovered_tot/total_age_pop],T',total_infecteds, 'plot');
set(AX,{'ycolor'},{'b';'r'})  % Left color green, right color red...
%set(AX(1),'YScale', 'linear', 'YTick',[0; 1]); 
%set(AX(2),'YScale', 'linear', 'YTick',[0; 5*10^-4]);
 
    ylim('manual');
    %set(AX(1),'ylim',[0:0.1:1]);
if exist('marker')
    set(AX(2),'ylim',[0,0.02]);
    set(AX(2),'YTick',[0:0.004:0.02]);
    set(AX(2),'YTickLabel','');
    set(H1(1),'Color','b');
    set(H1(2),'Color','b');
    set(H2,'Color','r');
    set(H1, 'linestyle', '--');
    set(H2, 'linestyle', '--');
else
    set(H1(1),'Color','b');
    set(H1(2),'Color','b');
    set(H1(2),'linewidth',2) 
    set(H2,'Color','r');
    set(AX(2),'ylim',[0,0.02]);%
    set(AX(2),'YTick',[0:0.004:0.02]);
    set(AX(2),'YTickLabel',[0:0.4:2]);
    %set(AX(2),'YTickLabel',[0:0.5:2]);
    %axes(AX(1)); ylabel('Susceptible, cumulated incidences (%)');
    set(get(AX(1),'Ylabel'),'String','Susceptible, cumulative incidence (%)', 'FontSize', 12)  
    set(get(AX(2),'Ylabel'),'String','Infected (%)', 'FontSize', 12)  
    %ylabel('Susceptible, cumulated incidences (%)');
end
    set(AX(1),'YTick',[0:0.2:1]);
    set(AX(1),'YTickLabel',[0:20:100]);

set(AX(1),'xlim',[0,30.5*10]);
set(AX(2),'xlim',[0,30.5*10]);
set(AX(1),'XTickLabel','');
set(AX(2),'XTickLabel','');
set(gca,'xtick',[]);
disp(['Seroprevalence:' num2str(recovered_tot(end)-recovered_tot(1))]);
disp done;
end

function [ titremat] = gen_strain_titres__( y, times, agestates, pars, strain)
  nots = times; %times
  nostates = pars.maxi;
  titremat = zeros(length(nots), nostates);
  X = strain;
  maxl = 1;
  maxm = 1;
  maxn = 1;  
  if X == 1     %X = 1; %only plot for the first strain
      maxl = pars.maxi;
  elseif X == 2 %X = 2; %only plot for the second strain
      maxm = pars.maxj;
  elseif X == 3 %X = 3; %only plot for the third strain
      maxn = pars.maxk;
  end
  
  for i=1:length(nots)
          for l=1:maxl
              a = agestates;
              for m = maxm; %1:pars.maxj
                 for n = maxn; %1:pars.maxk
                    %i=times; l=immune status for strain1
                    if X == 1
                        titremat(i,l) = titremat(i,l) + y(i,pars.arrSlu(a,l,m,n));
                    %elseif X == 2
                    %    titremat(i,m) = titremat(i,m) + y(i,pars.arrSlu(a,l,m,n));
                    %elseif X == 3
                    %    titremat(i,n) = titremat(i,n) + y(i,pars.arrSlu(a,l,m,n));
                    end
                        %for X=1:pars.maxX
                    %    titremat(i,l) = titremat(i,l) + y(i,pars.arrIlu(X,a,l,m,n));
                    %end
                 end
              end
          end

  end
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

