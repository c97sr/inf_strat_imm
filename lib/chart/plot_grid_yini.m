function [ serodist_mat ] = plot_grid_yini( y, pars, times, agegroup, yini )
%plot_grid Summary of this function goes here
% plot the antibody titres for susceptible individuals along the time
% Written by Sean Yuan (hyuan@imperial.ac.uk)

arrSlu = pars.arrSlu;
arrIlu = pars.arrIlu;

%Setup some standard auxilliary functions
tmin = min(times);
tmax = max(times);

%define age groups
arrIlu = pars.arrIlu; %I lookup
arrSlu = pars.arrSlu; %S lookup
serodist = zeros(length(times),1);
serodist_mat = []; % [age x times x titre]
for a = 1:pars.maxa
    %serodist = gen_strain_titres(y(times,:), 1:length(times), a, pars); % [times x titre]
    serodist = gen_strain_titres(y, times, a, pars); % [times x titre]
    serodist_ab = gen_strain_titres(yini, 1, a, pars);
    serodist_ab_t = repmat(serodist_ab,length(serodist),1);
    serodist_t = serodist - serodist_ab_t;
    serodist_mat(a,:,:) = reshape(serodist - serodist_ab_t, [1, size(serodist)]);  % [age x times x titre]
end

% Total
if ~exist('agegroup')
    agegroup = 1:pars.maxa;
end
serodist_tot3d = sum(serodist_mat(agegroup,:,:),1);  % [1 x times x titre]
serodist_tot = reshape(serodist_tot3d, size(serodist_tot3d,2), []); % [times x titre]
%normialze 
%nf = sum(serodist_tot,2);
%serodist_tot = serodist_tot./nf(1);

% example data
if pars.maxi>2
titres = length(serodist_tot(1,:)); % number of titre levels 
x1 = linspace(1,length(times),length(times)); % return the index of x on 2D grid
y1 = linspace(1,titres,titres); % return the index of y on 2D grid
%y = [0.1 0.5 3 3.2 3.5 6 7 8 9];
[X,Y] = meshgrid(x1,y1);
%heat = X.*Y; % titre values
serodist_tot(:,end+1) = serodist_tot(:,end); %for using surf to plot, the last column should duplicate otherwise
                                             %will miss the last column    
heat = serodist_tot(:,2:end)';
h = surf(X,Y,heat);
%h = mesh(X,Y,heat);
set(h, 'edgecolor','none')
view(0,90);
colormap(hot);
%colorbar;
set(gca,'xlim',[1 30.5*10]);
set(gca,'ylim',[1 titres]);
%set(gca,'YTick',[1.5:1:9.5]);
set(gca,'YTick',[1.5:1:pars.maxt+0.5]);
set(gca,'YTickLabel',['1';'2';'3';'4';'5';'6';'7';'8';'9']);
end

if pars.maxi == 2
titres = length(serodist_tot(1,:)); % number of titre levels 
x1 = linspace(1,length(times),length(times)); % return the index of x on 2D grid
y1 = linspace(1,pars.maxi+1,pars.maxi+1); % return the index of y on 2D grid
[X,Y] = meshgrid(x1,y1);
%heat = X.*Y; % titre values
serodist_tot(:,end+1) = serodist_tot(:,end); %for using surf to plot, the last column should duplicate otherwise
                                             %will last column won't be plotted   
heat = serodist_tot(:,1:end)';

% the plot
%figure;
h = surf(X,Y,heat);
set(h, 'edgecolor','none')
view(0,90);
colormap(hot);
%colorbar;
set(gca,'xlim',[1 30.5*10]); %xlim = 10 months
set(gca,'ylim',[1 titres+1]);
%set(gca,'YTick',[1.5:1:9.5]);
set(gca,'YTick',[1.5:1:pars.maxi+0.5]);
set(gca,'YTickLabel',['0';'1']);
end
    
    
  

%xlabel('Days');
ylabel('Antibody titres');
end

function [ titremat] = gen_strain_titres__( y, times, agestates, pars )
  nots = times;
  nostates = pars.maxi;
  titremat = zeros(length(nots), nostates);
  %X = 1; %only plot for the first strain
  for i=1:length(nots)
      for l=1:nostates
          a = agestates;
          for m = 1:pars.maxj
              for n = 1:pars.maxk
                titremat(i,l) = titremat(i,l) + y(i,pars.arrSlu(a,l,m,n));
                for X=1:pars.maxX
                    titremat(i,l) = titremat(i,l) + y(i,pars.arrIlu(X,a,l,m,n));
                end
              end
          end
      end
  end
end
