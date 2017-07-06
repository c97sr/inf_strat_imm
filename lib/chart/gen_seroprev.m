function [recovered_tot] = gen_seroprev( y, pars, times, strain, agegroup, display, marker )
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
    serodist= gen_strain_titres(y, times, a, pars); % [times x titre]
    serodist_mat(a,:,:) = reshape(serodist, [1, size(serodist)]);  % [age x times x titre]
end

% Total
if ~exist('agegroup')
    agegroup = 1:pars.maxa;
end
serodist_tot3d = sum(serodist_mat(agegroup,:,:),1);  % [1 x times x titre]
serodist_tot = reshape(serodist_tot3d, size(serodist_tot3d,2), []); % [times x titre]
if times == 1
    serodist_tot = sum(serodist_mat);
    susceptible_tot = sum(serodist_tot(1:abl-1)); % [times x 1], sum of individuals for titres <1:40; check abl=4
    recovered_tot = sum(serodist_tot(abl:pars.maxi)); % [times x 1], sum of individuals for titres >=1:40
else
    susceptible_tot = sum(serodist_tot(:,1:abl-1),2); % [times x 1], sum of individuals for titres <1:40; check abl=4
    recovered_tot = sum(serodist_tot(:,abl:pars.maxi),2); % [times x 1], sum of individuals for titres >=1:40
    recovered_tot = recovered_tot - recovered_tot(1); % extract the baseline
end
end