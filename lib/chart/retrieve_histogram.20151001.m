function [ Y ] = retrieve_histogram( y, pars, times, sampletime, agegroup )
% plot_infecteds_distribution_deter for deterministic siluation. Summary of this function goes here
% plot serological distribution
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
    serodist = gen_strain_titres(y, times+1, a, pars); % [times x titre]
    %serodist_mat(a,:,:) = reshape(serodist, [1, size(serodist)]);  % [age x times x titre]
    serodist_mat(a,:,:) = serodist;
end

% Total
if ~exist('agegroup')
    agegroup = 1:pars.maxa;
end
serodist_tot3d = sum(serodist_mat(agegroup,:,:),1);  % [1 x times x titre]
if ndims(serodist_tot3d)>2
    serodist_tot = reshape(serodist_tot3d, size(serodist_tot3d,2), []); % [times x titre]
else
    serodist_tot = serodist_tot3d;
end
serodist_pro_tmp = serodist_tot(sampletime,:)./sum(serodist_tot(sampletime,:));
%serodist_pro_tmp = serodist_tot(sampletime,:) %without normalizing
serodist_pro = [];
for i=1:length(serodist_pro_tmp)
    
    if i > 1
        serodist_pro(i) = serodist_pro(i-1)+serodist_pro_tmp(i);
    else
        serodist_pro(i) = serodist_pro_tmp(i);
    end
end

%figure;
Y = serodist_pro_tmp;
end


function [ infecteds] = gen_total_infecteds( y, times, pars)
    dims_Slu = prod(size(pars.arrSlu));
    dims_Ilu = prod(size(pars.arrIlu));
    infecteds = sum(y(:,dims_Slu+1:dims_Slu+dims_Ilu),2);
end

