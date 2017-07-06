function [ multi_p multi_p_byage] = get_titres_prob( y, pars, times, sampletime )
%plot_line Summary of this function goes here
% plot the total susceptible individuals along the time
% return overall multi_p for antibody titres without age structures
% return multi_p_byage for antibody titres with age structures
% Written by Sean Yuan (hyuan@imperial.ac.uk)
% 3 Jul 2014

arrSlu = pars.arrSlu;
arrIlu = pars.arrIlu;
multi_p = zeros(1,1);
multi_p_byage = zeros(pars.maxa,pars.maxi);

%define age groups
arrIlu = pars.arrIlu; %I lookup
arrSlu = pars.arrSlu; %S lookup
serodist = zeros(length(times),1);
serodist_mat = []; % [age x times x titre]
for a = 1:pars.maxa
    serodist = gen_strain_titres(y, times, a, pars); % [times x titre]
    serodist_mat(a,:,:) = reshape(serodist, [1 size(serodist)]);  % [age x times x titre]
    multi_p_byage(a,:) = serodist(sampletime,:)./sum(serodist(sampletime,:));
    multi_p_byage(a,find(multi_p_byage(a,:)<0)) = 0;
end

%no age groups
serodist_tot3d = sum(serodist_mat,1);  % [1 x times x titre]
serodist_tot = reshape(serodist_tot3d, [size(serodist_tot3d,2) size(serodist_tot3d,3)]); % [times x titre]
multi_p = serodist_tot(sampletime,:)./sum(serodist_tot(sampletime,:));
%disp done;
end


function [ infecteds] = gen_total_infecteds( y, times, pars)
    dims_Slu = prod(size(pars.arrSlu));
    dims_Ilu = prod(size(pars.arrIlu));
    infecteds = sum(y(:,dims_Slu+1:dims_Slu+dims_Ilu),2);
end

