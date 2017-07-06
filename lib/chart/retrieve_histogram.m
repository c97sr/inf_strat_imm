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
if length(agegroup)>1

    for a = 1:length(agegroup)
        serodist_age(a,:,:) = gen_strain_titres(y, sampletime, agegroup(a), pars); % [times x titre]
    end
    serodist = squeeze(sum(serodist_age));
else
    serodist = gen_strain_titres(y, sampletime, agegroup, pars); % [times x titre]
end

if length(sampletime) > 1
    sumsero = repmat(sum(serodist,2), 1, length(serodist(1,:)));
    serodist_pro_tmp = (serodist./sumsero);
else
    serodist_pro_tmp = (serodist/sum(serodist))';
end
%figure;
Y = serodist_pro_tmp;
end


function [ infecteds] = gen_total_infecteds( y, times, pars)
    dims_Slu = prod(size(pars.arrSlu));
    dims_Ilu = prod(size(pars.arrIlu));
    infecteds = sum(y(:,dims_Slu+1:dims_Slu+dims_Ilu),2);
end

