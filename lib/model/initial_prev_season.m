function [ init_prev ] = initial_prev_season( pa)
% Summary of this function goes here
%   Detailed explanation goes here
% Used for seasonal influenza serological data (based on season)
% The function is called in make_ics_naive() to generate serological titre 
        init_count = [];
        init_prev = []; 
        age = [0 20 40 65 110];
        for ai=1:(4)
            abl = pa.Antibody.K(pa.initK).Abl
            age_idx = find(pa.Antibody.K(pa.initK).age >= age(ai) & pa.Antibody.K(pa.initK).age < age(ai+1)); 
            Ablage = abl(age_idx);
            for ti = 1:pa.maxi
                init_count(ai,ti) = length(find(Ablage == ti-1));
            end
        end    
        %normalize
        for ai=1:4
          init_prev(ai,:) = init_count(ai,:)./sum(init_count(ai,:));
        end
end

