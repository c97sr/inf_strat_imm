function [ init_prev ] = initial_prev_seasonbytime( pa, imm )
%INITIAL_PREV_XU_UNIFORM Summary of this function goes here
%   Detailed explanation goes here
% Used serological data based on number of days. Replaced by initial_prev_season.   
        init_count = [];
        init_prev = []; 
        age = [0 20 40 65 110];
        for ai=1:(4)
            idx = find((pa.Antibody.numdays-pa.OutbreakStartingDay) < 0 & (pa.Antibody.numdays-pa.OutbreakStartingDay) > -150 & pa.Antibody.age >= age(ai) & pa.Antibody.age < age(ai+1)); 
            Ablage = pa.Antibody.Abl(idx);
            for ti = 1:pa.maxi
                init_count(ai,ti) = length(find(Ablage == ti-1));
            end
        end    
        %normalize
        for ai=1:4
         init_prev(ai,:) = init_count(ai,:)./sum(init_count(ai,:));
        end
end

