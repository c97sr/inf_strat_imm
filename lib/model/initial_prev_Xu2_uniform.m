function [ init_prev ] = initial_prev_Xu2_uniform( pa )
%INITIAL_PREV_XU_UNIFORM Summary of this function goes here
%   Detailed explanation goes here
        init_prev = []; 
        for a=1:4
        if a<4
            seroconvert_obs = 0.033;
            naive = 1 - 2*seroconvert_obs;
            init_prev(a,:) = [naive [0.5 0.5 0.5 0.5]*seroconvert_obs zeros(1,pa.maxi-3)];
        end
        if a==4
            seroconvert_obs = 0.033*2;
            naive = 1 - 2*seroconvert_obs;
            init_prev(a,:) = [naive [0.5 0.5 0.5 0.5]*seroconvert_obs zeros(1,pa.maxi-3)];
        end
        end
end

