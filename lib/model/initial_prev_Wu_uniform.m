function [ init_prev ] = initial_prev_Wu_uniform( pa, imm )
%INITIAL_PREV_XU_UNIFORM Summary of this function goes here
%   Detailed explanation goes here
        init_prev = []; 
        for a=1:4
        if a<4
            %seroconvert_obs = 0.033; %because MN1:40 = 3.3% http://cid.oxfordjournals.org/content/51/10/1184.full
            prev = [0.5 0.5 0.5 0.5];
            if isfield(pa,'inittitres') 
                seroconvert_obs = pa.inittitres;   %because HAI1:40 will be less than MN1:40
            else
                seroconvert_obs = 0.033;
            end
            naive = 1 - sum(prev)*seroconvert_obs;
            init_prev(a,:) = [naive prev*seroconvert_obs zeros(1,pa.maxi-length(prev)-1)];
            if pa.inittitres_flag == 3
                init_prev(a,:) = pa.init_prev;
                %init_prev(a,:) = [0.7 0.7^2 0.7^3 0.7^5 0.7^5 0.7^6 0.7^7 0.7^8 0.7^9 0.7^10];
                %init_prev(a,10) = 1-sum(init_prev(a,1:9));
            end
        end
        if a==4
            prev = [0.5 0.5 0.5 0.5];
            %seroconvert_obs = 0.033*2;
            seroconvert_obs = 0.02*2;
            if isfield(pa,'inittitres') 
                seroconvert_obs = pa.inittitres*2;   %because HAI1:40 will be less than MN1:40
                    if exist('imm')
                        seroconvert_obs = pa.inittitres*imm;
                    end
                    if pa.inittitres_flag == 2
                        seroconvert_obs = pa.inittitres;
                    end
            else
                seroconvert_obs = 0.033*2;
            end
            %seroconvert_obs
            naive = 1 - sum(prev)*seroconvert_obs;
            init_prev(a,:) = [naive prev*seroconvert_obs zeros(1,pa.maxi-length(prev)-1)];
            if pa.inittitres_flag == 3
                init_prev(a,:) = pa.init_prev;
            end
        end
        end
end

