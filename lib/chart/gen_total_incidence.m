function [ naive_inc imm_inc] = gen_total_incidence( y, times, pars)
%gen_total_infecteds Summary of this function goes here
%   Detailed explanation goes here
% y: model output
% times: sampling time
% pars: parameters
% a: age group

    dims_Slu = prod(size(pars.arrSlu)); %[a i j x]
    dims_Ilu = prod(size(pars.arrIlu)); %[X a i j k]
    %acc_incidence = sum(y(:,dims_Slu+1:dims_Slu+dims_Ilu),2); % total infected number, regardless which strains
    % add the number of infecteds by strain

    naive_inc = zeros(1,pars.maxa); %[1 x maxa]
    imm_inc = zeros(1,pars.maxa);   %[1 x maxa]
    
    %assign the number for accumulated infected
    yt = y(times,:);
    for X=1:pars.maxX
        for a=1:pars.maxa
            for i=1:pars.maxi
                for j=1:pars.maxj
                    for k=1:pars.maxk
                        if i>1  %titres > 1:5
                           imm_inc(1,a) = imm_inc(1,a) + yt(pars.arrCIlu(X,a,i,j,k));
                        end
                        if i==1 %titres = 1:5
                           naive_inc(1,a) = naive_inc(1,a) + yt(pars.arrCIlu(X,a,i,j,k));
                        end
                    end
                end
            end
        end
    end
    
    
end

