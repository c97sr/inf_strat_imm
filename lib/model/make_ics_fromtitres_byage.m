function [ yini age_arr prevalence] = make_ics_fromtitres_byage(  pars, arrSlu, arrIlu, arrCIlu, Abl, age )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% same as make_ics_fromtitres but with one extra age input variable 

 yini = zeros(1,pars.novars);
 pop = pars.demographic*pars.N;

age_arr = get_popweights_4_hk_2013(pars,pop);
seed_arr = (age_arr.*pars.seed)./pars.N; 

%serology prevalance 
prev = [];
agepop = [];
for a=1:pars.maxa
    if pars.mint == 0
        if length(find(age>=pars.ages(a,1) & age<pars.ages(a,2))) == 0
            data_bar = 0;
        else
            data_bar = accumarray(Abl(find(age>=pars.ages(a,1) & age<pars.ages(a,2)))+1,1);
            
            ab_init = Abl(find(age>=pars.ages(a,1) & age<pars.ages(a,2)))+1;
            for i = 1:pars.maxi
                data_bar(i) = sum(ab_init==i);
            end
            data_bar_k1 = data_bar./sum(data_bar);
        end
    elseif pars.mint == 1
        data_bar = accumarray(Abl(find(age>=pars.ages(a,1) & age<pars.ages(a,2))),1);
    end
    
    prev = data_bar/sum(data_bar);
    
    while pars.maxi > length(prev)
        prev(end+1,1) = 0;
    end

% if titres larger than predefined maximum titre, move them into
% the maximum group.
    while pars.maxi < length(prev)
            prev(pars.maxi) = prev(pars.maxi)+prev(end);
            prev(end) = [];
            %disp('Some antibody levels are larger than maximum antibody defined in the model.');
    end
    agepop(a).prev = prev;
    prev_byage(:,a) = prev*age_arr(a);
end


% Same sero prevalence distribution amoung each age groups
%prev_byage = [agepop.prev].*repmat(age_arr,pars.maxi,1);  %[titres x age] 
prevalence.age_arr = age_arr;
prevalence.prev = prev;
prevalence.prev_byage = prev_byage;

%determine the sero prevalence
for a=1:pars.maxa
    for i=1:pars.maxi
        for j=1:pars.maxj
            for k=1:pars.maxk %20140714
                %This is wrong, cause the total incidences less than the
                %correct number
                %yini(arrSlu(a,i,j,k)) = (age_arr(a) - seed_arr(a))*initial_prev(pars,a,i,j,k,prev_byage(:,a)); %not finished yet
                yini(arrSlu(a,i,j,k)) = initial_prev(pars,a,i,j,k,prev_byage(:,a)); %not finished yet
                if (i==1 & j==1 & k==1)
                    yini(arrSlu(a,i,j,k)) =  yini(arrSlu(a,i,j,k)) - seed_arr(a);
                end
            end
        end
    end
end

%add infectous seeds
for a=1:pars.maxa
    for X=1:pars.maxX % this will cause total number not 1
        %X = 1;
        %yini(arrIlu(X,a,1,1,1)) = seed_arr(a)./3;
        if X == 1
            yini(arrIlu(X,a,1,1,1)) = seed_arr(a);
        end
        %yini(arrIlu(X,a,1,1,1)) = 0;
    end
end

%% Below are subfunctions
function popweights = get_popweights_5_uk_2003(pa)
    popweights = zeros(1,pa.maxa);
    %pop = [3009303, 8731641, 19889635, 12707659, 8453916]; 
    if pa.maxa~=5
      %disp(['max number of age classes is ' num2str(pa.maxa)]);
      pop = pop(1:pa.maxa);
      popweights = pop./sum(pop);
    else
      %pop = [3009303, 8731641, 19889635, 12707659, 8453916]; 
      popweights = pop./sum(pop);
    end
end




function init_prev = initial_prev(pa,a,i,j,k,prv)
        init_prev = [];
        %no age structure
        if(a>=1 && a<=pa.maxa)
            init_prev = prv(i);
        else
            error('Problem in initial_prev');
        end
end

end

