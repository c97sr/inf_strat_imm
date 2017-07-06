function [ yini age_arr prevalence] = make_ics_fromtitres(  pars, arrSlu, arrIlu, arrCIlu, Abl )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 %proj = 'out/20140317/herdimmunity_hk/';
 %load([proj 'titres.mat']);
% Have to input antibody data with age structure
yini = zeros(1,pars.novars);
pop = [4.1 13.5 69.1 13.2];

age_arr = get_popweights_4_hk_2013(pars);
seed_arr = (age_arr.*pars.seed)./pars.N; 

%serology prevalance 
data_bar = accumarray(Abl,1);
prev = data_bar/sum(data_bar);
while pars.maxi > length(prev)
  prev(end+1) = 0;
end

% if titres larger than predefined maximum titre, move them into
% the maximum group.
while pars.maxi < length(prev)
            prev(pars.maxi) = prev(pars.maxi)+prev(end);
            prev(end) = [];
            %disp('Some antibody levels are larger than maximum antibody defined in the model.');
end


% Same sero prevalence distribution amoung each age groups
prev_byage = prev*age_arr; %[titres x age] 
prevalence.age_arr = age_arr;
prevalence.prev = prev;
prevalence.prev_byage = prev_byage;

%determine the sero prevalence
for a=1:pars.maxa
    for i=1:pars.maxi
        for j=1:pars.maxj
            for k=1:pars.maxk
                yini(arrSlu(a,i,j,k)) = (age_arr(a) - seed_arr(a))*initial_prev(pars,a,i,j,k); %not finished yet
            end
        end
    end
end

%add infectous seeds
for a=1:pars.maxa
    for X=1:pars.maxX % this will cause total number not 1
        %X = 1;
        yini(arrIlu(X,a,1,1,1)) = seed_arr(a)./3;
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


function popweights = get_popweights_4_hk_2013(pa)
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



function init_prev = initial_prev(pa,a,i,j,k)
        init_prev = [];
        %no age structure
        if(a>=1 && a<=pa.maxa)
            init_prev = prev(i);
        else
            error('Problem in initial_prev');
        end
end

end

