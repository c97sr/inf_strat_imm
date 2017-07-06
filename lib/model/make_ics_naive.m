function [ yini age_arr] = make_ics_naive( pars)
%make_ics Summary of this function goes here
% Setep the initial conditions with everyone susceptible
% Written by Sean Yuan (hyuan@imperial.ac.uk) 
% Complete naive individuals
% initial values are set to be 0 <<-??
yini = zeros(1,pars.novars);
%pop = [16.42 30.29 39.4 13.2]*pars.N;
pop = pars.demographic*pars.N;
age_arr =  get_popweights_4_hk(pars,pop); %age groups ratio
seed_arr = (age_arr.*pars.seed)./pars.N; 


arrSlu = pars.arrSlu;
arrIlu = pars.arrIlu;
arrCIlu = pars.arrCIlu;

% use pre season (initK) as the pre-existing antibody
initS = initial_prev_season(pars);

% initialize susceptible numbers
for a=1:pars.maxa
    for i=1:pars.maxi
        for j=1:pars.maxj
            for k=1:pars.maxk
                yini(arrSlu(a,i,j,k)) = (age_arr(a) - seed_arr(a))*initS(a,i); %only for 1st strain
            end
        end
    end
end




% initialize infected numbers
for a=1:pars.maxa
    for X=1:pars.maxX % this will cause total number not 1
        %X = 1;
        yini(arrIlu(X,a,1,1,1)) = seed_arr(a);
    end
end


% naive seroprevalence 
function inits = initial_naive(par)
        inits = zeros(par.maxa,par.maxi);
        for a1=1:par.maxa
          for i1=1:par.maxi
             if i1==1
               inits(a1,i1) = 1;
             end
          end
        end
end


end