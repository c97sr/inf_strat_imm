function [ yini age_arr S0_imm] = make_ics_naive2titres( pars, arrSlu, arrIlu, arrCIlu, age)
%make_ics Summary of this function goes here
% Setep the initial conditions with everyone susceptible
% Written by Sean Yuan (hyuan@imperial.ac.uk) 
% Complete naive individuals
% initial values are set to be 0 <<-??
yini = zeros(1,pars.novars);
pop = pars.demographic*pars.N;
age_arr =  get_popweights_4_hk(pars,pop); %age groups ratio
seed_arr = (age_arr.*pars.seed)./pars.N; 

% Setup s0_imm for PUAb > 0
initS = initial_prev_Wu_uniform(pars);
initS0_imm = zeros(1,pars.maxa);
S0_imm = zeros(1,pars.maxa);
for a = 1:pars.maxa
  initS0_imm(a) = initS(a,1).*pars.PUAb;
  S0_imm(a) = (age_arr(a) - seed_arr(a))*initS0_imm(a);  %S(age)*Ti(i)*PUAB
end

if pars.inittitres_flag == 1
    initS = zeros(pars.maxa, pars.maxi);
    %S0_imm = zeros(1,pars.maxa);
    sero_pos = pars.inittitres;
    initS = [1-sero_pos sero_pos; 1-sero_pos sero_pos; 1-sero_pos sero_pos; 1-sero_pos*2 sero_pos*2];
end

%if pars.inittitres_flag == 1
%    if pars.PUAb == 0.3
%           S0_imm = [0.0493    0.0832    0.1104    0.0346]; %same from make_ics_naive when PUAb=0.3 (copy the estimated value)
%    elseif pars.PUAb == 0
%           S0_imm = [0 0 0 0];
%    end
%end
if pars.inittitres_flag == 3
    initS = zeros(pars.maxa, pars.maxi);
    %S0_imm = zeros(1,pars.maxa);
    sero_pos = pars.inittitres;
    initS = [1-sero_pos sero_pos; 1-sero_pos sero_pos; 1-sero_pos sero_pos; 1-sero_pos sero_pos];
end

if pars.inittitres_flag == 0
    sero_pos = 0;
    initS = [1-sero_pos sero_pos; 1-sero_pos sero_pos; 1-sero_pos sero_pos; 1-sero_pos*2 sero_pos*2];
%    if pars.PUAb == 0.3
%           S0_imm = [0.0528    0.0891    0.1182    0.0399]; %same from make_ics_naive when PUAb=0.3 (copy the estimated value)
%    elseif pars.PUAb == 0
%           S0_imm = [0 0 0 0];
%    end
end

if pars.inittitres_flag == 2
    initS = zeros(pars.maxa, pars.maxi);
    %S0_imm = zeros(1,pars.maxa);
    sero_pos = pars.inittitres;
    initS = [1-sero_pos sero_pos; 1-sero_pos sero_pos; 1-sero_pos sero_pos; 1-sero_pos sero_pos];
end

% initialize susceptible numbers
for a=1:pars.maxa
    for i=1:pars.maxi
        for j=1:pars.maxj
            for k=1:pars.maxk
                %yini(arrSlu(a,i,j,k)) = (age_arr(a) - seed_arr(a))*initial_prev_naive(pars,a,i,j,k); %not finished yet
                yini(arrSlu(a,i,j,k)) = (age_arr(a) - seed_arr(a))*initS(a,i); %only for 1st strain
            end
        end
    end
end


% initialize infected numbers
for a=1:pars.maxa
    for X=1:pars.maxX % this will cause total number not 1
        %X = 1;
        %yini(arrIlu(X,a,1,1,1)) = seed_arr(a)./pars.maxa; %%%should I divided by pars.maxa???
        yini(arrIlu(X,a,1,1,1)) = seed_arr(a);
    end
end

end