function [ par Par_stat ] = initialModel( model, par, nsteps)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Par_stat = [];
PriorsPDF = {};
mode = num2str(model);
propdistr  = {''};  %proposal distribution
discrete = 0;       %index of dicrete variable

mint = 0;           %Starting titres = 0 -> undetectable 
maxt = 9;           %Maximum titres level
maxi = 10;          %Total titres levels for strainX
seed = 1;           %Starting infected seed
OutbreakNDA = 0;    %Disease starting day   
age_mix_flag = 1;
age_flag = 1;       %Age structure flag
frac_flag = 1;      %Age specific contact proportion flag
ContFrac1 = 4;  
immune_flag = 1;    %Age specific immune protection flag
%PUAb = 0.165;      %Assuming susceptibility of S0 drops 0.165 from S1 
PUAb = 0.3;
inittitres_flag = 1;%Age specific initial titres flag; add after 20150312; used in make_ics_naive.m
inittitres = 0.01;  %1 % seropositive
errp = 0;           %measurement error probability
    
%My ideal full model
%Doesn't work, immune_beta keep increasing
if strcmp(mode,'0.1')
    x0 =     [0.06,  6,     5,     4,     5,     2,      2,         2,        2     5];
    var =    [0.005  0.15   0.1    0.1    0.18   0.2     0.2        0.2       0.3   0.05];
    lb =     [0.01   1      1      1      1      0.5     0.5        0.5       0.5   1];
    ub =     [1      8      8      8      8      10      10         10        10    35];
    x_name = {'beta','AbB1','AbB2','AbB3','AbB4','immune_alpha1','immune_alpha2','immune_alpha3','immune_alpha4','immune_beta'};
    distr  = {'uniform','uniform','uniform','uniform','uniform','uniform','uniform','uniform','uniform','uniform'};
    optb =   ub;
    
    age_flag = 1;       %Age structure flag
    frac_flag = 1;      %Age specific contact proportion flag
    ContFrac1 = 4;  
    immune_flag = 1;    %Age specific immune protection flag
    PUAb = 0.165;       %Assuming susceptibility of S0 drops 0.165 from S1 
    inittitres_flag = 1;%Age specific initial titres flag; add after 20150312; used in make_ics_naive.m
end

%% The simplest model with no age mixing
%SIR
%-age strcuture

%% The main model in paper
%SIR
%+age strcuture
%+age specific boosting
%with PUAb 0.3
%with fixed protection curve shape

%predicting
if strcmp(mode,'1.3')
    setTitresModelT1();
    seed = 10;
end

if strcmp(mode,'1.4')
    setTitresModelT1_2();
    seed = 10;
end

if strcmp(mode,'1.22')
    setTitresModelT1_3();
    seed = 10;
end


%full titres model
if strcmp(mode,'1')
    setTitresModel();
    %setTitresModelFix();
    %setTitresModelImmune();
    %setTitresModelBoost();
    seed = 10;
end

%full titres model
if strcmp(mode,'3')
    %setTitresModel();
    %setTitresModelFix();
    %setTitresModelImmune();
    setTitresModelBoost();
    seed = 10;
end

%full titres model
if strcmp(mode,'4')
    %setTitresModel();
    %setTitresModelFix();
    setTitresModelImmune();
    %setTitresModelBoost();
    seed = 10;
end


if strcmp(mode,'1.11') %Change seed -> 1
    setTitresModel();
    seed = 1;
end

if strcmp(mode,'1.12') %Change seed -> 3
    setTitresModel();
    seed = 3;
end

if strcmp(mode,'1.13') %Change seed -> 30
    setTitresModel();;
    seed = 30;
end

if strcmp(mode,'1.14') %Change seed -> 100
    setTitresModel();
    seed = 100;
end

if strcmp(mode,'1.15') %Change seed -> 100
    setTitresModel();
    seed = 1000;
end

if strcmp(mode,'1.16') %Change seed -> 100
    setTitresModel();
    seed = 5;
end

if strcmp(mode,'1.17') %Change seed -> 100
    setTitresModel();
    seed = 20;
end

if strcmp(mode,'1.2') %Change seed -> 100
    setTitresModelF1();
    seed = 10;
end



if strcmp(mode,'1.5')
    setTitresModel();
    seed = 10;
    inittitres_flag = 2;
    inittitres = 0.02;
end

if strcmp(mode,'1.6')
    setTitresModel();
    seed = 10;
    inittitres = 0.08;
end


%% The alternative model with age mixing
%SIR
%+age strcuture
%with PUAb 0.3
%with fixed protection curve shape

if strcmp(mode,'2.2')
    setThresholdModelT1();
    seed = 10;
end

if strcmp(mode,'2')
    %setThresholdModel();
    setThresholdModelFix();
    seed = 10;
end

if strcmp(mode,'2.11')
    setThresholdModel();
    seed = 1;
end

if strcmp(mode,'2.12')
    setThresholdModel();
    seed = 3;
end

if strcmp(mode,'2.13')
    setThresholdModel();
    seed = 30;
end

if strcmp(mode,'2.14')
    setThresholdModel();
    seed = 100;
end

if strcmp(mode,'2.15')
    setThresholdModel();
    seed = 1000;
end

if strcmp(mode,'2.151')
    setThresholdModel();
    seed = 20000;
end

if strcmp(mode,'2.16')
    setThresholdModel();
    seed = 5;
end

if strcmp(mode,'2.17')
    setThresholdModel();
    seed = 20;
end

if strcmp(mode,'2.3')
    setThresholdModel();
    seed = 10;
    PUAb = 0.05;
end
if strcmp(mode,'2.311')
    setThresholdModel();
    seed = 1;
    PUAb = 0.3;
end

if strcmp(mode,'2.312')
    setThresholdModel();
    seed = 3;
    PUAb = 0.3;
end

if strcmp(mode,'2.313')
    setThresholdModel();
    seed = 30;
    PUAb = 0.3;
end

if strcmp(mode,'2.314')
    setThresholdModel();
    seed = 100;
    PUAb = 0.3;
end

if strcmp(mode,'2.4')
    setThresholdModelT1();
    seed = 10;
    PUAb = 0.3;
end

if strcmp(mode,'2.5')
    setThresholdModel();
    seed = 10;
    inittitres_flag = 2;
    inittitres = 0.02;
end


%% Update date parameters flags
%serological parameters
    par = setParameters(par,'mint',mint);          %Starting titres = 0 -> undetectable 
    par = setParameters(par,'maxt',maxt);          %Maximum titres level
    par = setParameters(par,'maxi',maxi);          %Total titres levels for strainX
    par = setParameters(par,'seed',seed);          %Starting infected seed
    par = setParameters(par,'OutbreakNDA',OutbreakNDA);

%age specefic 
    par = setParameters(par,'age_mix_flag',age_mix_flag);
    par = setParameters(par,'age_flag',age_flag);    %Age structure flag
    par = setParameters(par,'frac_flag',frac_flag);   %Age specific contact proportion flag
    par = setParameters(par,'ContFrac1',ContFrac1);  
    par = setParameters(par,'immune_flag',immune_flag); %Age specific immune protection flag
    par = setParameters(par,'PUAb',PUAb);   %Assuming susceptibility of S0 drops 0.165 from S1 
    par = setParameters(par,'inittitres_flag',inittitres_flag); %Age specific initial titres flag; add after 20150312; used in make_ics_naive.m
    par = setParameters(par,'inittitres',inittitres);
    
%measurement error
    par = setParameters(par,'errp',errp);
    par.matY = observe_matrix(par.maxi, par.errp);
    
%% Update model metadata
    Par_stat.model = model;
    Par_stat.ode.opt.A = diag(ones(1,length(x_name)));
    Par_stat.ode.opt.b = optb;
    Par_stat.ode.opt.x0 = x0;
    Par_stat.ode.opt.x_name = x_name;
    Par_stat.ode.opt.lb = lb;
    Par_stat.ode.opt.ub = ub;

    Par_stat.ode.mcmc.x0 =     x0;
    Par_stat.ode.mcmc.var =    var;
    Par_stat.ode.mcmc.lb =     lb;
    Par_stat.ode.mcmc.ub =     ub;
    Par_stat.ode.mcmc.x_name = x_name;
    Par_stat.ode.mcmc.distr  = distr;
    

    %Define Prior and Log Likelihood 
    for i=1:length(var)
      Prior(i).init = x0(i);
      Prior(i).name = char(x_name(i));
      Prior(i).distribution = char(distr(i));
      Prior(i).var = var(i);
      Prior(i).lb = lb(i);
      Prior(i).ub = ub(i);
      if(strcmp(Prior(i).distribution,'uniform'))
        Prior(i).pdf = @unifpdf;
        PriorPDF = @(x) Prior(i).pdf(x,Prior(i).lb,Prior(i).ub);
      end
      if(strcmp(Prior(i).distribution,'Normal'))
        Prior(i).norm_mu = mu(i);
        Prior(i).norm_sd = va(i);
        Prior(i).pdf = @truncate_normpdf;
        PriorPDF = @(x) Prior(i).pdf(x,Prior(i).norm_mu,Prior(i).norm_sd,Prior(i).lb,Prior(i).ub);
      end
      %PriorPDF = @(x) Prior(i).pdf(x,Prior(i).lb,Prior(i).ub);
      %%Call priorpdf in lib/model
      PriorsPDF(i) = {PriorPDF};
    end
    Par_stat.ode.mcmc.Prior = Prior;
    Par_stat.ode.mcmc.PriorsPDF = PriorsPDF;
    
    



function [] = setTitresModel()
        x0 =     [0.06,  6,     5,     4,     5,     2,      2,         2,        2       4   ];
        var =    [0.005  0.15   0.15   0.15   0.15   0.2     0.2        0.2       0.2     0.5 ];
        lb =     [0.01   1      1      1      1      0.5     0.5        0.5       0.5     1];
        ub =     [1      8      8      8      8      10      10         10        10      10];
        mu = x0;
        %va = var;
        va =     [0.005  0.15   0.15   0.15   0.15   0.2     0.2        0.2       0.2     0.5]; %used for normal prior
        x_name = {'beta','AbB1','AbB2','AbB3','AbB4','immune_alpha1','immune_alpha2','immune_alpha3','immune_alpha4','ContFrac1'};
        distr  = {'uniform','uniform','uniform','uniform','uniform','uniform','uniform','uniform','uniform','Normal'};
        optb =   ub;
        par = setParameters(par,'immune_beta',2.102); %still use 2.102 because the shape is reasonable
        mint = 0;                                   %starting titres = 0 -> undetectable 
        maxt = 9;                                   %maximum titres level
        maxi = 10;                                  %total titres levels for strainX
        seed = 10;                                   %Starting infected seed
        OutbreakNDA = 0;
        age_mix_flag = 1;
        age_flag = 1;
        frac_flag = 1;
        ContFrac1 = 4; %something is wrong when I use contfrac1 = 1 => no dynamics
        immune_flag = 1;
        PUAb = 0;
        inittitres_flag = 1;
        inittitres = 0.02;
        errp = 0.005;     
end

function [] = setTitresModelFix()
        x0 =     [0.06,  6,     5,     4,     5,     2,      2,         2,        2       ];
        var =    [0.005  0.15   0.15   0.15   0.15   0.2     0.2        0.2       0.2     ];
        lb =     [0.01   1      1      1      1      0.5     0.5        0.5       0.5     ];
        ub =     [1      8      8      8      8      10      10         10        10      ];
        mu = x0;
        %va = var;
        va =     [0.005  0.15   0.15   0.15   0.15   0.2     0.2        0.2       0.2     ]; %used for normal prior
        x_name = {'beta','AbB1','AbB2','AbB3','AbB4','immune_alpha1','immune_alpha2','immune_alpha3','immune_alpha4'};
        distr  = {'uniform','uniform','uniform','uniform','uniform','uniform','uniform','uniform','uniform'};
        optb =   ub;
        par = setParameters(par,'immune_beta',2.102); %still use 2.102 because the shape is reasonable
        mint = 0;                                   %starting titres = 0 -> undetectable 
        maxt = 9;                                   %maximum titres level
        maxi = 10;                                  %total titres levels for strainX
        seed = 10;                                   %Starting infected seed
        OutbreakNDA = 0;
        age_mix_flag = 1;
        age_flag = 1;
        frac_flag = 1;
        ContFrac1 = 1; %%%something is wrong when I use contfrac1 = 1 => no dynamics
        immune_flag = 1;
        PUAb = 0;
        inittitres_flag = 1;
        inittitres = 0.02;
        errp = 0.005;     
end



function [] = setTitresModelImmune()
    x0 =     [0.06,  6,     5,     4,     5,     2    ];
    var =    [0.005  0.15   0.15   0.15   0.15   0.2  ];
    lb =     [0.01   1      1      1      1      0.5  ];
    ub =     [1      8      8      8      8      10   ];
    
    mu = x0;
    va =     [0.005  0.15   0.15   0.15   0.15   0.5   0.5]; %used for normal prior  
    x_name = {'beta','AbB1','AbB2','AbB3','AbB4','immune_alpha'};
    distr  = {'uniform','uniform','uniform','uniform','uniform','uniform'};
    optb =   ub;
    par = setParameters(par,'immune_beta',2.102); %still use 2.102 because the shape is reasonable
    mint = 0;                                   %starting titres = 0 -> undetectable 
    maxt = 9;                                   %maximum titres level
    maxi = 10;                                  %total titres levels for strainX
    seed = 10;                                   %Starting infected seed
    OutbreakNDA = 0;
    age_mix_flag = 1;
    age_flag = 1;
    frac_flag = 1;
    ContFrac1 = 4; %something is wrong when I use contfrac1 = 1 => no dynamics
    immune_flag = 0;
    PUAb = 0;
    inittitres_flag = 1;
    inittitres = 0.02;
    errp = 0.005;    
end

function [] = setTitresModelBoost()
        x0 =     [0.06,  6,     2,      2,         2,        2       4   ];
        var =    [0.005  0.15   0.2     0.2        0.2       0.2     0.01 ];
        lb =     [0.01   1      0.5     0.5        0.5       0.5     1];
        ub =     [1      8      10      10         10        10      10];
        mu = x0;
        %va = var;
        va =     [0.005  0.15     0.2     0.2        0.2       0.2     0.5]; %used for normal prior
        x_name = {'beta','AbB','immune_alpha1','immune_alpha2','immune_alpha3','immune_alpha4','ContFrac1'};
        distr  = {'uniform','uniform','uniform','uniform','uniform','uniform','Normal'};
        optb =   ub;
        par = setParameters(par,'immune_beta',2.102); %still use 2.102 because the shape is reasonable
        mint = 0;                                   %starting titres = 0 -> undetectable 
        maxt = 9;                                   %maximum titres level
        maxi = 10;                                  %total titres levels for strainX
        seed = 10;                                   %Starting infected seed
        OutbreakNDA = 0;
        age_mix_flag = 1;
        age_flag = 0;
        frac_flag = 1;
        ContFrac1 = 4; %something is wrong when I use contfrac1 = 1 => no dynamics
        immune_flag = 1;
        PUAb = 0;
        inittitres_flag = 1;
        inittitres = 0.02;
        errp = 0.005;     
end


function [] = setTitresModelF1()
        x0 =     [0.06,  6,     5,     4,     5,     2,      2,         2,        2       1   ];
        var =    [0.005  0.15   0.15   0.15   0.15   0.2     0.2        0.2       0.2     0.5 ];
        lb =     [0.01   1      1      1      1      0.5     0.5        0.5       0.5     1];
        ub =     [1      8      8      8      8      10      10         10        10      10];
        mu = x0;
        %va = var;
        va =     [0.005  0.15   0.15   0.15   0.15   0.2     0.2        0.2       0.2     0.5]; %used for normal prior
        x_name = {'beta','AbB1','AbB2','AbB3','AbB4','immune_alpha1','immune_alpha2','immune_alpha3','immune_alpha4','ContFrac1'};
        distr  = {'uniform','uniform','uniform','uniform','uniform','uniform','uniform','uniform','uniform','Normal'};
        optb =   ub;
        par = setParameters(par,'immune_beta',2.102); %still use 2.102 because the shape is reasonable
        mint = 0;                                   %starting titres = 0 -> undetectable 
        maxt = 9;                                   %maximum titres level
        maxi = 10;                                  %total titres levels for strainX
        seed = 10;                                   %Starting infected seed
        OutbreakNDA = 0;
        age_mix_flag = 1;
        age_flag = 1;
        frac_flag = 1;
        ContFrac1 = 4; %something is wrong when I use contfrac1 = 1 => no dynamics
        immune_flag = 1;
        PUAb = 0;
        inittitres_flag = 1;
        inittitres = 0.02;
        errp = 0.005;     
end

function [] = setTitresModelT1()   %%%Maybe TP50 fixed to 1:40
        x0 =     [0.06    ];
        var =    [0.005   ];
        lb =     [0.01    ];
        ub =     [1       ];
        mu = x0;
        va =     [0.005 ]; %used for normal prior
        x_name = {'beta'};
        distr  = {'uniform'};
        optb =   ub;
        par = setParameters(par,'immune_beta',2.102); %still use 2.102 because the shape is reasonable
        %par = setParameters(par,'immune_alpha',3);    %still use 2.102 because the shape is reasonable
        mint = 0;                                     %starting titres = 0 -> undetectable 
        maxt = 9;                                     %maximum titres level
        maxi = 10;                                    %total titres levels for strainX
        seed = 10;                                    %Starting infected seed
        OutbreakNDA = 0;
        age_mix_flag = 1;
        age_flag = 1;
        frac_flag = 1;
        ContFrac1 = 4;   %something is wrong when I use contfrac1 = 1 => no dynamics
        immune_flag = 1; %different TP50 among each age 
        PUAb = 0;
        inittitres_flag = 1;
        inittitres = 0.02;
        errp = 0.005;     
        par = setParameters(par,'age_flag',1);
        par = setParameters(par,'immune_flag',1);
        par = setParameters(par,'AbB1',5.96);
        par = setParameters(par,'AbB2',4.97);
        par = setParameters(par,'AbB3',3.78);
        par = setParameters(par,'AbB4',4.79);
        par = setParameters(par,'immune_alpha1',2.15);
        par = setParameters(par,'immune_alpha2',3.40);
        par = setParameters(par,'immune_alpha3',2.80);
        par = setParameters(par,'immune_alpha4',5.76);
        ContFrac1 = 5.01;
end

function [] = setTitresModelT1_2()   %%%Maybe TP50 fixed to 1:40
        x0 =     [0.06    ];
        var =    [0.005   ];
        lb =     [0.01    ];
        ub =     [1       ];
        mu = x0;
        va =     [0.005 ]; %used for normal prior
        x_name = {'beta'};
        distr  = {'uniform'};
        optb =   ub;
        par = setParameters(par,'immune_beta',2.102); %still use 2.102 because the shape is reasonable
        %par = setParameters(par,'immune_alpha',3);    %still use 2.102 because the shape is reasonable
        mint = 0;                                     %starting titres = 0 -> undetectable 
        maxt = 9;                                     %maximum titres level
        maxi = 10;                                    %total titres levels for strainX
        seed = 10;                                    %Starting infected seed
        OutbreakNDA = 0;
        age_mix_flag = 1;
        age_flag = 1;
        frac_flag = 1;
        ContFrac1 = 4;   %something is wrong when I use contfrac1 = 1 => no dynamics
        immune_flag = 1; %different TP50 among each age 
        PUAb = 0;
        inittitres_flag = 1;
        inittitres = 0.02;
        errp = 0.005;
        par = setParameters(par,'age_flag',1);
        par = setParameters(par,'immune_flag',1);
        par = setParameters(par,'AbB1',6.0548);
        par = setParameters(par,'AbB2',4.8831);
        par = setParameters(par,'AbB3',4.0366);
        par = setParameters(par,'AbB4',3.7645);
        par = setParameters(par,'immune_alpha1',0.67587);
        par = setParameters(par,'immune_alpha2',3.0789);
        par = setParameters(par,'immune_alpha3',1.0268);
        par = setParameters(par,'immune_alpha4',4.5486);
        ContFrac1 = 6.5814;
end


function [] = setTitresModelT1_3()   %%%Maybe TP50 needs to be fixed
        x0 =     [0.06,  6,     5,     4,     5,     2,      2,         2,        2       4   ];
        var =    [0.005  0.15   0.15   0.15   0.15   0.2     0.2        0.2       0.2     0.5 ];
        lb =     [0.01   1      1      1      1      0.5     0.5        0.5       0.5     1];
        ub =     [1      8      8      8      8      10      10         10        10      10];
        mu = x0;
        va =     [0.005  0.15   0.15   0.15   0.15   0.2     0.2        0.2       0.2     0.5]; %used for normal prior
        x_name = {'beta','AbB1','AbB2','AbB3','AbB4','immune_alpha1','immune_alpha2','immune_alpha3','immune_alpha4','ContFrac1'};
        distr  = {'uniform','uniform','uniform','uniform','uniform','uniform','uniform','uniform','uniform','Normal'};
        optb =   ub;
        par = setParameters(par,'immune_beta',2.102); %still use 2.102 because the shape is reasonable
        mint = 0;                                   %starting titres = 0 -> undetectable 
        maxt = 9;                                   %maximum titres level
        maxi = 10;                                  %total titres levels for strainX
        seed = 10;                                   %Starting infected seed
        OutbreakNDA = 0;
        age_mix_flag = 1;
        age_flag = 1;
        frac_flag = 1;
        ContFrac1 = 4;   %something is wrong when I use contfrac1 = 1 => no dynamics
        immune_flag = 0; %same TP50 among each age 
        PUAb = 0;
        inittitres_flag = 1;
        inittitres = 0.02;
        errp = 0.005;     
end

        
function [] = setThresholdModel()
        x0 =     [0.05      4     ];
        var =    [0.001     0.1   ];
        lb =     [0.01      1     ];
        ub =     [1         10    ];
        %x0 =     [0.05  ];
        %var =    [0.01  ];
        %lb =     [0.01  ];
        %ub =     [1     ];
        mu = x0;
        va =     [0.001     0.5   ];;
        
        x_name = {'beta' 'ContFrac1'};
        distr  = {'uniform' 'Normal'};
        optb =   ub;
    
        par = setParameters(par,'immune_beta',50);  %step wise protection
        par = setParameters(par,'immune_alpha',0.5);%0/1 immunity
        par = setParameters(par,'AbB',5);
        
        mint = 0;                                   %starting titres = 0 -> undetectable 
        maxt = 1;                                   %maximum titres level
        maxi = 2;                                   %total titres levels for strainX
        seed = 10;                                   %Starting infected seed
        OutbreakNDA = 0;
        age_mix_flag = 1;                           %With Age mixing
        age_flag = 0;                               %No Age specific boosing
        frac_flag = 1;
        ContFrac1 = 4; 
        immune_flag = 0;
        PUAb = 0;                                 %For SIR model, set PUAb = 0
        inittitres_flag = 1;                        %For SIR model, set inittitres_flag = 0  
        inittitres = 0.02;
        errp = 0.005;
end

function [] = setThresholdModelFix()
        x0 =     [0.05      ];
        var =    [0.001     ];
        lb =     [0.01      ];
        ub =     [1         ];
        %x0 =     [0.05  ];
        %var =    [0.01  ];
        %lb =     [0.01  ];
        %ub =     [1     ];
        mu = x0;
        va =     [0.001     ];;
        
        x_name = {'beta' };
        distr  = {'uniform' };
        optb =   ub;
    
        par = setParameters(par,'immune_beta',50);  %step wise protection
        par = setParameters(par,'immune_alpha',0.5);%0/1 immunity
        par = setParameters(par,'AbB',5);
        
        mint = 0;                                   %starting titres = 0 -> undetectable 
        maxt = 1;                                   %maximum titres level
        maxi = 2;                                   %total titres levels for strainX
        seed = 10;                                   %Starting infected seed
        OutbreakNDA = 0;
        age_mix_flag = 1;                           %With Age mixing
        age_flag = 0;                               %No Age specific boosing
        frac_flag = 1;
        ContFrac1 = 4; 
        immune_flag = 0;
        PUAb = 0;                                 %For SIR model, set PUAb = 0
        inittitres_flag = 1;                        %For SIR model, set inittitres_flag = 0  
        inittitres = 0.02;
        errp = 0.005;
end


function [] = setThresholdModelT1()
        setThresholdModel();
end

function [] = setModel1()
        x0 =     [0.1  ];
        var =    [0.05 ];
        lb =     [0.01  ];
        ub =     [1     ];
        x_name = {'beta'};
        distr  = {'uniform'};
        optb =   ub;
    
        par = setParameters(par,'immune_beta',50);  %Step wise protection
        par = setParameters(par,'immune_alpha',0.5);%0/1 immunity
        par = setParameters(par,'AbB',5);
        
        mint = 0;                                   %Starting titres = 0 -> undetectable 
        maxt = 1;                                   %Maximum titres level
        maxi = 2;                                   %Total titres levels for strainX
        seed = 10;                                   %Starting infected seed
        OutbreakNDA = 0;
        age_mix_flag = 0;                           %No Age mixing
        age_flag = 0;                               %No Age specific boosing
        frac_flag = 0;                              %No Age specific virulence
        ContFrac1 = 1;                              %Children virulence=1
        immune_flag = 0;                            %No Age specific antibody protection
        PUAb = 0; % For SIR model, set PUAb = 0     %No protection for Ab undetected group
        inittitres_flag = 1;                        %Set initial titres for SIR model (set inittitres_flag=0)  
end

end


