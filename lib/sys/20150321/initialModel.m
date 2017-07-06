function [ par Par_stat ] = initialModel( model, par, nsteps)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Par_stat = [];
PriorsPDF = {};
mode = num2str(model);

%My ideal model
%Doesn't work, immune_beta keep increasing
if strcmp(mode,'1')
    x0 =     [0.06,  6,     5,     4,     5,     2,      2,         2,        2     5];
    var =    [0.005  0.15   0.1    0.1    0.18   0.2     0.2        0.2       0.3   0.05];
    lb =     [0.01   1      1      1      1      0.5     0.5        0.5       0.5   1];
    ub =     [1      8      8      8      8      10      10         10        10    35];
    x_name = {'beta','AbB1','AbB2','AbB3','AbB4','immune_alpha1','immune_alpha2','immune_alpha3','immune_alpha4','immune_beta'};
    distr  = {'uniform','uniform','uniform','uniform','uniform','uniform','uniform','uniform','uniform','uniform'};
    optb =   ub;
    par = setParameters(par,'age_flag',1);
    par = setParameters(par,'frac_flag',1);
    par = setParameters(par,'ContFrac1', 4); %contfrac1 = 4 still won't produce good immune_alpha
    par = setParameters(par,'immune_flag',1);
    par = setParameters(par,'PUAb', 0.165); % Assuming susceptibility of S0 drops 0.165 from S1 
    par = setParameters(par,'inittitres_flag', 1); %add after 20150312
end

%My 2nd ideal model
%with inital condition
%with PUAb 0
%with fixed protection curve shape
if strcmp(mode,'2')
    x0 =     [0.06,  6,     5,     4,     5,     2,      2,         2,        2    ];
    var =    [0.005  0.15   0.15   0.15   0.15   0.2     0.2        0.2       0.2  ];
    lb =     [0.01   1      1      1      1      0.5     0.5        0.5       0.5  ];
    ub =     [1      8      8      8      8      10      10         10        10   ];
    x_name = {'beta','AbB1','AbB2','AbB3','AbB4','immune_alpha1','immune_alpha2','immune_alpha3','immune_alpha4'};
    distr  = {'uniform','uniform','uniform','uniform','uniform','uniform','uniform','uniform','uniform'};
    optb =   ub;
    par = setParameters(par,'age_flag',1);
    par = setParameters(par,'frac_flag',1);
    par = setParameters(par,'ContFrac1', 4); %contfrac1 = 4 still won't produce good immune_alpha
    par = setParameters(par,'immune_flag',1);
    par = setParameters(par,'PUAb', 0);
    par = setParameters(par,'immune_beta',2.102);
    par = setParameters(par,'inittitres_flag', 1); %add after 20150312
end

if strcmp(mode,'3')
    x0 =     [0.06,  6,     5,     4,     5,     2,      2,         2,        2     5];
    var =    [0.005  0.15   0.15   0.15   0.15   0.2     0.2        0.2       0.2   0.05];
    lb =     [0.01   1      1      1      1      0.5     0.5        0.5       0.5   1];
    ub =     [1      8      8      8      8      10      10         10        10    35];
    x_name = {'beta','AbB1','AbB2','AbB3','AbB4','immune_alpha1','immune_alpha2','immune_alpha3','immune_alpha4','immune_beta'};
    distr  = {'uniform','uniform','uniform','uniform','uniform','uniform','uniform','uniform','uniform','uniform'};
    optb =   ub;
    par = setParameters(par,'age_flag',1);
    par = setParameters(par,'frac_flag',1);
    par = setParameters(par,'ContFrac1', 4); %contfrac1 = 4 still won't produce good immune_alpha
    par = setParameters(par,'immune_flag',1);
    par = setParameters(par,'PUAb', 0.3); % Assuming susceptibility of S0 drops 0.3 from 100% 
    par = setParameters(par,'inittitres_flag', 1); %add after 20150312
end

%% The main model in paper
if strcmp(mode,'4')
    x0 =     [0.06,  6,     5,     4,     5,     2,      2,         2,        2    ];
    var =    [0.005  0.15   0.15   0.15   0.15   0.2     0.2        0.2       0.2  ];
    lb =     [0.01   1      1      1      1      0.5     0.5        0.5       0.5  ];
    ub =     [1      8      8      8      8      10      10         10        10   ];
    x_name = {'beta','AbB1','AbB2','AbB3','AbB4','immune_alpha1','immune_alpha2','immune_alpha3','immune_alpha4'};
    distr  = {'uniform','uniform','uniform','uniform','uniform','uniform','uniform','uniform','uniform'};
    optb =   ub;
    par = setParameters(par,'age_flag',1);
    par = setParameters(par,'frac_flag',1);
    par = setParameters(par,'ContFrac1', 4); %contfrac1 = 4 still won't produce good immune_alpha
    par = setParameters(par,'immune_flag',1);
    par = setParameters(par,'PUAb', 0.3);
    par = setParameters(par,'immune_beta',2.102); %still use 2.102 because the shape is reasonable
    par = setParameters(par,'inittitres_flag', 1); %add after 20150312
end

%% The alternative model without age specific rates
if strcmp(mode,'5')
    x0 =     [0.06,  6,         2    ];
    var =    [0.005  0.15       0.2  ];
    lb =     [0.01   1          0.5  ];
    ub =     [1      8          10   ];
    x_name = {'beta','AbB','immune_alpha'};
    distr  = {'uniform','uniform','uniform'};
    optb =   ub;
    
    %% The following need to be changed 20150321
    par = setParameters(par,'age_flag',1);
    par = setParameters(par,'frac_flag',1);
    par = setParameters(par,'ContFrac1', 4); %contfrac1 = 4 still won't produce good immune_alpha
    par = setParameters(par,'immune_flag',1);
    par = setParameters(par,'PUAb', 0.3);
    par = setParameters(par,'immune_beta',2.102); %still use 2.102 because the shape is reasonable
    par = setParameters(par,'inittitres_flag', 1); %add after 20150312
end

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
    

    % longer test
    Par_stat.ode.mcmc.nSamples = 1000;
    Par_stat.ode.mcmc.burnIn = 100;
    if exist('nsteps')
       Par_stat.ode.mcmc.nSamples = nsteps;
       Par_stat.ode.mcmc.burnIn = 100;
    end

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
    
end

