function [Table2Left Table2Right p Inc] = table2( PosteriorSamples, pars, burnIn, samplesize) 
% Summary of the function table2
% Display the seroprevalence at T1 and T2 for observed and model output
% Written by Sean Yuan (hyuan@imperial.ac.uk) 
% Plot all age group information.


global proj Antibody;
init_collect = 1;
second_collect = 2;
third_collect = 3;
k = 2;
Ab = Antibody;
posterior = table2array(PosteriorSamples);

if exist('samplesize') == 0
    samplesize = 3;
end
if exist('burnIn') == 0
    burnIn = 1000;
end
post = mean(posterior(burnIn:end,:));
total = length(posterior(:,1))-burnIn;
idx = burnIn + round(rand(1, samplesize) * total);

for i=1:samplesize
%model:3
    vars = PosteriorSamples.Properties.VariableNames;
    for p=1:length(vars)
        if strcmpi('LLH',vars(p))
        else
           [pars] = setParameters(pars,char(vars(p)),posterior(idx(i),p));
        end
    end

%set parameters
beta = pars.beta;
AbB = [pars.AbB1 pars.AbB2 pars.AbB3 pars.AbB4];
immune_alpha = [pars.immune_alpha1 pars.immune_alpha2 pars.immune_alpha3 pars.immune_alpha4];
lastsamplingday = pars.SamplingLastDay;



%setup initial condition
%setup initial condition
if pars.maxi == 2 % only 2 titres
       [yini age_arr s0_imm] = make_ics_naive2titres( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu);
else
       [yini age_arr] = make_ics_naive( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu);
end
    
[yini_k1 age_arr_k1] = make_ics_fromtitres_byage( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, Ab.K(init_collect).Abl, Ab.K(init_collect).age);
[yini_k2 age_arr_k2] = make_ics_fromtitres_byage( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, Ab.K(k).Abl, Ab.K(k).age);


    %setep simulation time
    T0 = pars.OutbreakStartingDay;
    meanKdays(1) = mean(pars.Antibody.K(1).numdays - T0);
    meanKdays(2) = mean(pars.Antibody.K(2).numdays - T0);

    sample_time_K1 = round(meanKdays(1));
    sample_time_K2 = round(meanKdays(2));
    times = 0:1:lastsamplingday;
    sample_size_K1 = Ab.K(1).samplesize;
    for i1=1:pars.maxa
         lage = pars.ages(i1,1);
         uage = pars.ages(i1,2);
         sample_size_K1(end+1) = length(find(lage<=Ab.K(1).age & Ab.K(1).age<uage));
    end
    sample_size_K2 = Ab.K(2).samplesize;
    for i2=1:pars.maxa
         lage = pars.ages(i2,1);
         uage = pars.ages(i2,2);
         sample_size_K2(end+1) = length(find(lage<=Ab.K(2).age & Ab.K(2).age<uage));
    end

    %run simulation
    %initialize objects
    daysshift = 0;
    daysshift = pars.OutbreakNDA;
    
    %javaaddpath e:\workspace\MyJavaProject\bin\matlabjava.jar;
    javaaddpath(pars.javapath); %set ./java as default dir
    mepar_t2 = matlabjava.Parameters;
    meser_t2 = matlabjava.Serology;
    meser_t2.setParameters(mepar_t2);
    % set parameters
    meser_t2.updateParametersG(pars.arrg);
    meser_t2.updateParametersH(pars.arrh);
    meser_t2.updateParametersM(pars.matM);
    meser_t2.updateParametersBeta(pars.beta);  
    meser_t2.updateParameters('wan',pars.wan);
    meser_t2.updateParameters('s0_imm', pars.s0_imm);
    x0 = yini;  
    [t y] = ode23(@(t,x)odef_islmodjava(t,x, meser_t2), times, x0);  
    %[t y] = ode23(@(t,x)odef_islmod(t,x,pars), times, x0);
    clear 'mepar_t2';
    clear 'meser_t2';
    T = t;
    
    %mean titre changes
    for a=1:pars.maxa
    Xout_k1(a,:) = retrieve_histogram(y, pars, times, sample_time_K1+daysshift, a); % model output
    %Xoutput_k1_list(i,:) = Xout_k1;
    Xout_k2(a,:) = retrieve_histogram(y, pars, times, sample_time_K2+daysshift, a); % model output
    %Xoutput_k2_list(i,:) = Xout_k2;
    Yout_k0(a,:) = retrieve_histogram(yini, pars, times(1), 1, a);
    Yout_k1(a,:) = retrieve_histogram(yini_k1, pars, times(1), 1, a); % observed data
    Yout_k2(a,:) = retrieve_histogram(yini_k2, pars, times(1), 1, a); % observed data
    
    agetitres(a).Yout_k1(a,:) = retrieve_histogram(yini_k1, pars, times(1), 1, a); % observed data
    agetitres(a).Yout_k2(a,:) = retrieve_histogram(yini_k2, pars, times(1), 1, a); % observed data
    
    Xout_k1; % observed data
    Xout_k2; % observed data 
    seroprev(a).age(i,1) = sum(Yout_k0(a,4:end)); %T0 
    seroprev(a).age(i,2) = sum(Yout_k1(a,4:end)); %T1 Observed
    seroprev(a).age(i,3) = sum(Yout_k2(a,4:end)); %T2 Observed
    seroprev(a).age(i,4) = sum(Xout_k1(a,4:end)); %T1 Model output
    seroprev(a).age(i,5) = sum(Xout_k2(a,4:end)); %T2 Model output
    
    meantitre(a).age(i,1) = sum((Yout_k1(a,2:end).*[1:9]))./sum(Yout_k1(a,2:end)); %T1 Observed Mean Titres
    meantitre(a).age(i,2) = sum((Yout_k2(a,2:end).*[1:9]))./sum(Yout_k2(a,2:end)); %T2 Observed Mean Titres
    meantitre(a).age(i,3) = sum((Xout_k1(a,2:end).*[1:9]))./sum(Xout_k1(a,2:end)); %T1 Model output Mean Titres
    meantitre(a).age(i,4) = sum((Xout_k2(a,2:end).*[1:9]))./sum(Xout_k2(a,2:end)); %T2 Model output Mean Titres
    end
    
    
    %%%% why use atotal??? 20150316
    %seroconversion
    atotal=1:pars.maxa; %%%
    %total age group
    Xout_k1 = retrieve_histogram(y, pars, times, sample_time_K1+daysshift, atotal); % model output
    Xout_k2 = retrieve_histogram(y, pars, times, sample_time_K2+daysshift, atotal); % model output
    Yout_k0 = retrieve_histogram(yini, pars, times(1), 1, atotal);
    Yout_k1 = retrieve_histogram(yini_k1, pars, times(1), 1, atotal);     % observed data
    Yout_k2 = retrieve_histogram(yini_k2, pars, times(1), 1, atotal);     % observed data
    
    seroprev_total(i,1) = sum(Yout_k0(4:end));
    seroprev_total(i,2) = sum(Yout_k1(4:end)); %T1 Observed
    seroprev_total(i,3) = sum(Yout_k2(4:end)); %T2 Observed
    seroprev_total(i,4) = sum(Xout_k1(4:end)); %T1 Model output
    seroprev_total(i,5) = sum(Xout_k2(4:end)); %T2 Model output
    
    meantitre_total(i,1) = sum((Yout_k1(2:end).*[1:9]))./sum(Yout_k1(2:end)); %T1 Observed Mean Titres
    meantitre_total(i,2) = sum((Yout_k2(2:end).*[1:9]))./sum(Yout_k2(2:end)); %T2 Observed Mean Titres
    meantitre_total(i,3) = sum((Xout_k1(2:end).*[1:9]))./sum(Xout_k1(2:end)); %T1 Model output Mean Titres
    meantitre_total(i,4) = sum((Xout_k2(2:end).*[1:9]))./sum(Xout_k2(2:end)); %T2 Model output Mean Titres
    
end

    TotInc = seroprev_total(:,5) - seroprev_total(:,4);
    IncMean = mean(TotInc);
    IncUb = quantile(TotInc,0.975);
    IncLb = quantile(TotInc,0.025);
    Inc.Tot = TotInc;
    Inc.Mean = IncMean;
    Inc.UncUb = IncUb;
    Inc.UncLb = IncLb;
    
    %Seroprevalence
    %baseline = seroprev_total(1,1);
    precision = 1;
    Table(1,1:5) = mean(seroprev_total(:,:));% T0, initial, followup, modelT1, modelT2
    Table(1,6:10) = quantile(seroprev_total,0.025);% lower bound;
    Table(1,11:15) = quantile(seroprev_total,0.975);% upper bound;
    %Mean titres
    Table(1,16:19) = mean(meantitre_total(:,:));
    Table(1,20:23) = quantile(meantitre_total(:,:),0.025);
    Table(1,24:27) = quantile(meantitre_total(:,:),0.975);
    
    TableCI = []; % store 95 CI for sera data
    pbin = Table(1,2); %round1 prevalence
    %find lower bound
    y0 = pbin; %inital value
    sizeage = sample_size_K1(1); % total:523
    options = optimoptions('fsolve','Display','off'); % Option to display output
    lb = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.975','x','pbin','sizeage'),y0,[options],pbin,sizeage);
    %find upper bound
    ub = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.025','x','pbin','sizeage'),y0,[options],pbin,sizeage);
    TableCI(1,1:2) = [lb ub];
    
    pbin = Table(1,3); %round2 prevalence
    %find lower bound
    y0 = pbin; %inital value
    sizeage = sample_size_K2(1);
    lb = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.975','x','pbin','sizeage'),y0,[options],pbin,sizeage);
    %find upper bound
    ub = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.025','x','pbin','sizeage'),y0,[options],pbin,sizeage);
    TableCI(1,3:4) = [lb ub];
    
    %Geometric mean titres at round1
    meanT = meantitre_total(1,1);
    Fx = Yout_k1(2:end)/sum(Yout_k1(2:end));
    var = sum((((1:9)-meanT).^2).*Fx);
    SD = var.^0.5;
    TableCI(1,5) = SD;
    
    meanT = meantitre_total(1,2);
    Fx = Yout_k1(2:end)/sum(Yout_k2(2:end));
    var = sum((((1:9)-meanT).^2).*Fx);
    SD = var.^0.5;
    TableCI(1,6) = SD;
    
    
    
    
    for a=1:pars.maxa
        %Seroprevalence
        baseline_age = seroprev(a).age(1,1);
        Table(a+1,1:5) = mean(seroprev(a).age);             % - baseline_age;
        Table(a+1,6:10) = quantile(seroprev(a).age,0.025);  % - lower bound
        Table(a+1,11:15) = quantile(seroprev(a).age,0.975); % - upper bound;
        %Mean titres
        Table(a+1,16:19) = mean(meantitre(a).age(:,:));
        Table(a+1,20:23) = quantile(meantitre(a).age(:,:),0.025);
        Table(a+1,24:27) = quantile(meantitre(a).age(:,:),0.975);
        
        pbin = Table(a+1,2); %round1 prevalence
        %find lower bound
        y0 = pbin; %inital value
        sizeage = sample_size_K1(1+a); % sample size for each age group 
        lb = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.975','x','pbin','sizeage'),y0,[options],pbin,sizeage);
        %find upper bound
        ub = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.025','x','pbin','sizeage'),y0,[options],pbin,sizeage);
        TableCI(a+1,1:2) = [lb ub];
        
        pbin = Table(a+1,3); %round2 prevalence
        %find lower bound
        y0 = pbin; %inital value
        sizeage = sample_size_K2(1+a);  % sample size for each age group 
        lb = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.975','x','pbin','sizeage'),y0,[options],pbin,sizeage);
        %find upper bound
        ub = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.025','x','pbin','sizeage'),y0,[options],pbin,sizeage);
        TableCI(a+1,3:4) = [lb ub];
        
        %calculate SD for GMT at round 1
        meanT = meantitre(a).age(1,1);
        Fx = agetitres(a).Yout_k1(a,2:end)/sum(agetitres(a).Yout_k1(a,2:end));
        var = sum((((1:9)-meanT).^2).*Fx);
        SD = var.^0.5;
        TableCI(a+1,5) = SD;
        
        %calculate SD for GMT at round 2
        meanT = meantitre(a).age(1,2);
        Fx = agetitres(a).Yout_k2(a,2:end)/sum(agetitres(a).Yout_k2(a,2:end));
        var = sum((((1:9)-meanT).^2).*Fx);
        SD = var.^0.5;
        TableCI(a+1,6) = SD;
        
    end
    
 %Produce Table2 Left Seroprevalence
 Table2Left = array2table(zeros(5,9));
 x_name =  {'Initial','YR1','YR1CI','YR2','YR2CI','XT1','XT1CI','XT2','XT2CI'}; 
 pro = 100;
 precision = 1000;

 Table_sero = pro*(round(Table(:,1:15)*precision)/precision);
 TableCI_sero = pro*(round(TableCI(:,1:4)*precision)/precision);
 
 Table2Left.Properties.VariableNames = x_name;
 Table2Left.Properties.RowNames = {'Total','Age1','Age2','Age3','Age4'};
 Table2Left.Initial = Table_sero(:,1);
 Table2Left.YR1 = [num2str(Table_sero(:,2)) ];
 Table2Left.YR1CI = [repmat(' [',5,1) num2str(TableCI_sero(:,1)) repmat(' ',5,1) num2str(TableCI_sero(:,2)) repmat(']',5,1)];
 Table2Left.YR2 = [num2str(Table_sero(:,3)) ];
 Table2Left.YR2CI = [repmat(' [',5,1) num2str(TableCI_sero(:,3)) repmat(' ',5,1) num2str(TableCI_sero(:,4)) repmat(']',5,1)];
 Table2Left.YR2 = [num2str(Table_sero(:,3)) ]; 
 Table2Left.XT1 = [num2str(Table_sero(:,4)) ];
 Table2Left.XT1CI = [repmat(' [',5,1) num2str(Table_sero(:,9)) repmat(' ',5,1) num2str(Table_sero(:,14)) repmat(']',5,1)];
 Table2Left.XT2 = [num2str(Table_sero(:,5)) ]; 
 Table2Left.XT2CI = [repmat(' [',5,1) num2str(Table_sero(:,10)) repmat(' ',5,1) num2str(Table_sero(:,15)) repmat(']',5,1)]; 

 
 %Produce Table2 Right Meantitres
 pro = 1;
 precision = 10;
 Table_gmt = Table;
 TableCI_gmt = TableCI;
 Table_gmt(:,16:27) = pro*(round(Table(:,16:27)*precision)/precision);
 TableCI_gmt(:,5:6) = pro*(round(TableCI(:,5:6)*precision)/precision);
 Table2Right = array2table(zeros(5,8));
 x_name = {'YR1','YR1SD','YR2','YR2SD','XT1','XT1CI','XT2','XT2CI'};
 Table2Right.Properties.VariableNames = x_name;
 Table2Right.Properties.RowNames = {'Total','Age1','Age2','Age3','Age4'};
 Table2Right.YR1 = [num2str(Table_gmt(:,16)) ];
 Table2Right.YR1SD = [num2str(TableCI_gmt(:,5)) ];
 Table2Right.YR2 = [num2str(Table_gmt(:,17)) ];
 Table2Right.YR2SD = [num2str(TableCI_gmt(:,6)) ];
 Table2Right.XT1 = [num2str(Table_gmt(:,18)) ];
 Table2Right.XT1CI = [repmat(' [',5,1) num2str(Table_gmt(:,22)) repmat(' ',5,1) num2str(Table_gmt(:,26)) repmat(']',5,1)];
 Table2Right.XT2 = [num2str(Table_gmt(:,19)) ];
 Table2Right.XT2CI = [repmat(' [',5,1) num2str(Table_gmt(:,23)) repmat(' ',5,1) num2str(Table_gmt(:,27)) repmat(']',5,1)];


 Table2Left
 Table2Right
 writetable(Table2Left,'Table2a.csv');
 writetable(Table2Right,'Table2b.csv');

 disp 'keep walking';

end