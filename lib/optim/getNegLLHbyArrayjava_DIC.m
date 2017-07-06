function [nLLH LLH_Map_age]= getNegLLHbyArrayjava_DIC(theta, theta_name, par, meser, x0, yt_array)
%x: input parameter values
%x: input parameter names
%par: parameter object
%meser: sera model in java
%yini: initial data
%yt: array of yt (will improve efficiency)
    

    %Setup parameters
    OutbreakStartingDay = par.OutbreakStartingDay; %default 120 days 
    SamplingLastDay = par.SamplingLastDay;
    
    % Update parameters
    for i=1:length(theta)
        par = setParameters(par,theta_name{i},theta(i));
    end
    
    %Read Obserbed Yt
    %sampletime_k = par.currK;
    Antibody = par.Antibody;
    SampleSize = Antibody.samplesize;
    Abl = Antibody.Abl;
    age = Antibody.age;
    numdays = Antibody.numdays;
    Abl(find(Abl>par.maxt)) = par.maxt; %substitute Ab level >maxt
    %transform observed titres into 2D-array
    
    %SampleSize = length(Abl);
    %SampleTime = unique(numdays) - OutbreakStartingDay;                     

    
    
   %%
    % run ODE
    daysshift = 0;
    daysshift = par.OutbreakNDA;
    %if daysshift < 0 
    %  daysshift = 0;
    %end
    if daysshift< 0
        times = -daysshift + (0:1:SamplingLastDay); %Simulation time period; Set 120d = 0; Total 365d
    else
        %if days shift earlier, simulate longer time
        times = -daysshift + (0:1:SamplingLastDay+daysshift); %Simulation time period; Set 120d = 0; Total 365d
    end
    % set parameters
    meser.updateParametersG(par.arrg);
    meser.updateParametersH(par.arrh);
    meser.updateParametersM(par.matM);
    meser.updateParametersBeta(par.beta);
    meser.updateParameters('wan',par.wan);
    meser.updateParameters('s0_imm', par.s0_imm);
    %meser.updateParameters('seed', par.seed); 
    
    [t Xt] = ode45(@(t,x)odef_islmodjava(t,x, meser), times, x0);
    %adjusted by datsshift
    %Xt = Xt(times+1,:);   
    
    % define serology distribution matrix
    %serodist_mat = zeros(par.maxa,t(end),par.maxi); % [age x times x titre]
    %%T0:T0-60+515,:), pars, T_rel
    if daysshift< 0
      sh = abs(daysshift);
      Xt0 = [repmat([Xt(1,:)],sh,1); Xt(1:(length(Xt)-sh),:)];
    else
      %%% There is a bug, it doesn't change
      %Xt0 = [Xt(daysshift+1:end,:); repmat([Xt(end,:)],daysshift,1)]; %Change Xt->Xt0, status now starts from T0(120d)
      Xt0 = Xt(daysshift+1:end,:);
      
      
      %Xt0 = Xt(daysshift+times+1,:);             %Replace Xt->Xt0, status starts from T0(120d)
    end
    T_rel = (0:1:SamplingLastDay)+1; %%The relative days after 120d 
    for a = 1:par.maxa
        serodist = gen_strain_titres(Xt0, T_rel, a, par); % [times x titre]
        serodist(find(serodist< 0)) = 0; %replace negative values by zero
        multinomial(a).p = serodist./(diag(sum(serodist,2))*ones(size(serodist)));
        multi_p_byagebytime(a,:,:) = multinomial(a).p;
    end

    %calculate loglikelihood
    %yt_array first element is 120d
    AbCut = 3;
    for a=1:par.maxa %for each age group
        prob_x = multinomial(a).p;
        prob = prob_x * par.matY;
        if par.maxi == 10
            prob_corr = [sum(prob(:,1:AbCut),2) sum(prob(:,AbCut+1:end),2)];
        else
            prob_corr = prob;
        end
        obs_titres = yt_array(a).obs_titres;
        llh_array = log(prob_corr).*obs_titres;
        llh_array(isinf(llh_array)) = -10E2; %If -Inf, replace to -10E2 
        llh_array(~isfinite(llh_array))=0;  %If NaN, replace to zero
        LLH_Map_age(a).llh = sum(sum(llh_array));
    end
    
    LLH_Map_totalage.llh = sum([LLH_Map_age.llh]);
    nLLH = -LLH_Map_totalage.llh;
end
