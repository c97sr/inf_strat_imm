function [nLLH LLH_Map_age]= getNegLLHAgejava(theta, theta_name, par, meser, x0)
%x: input parameter values
%x: input parameter names
%par: parameter object
%meser: sera model in java
%yini: initial data
    

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
    Abl(find(Abl>par.maxt)) = par.maxt; %substitute Ab level >8 by 8;set max antibody level=8
    %SampleSize = length(Abl);
    %SampleTime = unique(numdays) - OutbreakStartingDay;                     

    
    
    %%
    % run ODE
    times = 0:1:365*1.2; %Make sure they reach equilibrium
    % set parameters
    meser.updateParametersG(par.arrg);
    meser.updateParametersH(par.arrh);
    meser.updateParametersM(par.matM);
    meser.updateParametersBeta(par.beta); 
    %meser.updateParameters('wan',par.wan);
    [t Xt] = ode23(@(t,x)odef_islmodjava(t,x, meser), times, x0);
    %[t Xt] = ode23(@(t,x)odef_islmod(t,x, par), 1:100, x0);
    %for a=1:par.maxa
    %    Xt(:,(a-1)*10+1) = Xt(:,(a-1)*10+1)+par.s0_imm(a);
    %end   
    %some elements have negative values, will cause errors
    
    % define serology distribution matrix
    %serodist_mat = zeros(par.maxa,t(end),par.maxi); % [age x times x titre]
    for a = 1:par.maxa
        serodist = gen_strain_titres(Xt, times, a, par); % [times x titre]
        serodist(find(serodist< 0)) = 0; %replace negative values by zero
        multinomial(a).p = serodist./(diag(sum(serodist,2))*ones(size(serodist)));
        multi_p_byagebytime(a,:,:) = multinomial(a).p;
        %multi_p_byagebytime(a,:,:) = serodist./(diag(sum(serodist,2))*ones(size(serodist)));
        %multinomial(a).p = serodist./(diag(sum(serodist,2))*ones(size(serodist)));
    end

    

    %%
    %calculate loglikelihood
    for a=1:par.maxa %for each age group
        Yt = Abl(find(age>=par.ages(a,1) & age<par.ages(a,2)));
        sampletime = numdays(find(age>=par.ages(a,1) & age<par.ages(a,2))) - OutbreakStartingDay;
        llh = 0;
        Yt = Yt(find(sampletime < SamplingLastDay));
        sampletime = sampletime(find(sampletime < SamplingLastDay));
        %analize individual sequentlially
        for i=1:length(Yt)
            abl_obs = Yt(i);
            time_obs = sampletime(i); 
            if time_obs > times(end)
                time_obs = times(end);
            end
            prob = multinomial(a).p(time_obs,:);
            llh = llh + calculateLogLikelihood(prob,abl_obs+1);
        end
        LLH_Map_age(a).llh = llh;
    end
    
    LLH_Map_totalage.llh = sum([LLH_Map_age.llh]);
    nLLH = -LLH_Map_totalage.llh
end
