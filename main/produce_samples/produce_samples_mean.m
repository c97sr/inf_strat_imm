function [sim_samples] = produce_samples_mean( PosteriorSamples, pars)
% Produce simulated samples at T1 and T2 using posterior mean
% Written by Sean Yuan (hyuan@imperial.ac.uk) 
% Plot each age group information for only 10 titres levels
% E:\working\Projects.IC\Projects\isl\mat\Misltr\isltr-1.4\out\mcmc\ph1n1\20150106

global proj Antibody;
init_collect = 1;
second_collect = 2;
third_collect = 3;
k = 2;
Ab = Antibody;
posterior = table2array(PosteriorSamples);
meanposterior = mean(posterior);
%[maxpost idx] = max(posterior(:,end));
idx = 1;
load('sampling_time');
sampling_num = zeros(366,1);
for a=1:pars.maxa
  age(a).sampling_num = sum(observe(a).obs_titres,2);
end

for i = 1:length(idx)
    vars = PosteriorSamples.Properties.VariableNames;
    for p=1:length(vars)
        if strcmpi('LLH',vars(p))
        else
           [pars] = setParameters(pars,char(vars(p)),meanposterior(idx(i),p));
        end
    end
    
    %set parameters
    beta = pars.beta;
    AbB = [pars.AbB1 pars.AbB2 pars.AbB3 pars.AbB4];
    immune_alpha = [pars.immune_alpha1 pars.immune_alpha2 pars.immune_alpha3 pars.immune_alpha4];
    lastsamplingday = pars.SamplingLastDay;

    %setup initial condition
    [yini age_arr] = make_ics_naive( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu);
    %[yini_k1 age_arr_k1] = make_ics_fromtitres_byage( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, Ab.K(init_collect).Abl, Ab.K(init_collect).age);
    %[yini_k2 age_arr_k2] = make_ics_fromtitres_byage( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, Ab.K(k).Abl, Ab.K(k).age);
    %[yini_k3 age_arr_k3] = make_ics_fromtitres_byage( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, Ab.K(k).Abl, Ab.K(k).age);


    %setep simulation time
    T0 = pars.OutbreakStartingDay;
    meanKdays(1) = mean(pars.Antibody.K(1).numdays - T0);
    meanKdays(2) = mean(pars.Antibody.K(2).numdays - T0);
    %meanKdays(3) = mean(pars.Antibody.K(3).numdays - T0);
    sample_time_K1 = round(meanKdays(1));
    sample_time_K2 = round(meanKdays(2));
    times = 0:1:lastsamplingday;
    sample_size_K1 = Ab.K(1).samplesize;
    sample_size_K2 = Ab.K(2).samplesize;


    
    % Update parameters
    %for i=1:length(theta)
    %    par = setParameters(par,theta_name{i},theta(i));
    %end
    
    %javaaddpath e:\workspace\MyJavaProject\bin\matlabjava.jar; %use jave working dir 
    javaaddpath(pars.javapath); %set ./java as default dir
    mepar_3b = matlabjava.Parameters;
    meser_3b = matlabjava.Serology;
    meser_3b.setParameters(mepar_3b);
    % set parameters
    meser_3b.updateParametersG(pars.arrg);
    meser_3b.updateParametersH(pars.arrh);
    meser_3b.updateParametersM(pars.matM);
    meser_3b.updateParametersBeta(pars.beta);  
    meser_3b.updateParameters('wan',pars.wan);
    meser_3b.updateParameters('s0_imm', pars.s0_imm);
    x0 = yini;  
    [t y] = ode23(@(t,x)odef_islmodjava(t,x, meser_3b), times, x0);  
    %[t y] = ode23(@(t,x)odef_islmod(t,x,pars), times, x0);
    clear('mepar_3b');
    clear('meser_3b');
    T = t;

    
    for a=1:4
        Xout_k1(:,:) = retrieve_histogram(y, pars, times, sample_time_K1, a); % model output
        %Xout_k1_2(:,:) = retrieve_histogram_20151001(y, pars, times, sample_time_K1, a); % model output
        Xoutput_k1_list(a,i,:) = Xout_k1;
        Xout_k2 = retrieve_histogram(y, pars, times, sample_time_K2, a); % model output
        Xoutput_k2_list(a,i,:) = Xout_k2;
    end
    
    for a=1:4
        for sample_time=1:366
            sample_time
            Xout_t(:,:) = retrieve_histogram(y, pars, times, sample_time, a); % model output
            Xoutput_list(a).p(sample_time,:) = Xout_t;
        end
    end
end

for a=1:4
    for day=1:366
       n = age(a).sampling_num(day);
       p = Xoutput_list(a).p(day,:);
       if n>1
            disp n;
       end
       sim_samples(a).obs_titres(day,:) = mnrnd(n,p,1);
    end
    sim_samples(a).obs_titres_numdays = (120:120+366-1)';
end


%%save age(a).R as ouput
obs_titres = sim_samples(1).obs_titres+sim_samples(2).obs_titres+sim_samples(3).obs_titres+sim_samples(4).obs_titres
obs_titres(39:41,:) = 0;
Abl_1 = 0;
for d = 1:180
  for j=1:10  % titre
     if obs_titres(d,j)>0
        for i=1:obs_titres(d,j)
          Abl_1(end+1,1) = j-1;
        end
     end
  end
end

Abl_2 = 0;
for d = 181:366
  for j=1:10  % titre
     if obs_titres(d,j)>0
        for i=1:obs_titres(d,j)
          Abl_2(end+1,1) = j-1;
        end
     end
  end
end

        

end