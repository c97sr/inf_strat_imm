function [ observe ] = getObsSeroprevalence( par )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    Antibody = par.Antibody;
    SampleSize = Antibody.samplesize;
    Abl = Antibody.Abl;
    age = Antibody.age;
    Abl(find(Abl>9)) = 9; %substitute Ab level >9

    
    %Only 0 and 1 immune status
    %par.maxt = 1;
    Abl(find(Abl<3)) = 0;
    Abl(find(Abl>=3)) = 1;

    Antibody = par.Antibody;
    OutbreakStartingDay = par.OutbreakStartingDay; %default 120 days 
    corrected_numdays = Antibody.numdays-OutbreakStartingDay;
    
%Transform observed titres Abl into 2D-array [time x titres]
for a=1:par.maxa
    obs_titres = zeros(par.SamplingLastDay+1, 2);
for i=1:length(Abl)
    titres = Abl(i);
    time = corrected_numdays(i);
    ind_age = age(i);
    if ind_age>=par.ages(a,1) & ind_age<par.ages(a,2)
        if time<par.SamplingLastDay+1
            obs_titres(time,titres+1) = obs_titres(time,titres+1)+1;
        end
    end
end
    observe(a).obs_titres = obs_titres;
    observe(a).obs_titres_numdays = OutbreakStartingDay+(0:par.SamplingLastDay)'; % from 1 Jan 2009
end
obs_titres = observe(1).obs_titres + observe(2).obs_titres + observe(3).obs_titres + observe(4).obs_titres;

end

