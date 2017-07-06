function [ Output ] = main_producetitres_bytime(strain)
% main_producetitres_bytime Summary of this function goes here
% serological data sampled from HongKong
% extract titres data from csv file and save into mat variables
% Written by Sean Yuan (hyuan@imperial.ac.uk)

column_date = {'blood_1_date','blood_2_date','blood_3_date','blood_4_date'};
if strcmp(strain,'H3N2')
    column_strain = {'H3N2.T1','H3N2.T2','H3N2.T3','H3N2.T4'};
end
if strcmp(strain,'H1N1')
    %column_strain = {'H1N1.T1','H1N1.T2','H1N1.T3','H1N1.T4'};
    column_strain = {'H1N1_T1','H1N1_T2','H1N1_T3','H1N1_T4'};
end
column_title = {'T1:07/2009-09/2009','T2:11/2009-01/2010','T3:12/2010-03/2011','08/2011-12/2011'};
column_age = {'Age_1','Age_2','Age_3','Age_4'};

filename = './dat/part_R14_2013_07_22_sr_vac.csv';
TitresTablePaired = extract_titres_table_paired(filename, column_date, column_strain, column_title, column_age); 



% plot antibody titres
titres_list = {};
age_list = {};
TitresTableList = [];
for i=1:length(column_date)
    %%%VVV No such function!!!!
    T = extract_titres_table(column_date{i}, column_strain{i}, column_title{i}, column_age{i}); 
    T1 = table(ones(height(T),1)*i);
    T1.Properties.VariableNames{'Var1'} = 'SamplingK'; %add one column SamplingK
    T = [T T1];
    TitresTableList(i).T = T;  
    %no such function
    %[titres ages] = plot_histogram_hk_new(column_date{i}, column_strain{i}, column_title{i}, column_age{i});    
    %titres_list(i) = {titres}; 
    %age_list(i) = {ages}; 
end
TitresTableTotal = [TitresTableList(1).T; TitresTableList(2).T; TitresTableList(3).T; TitresTableList(4).T];

%add one column number of days from 01/01/2009
nodays = [];
for i=1:height(TitresTableTotal)
   enddate = table2array(TitresTableTotal(i,'blood_date'));
   startdate = '01/01/2009';
   nodays(i) = datenum(enddate,'dd/mm/yyyy')-datenum(startdate,'dd/mm/yyyy');
end
T1 = table(nodays');
T1.Properties.VariableNames{'Var1'} = 'Nodays'; 
TitresTableTotal = [TitresTableTotal T1];


titres_total = TitresTableTotal.Levels;
age_total = TitresTableTotal.Age;
blood_date = TitresTableTotal.blood_date;
numdays = TitresTableTotal.Nodays;
samplingK = TitresTableTotal.SamplingK;
sr_index = TitresTableTotal.sr_index;

% save antibody titres
% Attributes: samplesize, Abl, age, blood_date, nodays
Antibody = [];
Antibody.samplesize = length(titres_total);
Antibody.Abl = titres_total;
Antibody.age = age_total;
Antibody.date = blood_date;
Antibody.numdays = numdays;
Antibody.samplingK = samplingK;
Antibody.sr_index = sr_index;

for k=1:length(column_strain)
    Antibody.K(k).Abl = Antibody.Abl(find(Antibody.samplingK == k));
    Antibody.K(k).age = Antibody.age(find(Antibody.samplingK == k));
    Antibody.K(k).date = Antibody.date(find(Antibody.samplingK == k),:);
    Antibody.K(k).numdays = Antibody.numdays(find(Antibody.samplingK == k));
    Antibody.K(k).sr_index = Antibody.sr_index(find(Antibody.samplingK == k));
    Antibody.K(k).samplesize = length(Antibody.K(k).Abl); 
end

if strcmp(strain,'H1N1')
    pars.filename = 'h1n1_titres.mat';
    pars.proj = 'hk_ph1n1';
end
if strcmp(strain,'H3N2')
    pars.filename = 'h3n2_titres.mat';
    pars.proj = 'hk_h3n2';
end

% create output directory
date_str = [datestr(now,10) datestr(now,5) datestr(now,7)];
proj_str = pars.proj;
pars.out_dir = ['out/' date_str '/' proj_str];
%pars.out_dir = ['out/' datestr(now,10) datestr(now,5) datestr(now,7) '_' pars.proj];
out_dir = pars.out_dir
if(exist(out_dir)==7)
else
    mkdir(out_dir)
end
    
params = pars;
save([params.out_dir './' params.filename],'Antibody','TitresTableTotal','TitresTablePaired','params');
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    

