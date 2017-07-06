function [ Output ] = main_producetitres_bytime_mrc(fname)
% main_producetitres_bytime Summary of this function goes here
% serological data sampled from HongKong
% extract titres data from csv file and save into mat variables
% Written by Sean Yuan (hyuan@imperial.ac.uk)

%Date
%Strains
%SamplingPeriod
%Ages



column_date = {'Year'};
column_month = {'Month'};
%if strcmp(strain,'H3N2')
column_subtype = {'Subtype'};
column_strain = {'Strain'};

strain_list ={'A/England/878/1969','A/Sydney/5/1997','A/Fujian/411/2002','A/USSR/90/1977','A/Beijing/262/1995','A/NewCaledonia/20/1999'};
prompt = ['which strain you want to extract: \n' ...
         'H3N2 strain [1]: A/England/878/1969 \n' ...
         'H3N2 strain [2]: A/Sydney/5/1997 \n' ...
         'H3N2 strain [3]: A/Fujian/411/2002 \n' ...
         'H1N1 strain [4]: A/USSR/90/1977 \n' ...
         'H1N1 strain [5]: A/Beijing/262/1995 \n' ...
         'H1N1 strain [6]: A/NewCaledonia/20/1999 \n'];

in_x = input(prompt);
tar_strain = in_x; %strain we want to extract
if in_x <= 3
  strain = 'H3N2';
else
  strin = 'H1N1';
end
disp(['strain ' strain_list{in_x} ' was selected. serological data is extracted...']);

%H3N2
%strain 1: 'A/England/878/1969'
%strain 2: 'A/Sydney/5/1997'
%strain 3: 'A/Fujian/411/2002'
%H1N1
%strain 4: 'A/USSR/90/1977'
%strain 5: 'A/Beijing/262/1995'
%strain 6: 'A/NewCaledonia/20/1999'


column_title = {'Sample_ID'};

column_age = {'Age'};

column_titre = {'Titres'};

column_sample_id = {'Sample_ID'};

filename = '../../dat/20160914/serological.csv';
filename = fname;

% plot antibody titres
titres_list = {};
age_list = {};
TitresTableList = [];
%for i=1:length(column_date)
    %%%VVV No such function!!!!
    i = 1;
    T = extract_titres_table_mrc(column_date{i}, column_strain{i}, column_title{i}, column_age{i}, column_subtype{i}, column_titre{i}, column_sample_id{i}); 
    Tcell = cellstr(T.Strain);
    str = unique(cellstr(T.Strain));
    Tmp_C1 = table(zeros(height(T),1)*1);
    Tmp_C1.Properties.VariableNames{'Var1'} = 'StrainJ';
    T = [T Tmp_C1];
    Tsub = cellstr(T.Subtype);
    %str = unique(cellstr(T.Subtype));
    index = find(ismember(Tsub,strain));
    %T = T(index,:);
    %for i=1:length(str)
             
            index = find(ismember(Tcell,'A/England/878/1969'));  %H3n2
            T(index,'StrainJ') = {1};
            index = find(ismember(Tcell,'A/Sydney/5/1997'));     %H3N2
            T(index,'StrainJ') = {2};
            index = find(ismember(Tcell,'A/Fujian/411/2002'));   %H3N2
            T(index,'StrainJ') = {3};
            index = find(ismember(Tcell,'A/USSR/90/1977'));      %H1N1
            T(index,'StrainJ') = {4};
            index = find(ismember(Tcell,'A/Beijing/262/1995'));  %H1N1
            T(index,'StrainJ') = {5};
            index = find(ismember(Tcell,'A/NewCaledonia/20/1999')); %H1N1 
            T(index,'StrainJ') = {6};
            
    %end
    
    T1 = table(ones(height(T),1)*1);
    T1.Properties.VariableNames{'Var1'} = 'SamplingK'; %add one column SamplingK
    T = [T T1];
   
nodays = [];
if tar_strain == 2
    startdate = '01/01/1997'; %same year as Syndey 1997 (Because there is no 1996 data)
end
if tar_strain == 3
    startdate = '01/01/2001'; %one year earlier than Fujian 2002
end

%% Need to create the data structure consistent to pH1N1 from here
% SamplingK
% numdays: not sure
% Data: not sure

%TitresTableTotal = [TitresTableList(1).T; TitresTableList(2).T; TitresTableList(3).T; TitresTableList(4).T];

%add one column number of days from 01/01/2009
%%% Starting from here  03/10/2016


%filtering out target strain
T(find(T.StrainJ ~= tar_strain),:) = [];

%find the maximum SamplingK
maxnod = max(datenum(table2array(T(:,'Date')),'dd/mm/yyyy') - datenum(startdate,'dd/mm/yyyy'));
maxK = 1;
for k = 1:10
    if (maxnod > (k-1)*365 & maxnod < k*365)
      maxK = k;  
    end
end

for i=1:height(T)
   T(i,'SamplingK') = table(0);
   enddate = table2array(T(i,'Date'));   
   nod = datenum(enddate,'dd/mm/yyyy')-datenum(startdate,'dd/mm/yyyy');
   nodays(i) = nod;
   for k = 1:maxK
      if nod > (k-1)*365 & nod < k*365
          T(i,'SamplingK') = table(k);
      end
   end
end

% date_k = 1:5;
T1 = table(nodays');
T1.Properties.VariableNames{'Var1'} = 'Nodays'; 
T = [T T1];
T(find(T.SamplingK < 1),:) = [];
TitresTableList(1).T = T;  
%% 2016/10/10


titres_total = str2num(TitresTableList(1).T.Titres);
titres_total(find(titres_total<0)) = 0;
age_total = TitresTableList(1).T.Age;
blood_date = char(TitresTableList(1).T.Date);
numdays = TitresTableList(1).T.Nodays;
samplingK = TitresTableList(1).T.SamplingK;
%sr_index = TitresTableTotal.sr_index;

% save antibody titres
% Attributes: samplesize, Abl, age, blood_date, nodays
Antibody = [];
Antibody.samplesize = length(titres_total);
Antibody.Abl = titres_total;
Antibody.age = age_total;
Antibody.date = blood_date;
Antibody.numdays = numdays;
Antibody.samplingK = samplingK;
%Antibody.sr_index = sr_index;



% K should change to strain J
for k=1:maxK
    Antibody.K(k).Abl = Antibody.Abl(find(Antibody.samplingK == k));
    Antibody.K(k).age = Antibody.age(find(Antibody.samplingK == k));
    Antibody.K(k).date = Antibody.date(find(Antibody.samplingK == k),:);
    Antibody.K(k).numdays = Antibody.numdays(find(Antibody.samplingK == k));
    %Antibody.K(k).sr_index = Antibody.sr_index(find(Antibody.samplingK == k));
    Antibody.K(k).samplesize = length(Antibody.K(k).Abl); 
end

if strcmp(strain,'H1N1')
    pars.filename = ['h1n1_titres_' num2str(in_x) '.mat'];
    pars.proj = 'world_h1n1';
end
if strcmp(strain,'H3N2')
    pars.filename = ['h3n2_titres_' num2str(in_x) '.mat'];
    pars.proj = 'world_h3n2';
end
if strcmp(strain,'B')
    pars.filename = 'b_titres.mat';
    pars.proj = 'world_b';
end
pars.strain = strain_list{tar_strain};

% create output directory
date_str = [datestr(now,10) datestr(now,5) datestr(now,7)];
proj_str = pars.proj;
pars.out_dir = ['out/' date_str '/' proj_str];
%pars.out_dir = ['out/' datestr(now,10) datestr(now,5) datestr(now,7) '_' pars.proj];
out_dir = pars.out_dir;
if(exist(out_dir)==7)
else
    mkdir(out_dir);
end
    
params = pars;
params.startdate = startdate;
TitresTableTotal = TitresTableList; 
save([params.out_dir '/' params.filename],'Antibody','TitresTableTotal','params');
disp(['save as ' params.out_dir '/' params.filename]);
clear all;
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    

