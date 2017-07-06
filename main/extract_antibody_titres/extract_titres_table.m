function [TSubset] = extract_titres_table(data_time, data_strain, data_title, data_age )
% main_dynamics_hk Summary of this function goes here
% plot time series population herd immunity
% serological data sampled from HongKong
% Written by Sean Yuan (hyuan@imperial.ac.uk)

%titres_level(:,1): Antibody Dilutions
%titres_level(:,2): Antibody levels
%titres_level(:,3): sr.index
%titres_level(:,4): Age

%setup parameters
pars = InitParameters();
pars.alpha = 0/365; %immune decay

%retrieve population data
[ titres_sheet metadata] = extract_titres_exvac;
TitresTable = cell2table(titres_sheet);
for i=1:length(metadata(:,1))
	newStr = regexprep(metadata(i,1),'\.','_');
	TitresTable.Properties.VariableNames(i) = newStr;
end

%create tmp subset
tmpTitresSubset = TitresTable(:,{data_time,data_strain,data_age,'sr_index'});

%filter out unwanted data
removeid = strcmp(table2cell(tmpTitresSubset(:,data_strain)),'-1'); %Remove records whose titres = -1
tmpTitresSubset(removeid,:) = [];
removeid = strcmp(table2cell(tmpTitresSubset(:,data_age)),'-1'); %Remove records whose titres = -1
tmpTitresSubset(removeid,:) = [];

%create subset
TSubset = table();
T1 = table(str2num(char(table2cell(tmpTitresSubset(:,'sr_index')))));
T1.Properties.VariableNames{'Var1'} = 'sr_index';
TSubset = [TSubset T1];
T1 = table(char(table2cell(tmpTitresSubset(:,data_time))));
T1.Properties.VariableNames{'Var1'} = 'blood_date';
TSubset = [TSubset T1];
T1 = table(str2num(char(table2cell(tmpTitresSubset(:,data_age)))));
T1.Properties.VariableNames{'Var1'} = 'Age';
TSubset = [TSubset T1];
T1 = table(str2num(char(table2cell(tmpTitresSubset(:,data_strain)))));
T1.Properties.VariableNames{'Var1'} = 'Titres';
TSubset = [TSubset T1];

titres = [];
for i=1:length(table2array(TSubset(:,'Titres')))
   titres(i) = round(log2(table2array(TSubset(i,'Titres')))-log2(5));
end
T1 = table(titres');
T1.Properties.VariableNames{'Var1'} = 'Levels';
TSubset = [TSubset T1];
end