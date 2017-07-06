function [TSubset] = extract_titres_table_paired(filename, data_time, data_strain, data_title, data_age )
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
[ titres_sheet metadata] = extract_titres_exvac(filename);
TitresTable = cell2table(titres_sheet);
for i=1:length(metadata(:,1))
	newStr = regexprep(metadata(i,1),'\.','_');
	TitresTable.Properties.VariableNames(i) = newStr;
end

%create tmp subset
kset = 1:2;
tmpTitresSubset = TitresTable(:,['sr_index',data_time(kset),data_strain(kset),data_age(kset)]);

%filter out unwanted data
removeid = strcmp(table2cell(tmpTitresSubset(:,data_strain(kset))),'-1'); %Remove records whose titres = -1
removeid = removeid(:,1) | removeid(:,2);
tmpTitresSubset(removeid,:) = [];
removeid = strcmp(table2cell(tmpTitresSubset(:,data_age(kset))),'-1'); %Remove records whose titres = -1
removeid = removeid(:,1) | removeid(:,2);
tmpTitresSubset(removeid,:) = [];

%create subset
TSubset = table();
T1 = table(str2num(char(table2cell(tmpTitresSubset(:,'sr_index')))));
T1.Properties.VariableNames{'Var1'} = 'sr_index';
TSubset = T1;

T1 = table(char(table2cell(tmpTitresSubset(:,data_time(1)))));
T1.Properties.VariableNames{'Var1'} = 'blood_1_date';
TSubset = [TSubset T1];
T1 = table(char(table2cell(tmpTitresSubset(:,data_time(2)))));
T1.Properties.VariableNames{'Var1'} = 'blood_2_date';
TSubset = [TSubset T1];
T1 = table(str2num(char(table2cell(tmpTitresSubset(:,data_age(1))))));
T1.Properties.VariableNames{'Var1'} = 'Age_1';
TSubset = [TSubset T1];
T1 = table(str2num(char(table2cell(tmpTitresSubset(:,data_age(2))))));
T1.Properties.VariableNames{'Var1'} = 'Age_2';
TSubset = [TSubset T1];
T1 = table(str2num(char(table2cell(tmpTitresSubset(:,data_strain(1))))));
T1.Properties.VariableNames{'Var1'} = 'Titres_T1';
TSubset = [TSubset T1];
T1 = table(str2num(char(table2cell(tmpTitresSubset(:,data_strain(2))))));
T1.Properties.VariableNames{'Var1'} = 'Titres_T2';
TSubset = [TSubset T1];

titres = [];
for i=1:length(table2array(TSubset(:,'Titres_T1')))
   titres(i) = round(log2(table2array(TSubset(i,'Titres_T1')))-log2(5));
end
T1 = table(titres');
T1.Properties.VariableNames{'Var1'} = 'Levels_T1';
TSubset = [TSubset T1];

titres = [];
for i=1:length(table2array(TSubset(:,'Titres_T2')))
   titres(i) = round(log2(table2array(TSubset(i,'Titres_T2')))-log2(5));
end
T1 = table(titres');
T1.Properties.VariableNames{'Var1'} = 'Levels_T2';
TSubset = [TSubset T1];
end