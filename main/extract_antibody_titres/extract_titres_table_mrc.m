function [TSubset] = extract_titres_table_mrc(data_time, data_strain, data_title, data_age, data_subtype, data_titre, data_sampleid  )
% main_dynamics_hk Summary of this function goes here
% plot time series population herd immunity
% serological data sampled from HongKong
% Written by Sean Yuan (hyuan@imperial.ac.uk)

%titres_level(:,1): Antibody Dilutions
%titres_level(:,2): Antibody levels
%titres_level(:,3): sr.index
%titres_level(:,4): Age

%setup parameters
%pars = InitParameters();
%pars.alpha = 0/365; %immune decay

%retrieve population data
[ titres_sheet metadata] = extract_titres_exvac_mrc('');
TitresTable = cell2table(titres_sheet);
for i=1:length(metadata(:,1))
	newStr = regexprep(metadata(i,1),'\.','_');
	if strcmp(newStr{1},'') 
      newStr = {'Index'};
    end
    TitresTable.Properties.VariableNames(i) = newStr;
end

%create tmp subset
tmpTitresSubset = TitresTable(:,{data_time,data_strain,data_age,data_sampleid,'Subtype','Titre','Month'});

%filter out unwanted data
%Select only one strain
removeid = strcmp(table2cell(tmpTitresSubset(:,data_strain)),'-1'); %Remove records whose titres = -1
tmpTitresSubset(removeid,:) = [];
removeid = strcmp(table2cell(tmpTitresSubset(:,data_age)),'-1'); %Remove records whose titres = -1
tmpTitresSubset(removeid,:) = [];

%create subset
TSubset = table();
%T1 = table(str2num(char(table2cell(tmpTitresSubset(:,'Index')))));
%T1.Properties.VariableNames{'Var1'} = 'Index';
%TSubset = [TSubset T1];

fulldate = {};
yyyy = tmpTitresSubset.Year;
mm = tmpTitresSubset.Month;
for i=1:height(tmpTitresSubset)

        m = mm{i};
        if length(m) < 2
            m = ['0' m];
        end
	fulldate(i,1) = {['15/' m '/' yyyy{i}]};
end

T1 = table(char(fulldate));
T1.Properties.VariableNames{'Var1'} = 'Date';
TSubset = [TSubset T1];
T1 = table(str2num(char(table2cell(tmpTitresSubset(:,data_age)))));
T1.Properties.VariableNames{'Var1'} = 'Age';
TSubset = [TSubset T1];
T1 = table(char(table2cell(tmpTitresSubset(:,data_strain))));
T1.Properties.VariableNames{'Var1'} = 'Strain';
TSubset = [TSubset T1];
T1 = table(char(table2cell(tmpTitresSubset(:,'Subtype'))));
T1.Properties.VariableNames{'Var1'} = 'Subtype';
TSubset = [TSubset T1];

T1 = table(char(table2cell(tmpTitresSubset(:,'Titre'))));
T1.Properties.VariableNames{'Var1'} = 'Titres';
TSubset = [TSubset T1];

T1 = table(char(table2cell(tmpTitresSubset(:,data_sampleid))));
T1.Properties.VariableNames{'Var1'} = 'Sample_ID';
TSubset = [TSubset T1];

end