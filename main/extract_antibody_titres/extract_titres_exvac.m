function [ titres_sheet metadata titres_sheet_exvac] = extract_titres_exvac( filename )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%fileID = fopen('./dat/all_R14_2013_07_22_sr.csv');
%fileID = fopen('./dat/part_R14_2013_07_22_sr.csv');
%fileID = fopen('./dat/part_R14_2013_07_22_sr_vac.csv');
fileID = fopen(filename);
if fileID == -1
    fileID = fopen('../../dat/part_R14_2013_07_22_sr_vac.csv');
end
vacid = 0;
count = 1;
while true
count;
tline = fgetl(fileID);
if ~ischar(tline)
    break;
end

if count == 1
    % first row contains metadata
    metadata = regexp(tline,',', 'split')';
    for i=1:length(metadata(:,1))
        term = metadata(i,1);
        indices = strcmp(metadata, term);
        metadata(i,2) = {find(indices == 1)};
    end
    %Metadata: column names of Viruses
    %indices = cellfun( @(x) strcmp(x, term),metadata);
    %colno = find(indices == 1);
    %metadata.name = colnames;
    %%%VVV still working on this
    cell2mat(metadata(strcmp(metadata, 'Old_id'),2))
    vacid = cell2mat(metadata(strcmp(metadata(:,1), 'vac0910_h1n1_cleaned'),2));
    count = count + 1;
else
    str = regexp(tline,',', 'split');
    newStr = regexprep(str,'NA','-1');
    if strcmp(newStr(vacid),'YES')
    else
    titres_sheet(count-1,:) = newStr;
    count = count + 1;
    end
    %%check newStr(17);
    %str  = textscan(fileID,'%s')
end
end

fclose(fileID);

