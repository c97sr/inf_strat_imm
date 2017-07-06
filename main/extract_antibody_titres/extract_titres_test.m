function [ titres_sheet metadata] = extract_titres_test( )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%fileID = fopen('./dat/all_R14_2013_07_22_sr.csv');
%fileID = fopen('./dat/part_R14_2013_07_22_sr.csv');
fileID = fopen('./dat/part_R14_2013_07_22_sr_vac.csv');

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
    strcmp(metadata, 'vac0910_h1n1_cleaned')
else
    str = regexp(tline,',', 'split');
    newStr = regexprep(str,'NA','-1');
    titres_sheet(count-1,:) = newStr;
    %str  = textscan(fileID,'%s')
end
count = count + 1;
end

fclose(fileID);

