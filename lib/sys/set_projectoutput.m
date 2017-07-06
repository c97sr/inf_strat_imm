function [ out_dir ] = set_projectoutput( mainoutdir, mainproj, datepath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if ~exist('datepath','var')
    datepath = [datestr(now,10) datestr(now,5) datestr(now,7)];
end
if ~strcmp(mainproj,'')
    out_dir = [mainoutdir '/' mainproj '/' datepath '/'];
else
    out_dir = [mainoutdir '/' datepath '/'];
end
%out_dir = proj;
if(exist(out_dir)==7)
else
   mkdir(out_dir)
end
end


function [] = checkParameter()
    %% Check parameters
if exist([ibmsdir 'dat/' metadata.ibms.filename])
  disp 'IBMS parameters are read.'
else
  disp 'Error: failed to read IBMS parameters.';
end
if exist(metadata.ibms.initFileVirus) 
  if exist(metadata.ibms.initFileSIR) 
    metadata.ibms.initFlag = 2; %IBM
  end
else
  if exist(metadata.ibms.initFileODE)   
    metadata.ibms.initFlag = 1; %ODE
  end 
end
end

