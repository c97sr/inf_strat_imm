function [ susc] = getSusc( titres2, mua, mub )
%getSusc Summary of this function goes here
%   input titres
%   return disease susceptibility
    %mua = 2.844; %From ALL
    %mub = 1.299; 
    if ~exist('mua') 
        mua = 3.385; %From Hobson
    end
    if ~exist('mub') 
        mub = 2.102;
    end
    a = mua;
    b = mub;
    %Titres = 2.^(titres2 + log2(5)); % dilution in excel sheet
    %susc = 1./(1+exp(b*(log(Titres)-a))); 20141029
    susc = 1./(1+exp(b*(titres2-a))); %a is 50% protection titres
    if titres2 == 0
      susc = 1; %set susc to be 1 one there is no antibody detected
    end
end

