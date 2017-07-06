function [Params] = InitParametersforSIR()
%Read and initialize Parameters
% InitParameters Summary of this function goes here
% Function to return all required parameters and lookup arrays in the model
% Written by Sean Yuan (hyuan@imperial.ac.uk) 

% Model path
pa.javapath = 'e:\workspace\MyJavaProject\bin\matlabjava.jar';
%pa.javapath = 'java/matlabjava.jar'

% Basic epidemiological parameters 
pa.beta = 0.5;
pa.Tg = 3.3;
%pa.Tg = 5;
pa.Tga = [3.3 3.3 3.3 3.3]; %age group1 
pa.wan = 1/(3*7);
pa.alpha = 0/420;
pa.mint = 0; %starting titres = 0 -> undetectable 
pa.maxt = 1;
pa.maxi = 2; %total titres levels for strainX
pa.maxj = 1; 
pa.maxk = 1;
pa.maxa = 1;  %number of age groups
pa.maxX = 1;  %max number of strains
pa.N = 0.7E+7;
pa.seed = 500;
pa.trickle = 0;

pa.PUAb = 0; % Protected ratio for undetectable Antibody titres group. 0.06 is nice value;
pa.s0_imm = zeros(1,pa.maxa);
pa.age_flag = false;
pa.semi_mechanistic = 0;
pa.NormBoosting_flag = 0;
pa.ProtectionSlope_flag = 0;

% Other population parameters
pa.demographic = [17.6 29.7 39.4 13.3] %proportion of 0-19; 20-39; 40-64; 65+
pa.ages = [0 20; 20 40; 40 65; 65 100]; % Same as 0-19, 20-39, 40-64, 65 100 if age is discrete
pa.startdate = '01/01/2009';
pa.firstcase = '01/05/2009'; %http://en.wikipedia.org/wiki/2009_flu_pandemic_in_Hong_Kong.151 days after 01/01/2009
pa.OutbreakStartingDay = datenum(pa.firstcase,'dd/mm/yyyy') - datenum(pa.startdate,'dd/mm/yyyy'); %disease starting day if first case on 11/06/2009
pa.SamplingLastDay = 365; %Don't simulate and sample for days after Lastday
pa.model = [];

% Antibody boosting level
pa.AbB = 4; 
pa.AbB1 = 4;
pa.AbB2 = 4;
pa.AbB3 = 4;
pa.AbB4 = 4;
% Contact proportion
pa.frac_flag = 0;
pa.ContFrac1 = 1;
pa.ContFrac2 = 1;
pa.ContFrac3 = 1;
pa.ContFrac4 = 1;
% Contact assortativity
pa.assort_flag = 0;
pa.ContAssort1 = 1;
pa.ContAssort2 = 1;
pa.ContAssort3 = 1;
pa.ContAssort4 = 1;
% Immune protection
pa.immune_flag = 0;
pa.immune_alpha = 3;
pa.immune_alpha1 = 3;
pa.immune_alpha2 = 3;
pa.immune_alpha3 = 3;
pa.immune_alpha4 = 3;
pa.immune_beta = 2.102; %Coudeville et al. Hobson data
%pa.immune_beta = 1.299; %Coudeville et al. All data 


[Slu Ilu Rlu CIlu ] = make_S_I_R_lookups();
pa.arrSlu = Slu;
pa.arrIlu = Ilu;
pa.arrRlu = Rlu;
pa.arrCIlu = CIlu;
pa.novars = max(max(max(pa.arrCIlu)));
Params = pa;


function [ f ] = make_f_simple( input_args )
%make_f_simple Summary of this function goes here
%   Simple function to return a null infectivity array of the correct dimensons
%   Written by Sean Yuan (hyuan@imperial.ac.uk) 
maxa = pa.maxa; %number of age groups
maxi = pa.maxi; %Ab level to strain A
maxj = pa.maxj; %Ab level to strain B
maxk = pa.maxk; %Ab level to strain C
f = zeros(maxa,maxi,maxj,maxk);
for a=1:maxa
    for i=1:maxi
        for j=1:maxj
            for k=1:maxk
                f(a,i,j,k)=1;
            end    
        end
    end
end
end



function [ S I R CI ] = make_S_I_R_lookups( input_args )
%make_S_I_simple Summary of this function goes here
% Function to make a lookup array for the linearized index of the isl model
% Can also be used to make list of column names for the solution table
% Written by Sean Yuan (hyuan@imperial.ac.uk) 
maxX = pa.maxX; %number of strains
maxa = pa.maxa; %number of age groups
maxi = pa.maxi; %Ab level to strain A
maxj = pa.maxj; %Ab level to strain B
maxk = pa.maxk; %Ab level to strain C
S = zeros(maxa,maxi,maxj,maxk);
I = zeros(maxX,maxa,maxi,maxj,maxk);
CI = zeros(maxX,maxa,maxi,maxj,maxk);
counter = 1;

%assign the number for susceptible
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                    S(a,i,j,k) = counter;
                    counter = counter + 1;
                end
            end
        end
    end


%assign the number for infected
for X=1:maxX
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                    I(X,a,i,j,k) = counter;
                    counter = counter + 1;
                end
            end
        end
    end
end

%assign the number for recovered
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                    R(a,i,j,k) = counter;
                    counter = counter + 1;
                end
            end
        end
    end


%assign the number for accumulated infected
for X=1:maxX
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                    CI(X,a,i,j,k) = counter;
                    counter = counter + 1;
                end
            end
        end
    end
end
end

end