function [Params] = setParameters(pa,par_name,par_value)
%Read and initialize Parameters
% setParameters Summary of this function goes here
% Function to update selected parameters used in the model
% Written by Sean Yuan (hyuan@imperial.ac.uk) 

%update the parameters values
if strcmp(par_name,'R0')
  pa.R0 = par_value;
  disp('The R0 value has been reset. Make sure whether it has conflict with transmision rate beta.');
end
if strcmp(par_name,'beta')
  pa.beta = par_value;
end
if strcmp(par_name,'AbB')
  pa.AbB = par_value;
  pa.semi_mechanistic = 0;
end

if strcmp(par_name,'model')
  pa.model = par_value;
end

if strcmp(par_name,'seed')
  pa.seed = [par_value];
end

if strcmp(par_name,'wan')
  pa.wan = [par_value];
end



if strcmp(par_name,'NormBoosting_flag')
  pa.NormBoosting_flag = par_value;
end
if strcmp(par_name,'ProtectionSlope_flag')
  pa.ProtectionSlope_flag = par_value;
end

if strcmp(par_name,'OutbreakNDA')
  pa.OutbreakNDA = par_value;
end

if strcmp(par_name,'age_mix_flag')
  pa.age_mix_flag = [par_value];
end

%--------------------------------
if strcmp(par_name,'age_flag')
  pa.age_flag = [par_value];
end
if strcmp(par_name,'AbB1')
  pa.AbB1 = par_value;
end
if strcmp(par_name,'AbB2')
  pa.AbB2 = par_value;
end
if strcmp(par_name,'AbB3')
  pa.AbB3 = par_value;
end
if strcmp(par_name,'AbB4')
  pa.AbB4 = par_value;
end
%--------------------------------
if strcmp(par_name,'frac_flag')
  pa.frac_flag = [par_value];
  if pa.frac_flag == 1
    pa.assort_flag = 0;
  end
end
if strcmp(par_name,'inittitres_flag')
  pa.inittitres_flag = par_value;
end
if strcmp(par_name,'inittitres')
  pa.inittitres = par_value;
end
if strcmp(par_name,'ContFrac1')
  pa.ContFrac1 = par_value;
end
if strcmp(par_name,'ContFrac2')
  pa.ContFrac2 = par_value;
end
if strcmp(par_name,'ContFrac3')
  pa.ContFrac3 = par_value;
end
if strcmp(par_name,'ContFrac4')
  pa.ContFrac4 = par_value;
end
%--------------------------------
if strcmp(par_name,'assort_flag')
  pa.assort_flag = [par_value];
  if pa.assort_flag == 1
    pa.frac_flag = 0;
  end
end
if strcmp(par_name,'ContAssort1')
  pa.ContAssort1 = par_value;
end
if strcmp(par_name,'ContAssort2')
  pa.ContAssort2 = par_value;
end
if strcmp(par_name,'ContAssort3')
  pa.ContAssort3 = par_value;
end
if strcmp(par_name,'ContAssort4')
  pa.ContAssort4 = par_value;
end


%--------------------------------
if strcmp(par_name,'immune_flag')
  pa.immune_flag = [par_value];
end
if strcmp(par_name,'immune_alpha')
  pa.immune_alpha = par_value;
end
if strcmp(par_name,'immune_alpha1')
  pa.immune_alpha1 = par_value;
end
if strcmp(par_name,'immune_alpha2')
  pa.immune_alpha2 = par_value;
end
if strcmp(par_name,'immune_alpha3')
  pa.immune_alpha3 = par_value;
end
if strcmp(par_name,'immune_alpha4')
  pa.immune_alpha4 = par_value;
end

if strcmp(par_name,'immune_beta')
  pa.immune_beta = par_value;
end
if strcmp(par_name,'immes')
  pa.semi_mechanistic = 1;
  pa.sa = par_value;
end
if strcmp(par_name,'alpha')
  pa.alpha = par_value;
end
if strcmp(par_name,'Tg')
  pa.Tg = par_value;
end
if strcmp(par_name,'currK')
  pa.currK = par_value;
end
if strcmp(par_name,'prevK')
  pa.prevK = [par_value];
end
if strcmp(par_name,'maxi')
  pa.maxi = par_value;
  %if update maxi, recalculate S I R lookups
  [Slu Ilu Rlu CIlu ] = make_S_I_R_lookups();
  pa.arrSlu = Slu;
  pa.arrIlu = Ilu;
  pa.arrRlu = Rlu;
  pa.arrCIlu = CIlu;
  pa.novars = max(max(max(pa.arrCIlu)));
end
if strcmp(par_name,'mint')
  pa.mint = par_value;
end
if strcmp(par_name,'maxt')
  pa.maxt = par_value;
end
if strcmp(par_name,'maxa')
  pa.maxa = par_value;
end
if strcmp(par_name,'PUAb')
  pa.PUAb = par_value;
end
if strcmp(par_name,'s0_imm')
  pa.s0_imm = [par_value];
end
if strcmp(par_name,'errp')
  pa.errp = [par_value];
end
if strcmp(par_name,'age_arr')
  pa.age_arr = [par_value];
end
if strcmp(par_name,'initS')
  pa.initS = [par_value];
end
if strcmp(par_name,'sera_sample_end')
  pa.sera_sample_end = [par_value];
end
if strcmp(par_name,'R0_list')
  if  ~isfield(pa, 'R0_list')
      pa.R0_list = [par_value];
  else
      pa.R0_list = [pa.R0_list par_value];
  end
end
if strcmp(par_name,'AbB_list')
  if  ~isfield(pa, 'AbB_list')
      pa.AbB_list = [par_value];
  else
      pa.AbB_list = [pa.AbB_list par_value];
  end
end

if pa.age_flag == 0
  pa.AbB1 = pa.AbB;
  pa.AbB2 = pa.AbB;
  pa.AbB3 = pa.AbB;
  pa.AbB4 = pa.AbB; 
end

if pa.immune_flag == 2
  pt = pa.immune_alpha1;
  pa.immune_alpha2 = pt;
  pa.immune_alpha3 = pt;
end

if pa.immune_flag == 0
  pt = pa.immune_alpha;
  pa.immune_alpha1 = pt;
  pa.immune_alpha2 = pt;
  pa.immune_alpha3 = pt;
  pa.immune_alpha4 = pt;
end

pa.matM = make_M(pa); 
pa.arrf = make_f_simple(pa);
pa.arrg = make_g(pa);
pa.arrh = make_h(pa);
pa.s0_imm = make_s0_imm(pa);

Params = pa;

function [ S I CI ] = make_S_I_lookups( input_args )
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
