function [ M ] = make_M(pa)
%make_M_simple # Function to return the most simple immune state boost:
% Simple function to return a null mixing matrix
% Written by Sean Yuan (hyuan@imperial.ac.uk) 
maxa = pa.maxa; %number of age groups
M = zeros(maxa,maxa);
a1=1;
b1=1;
a2=1;
b2=1;
a3=1;
b3=1;
a4=1;
b4=1;

if pa.frac_flag == 1
    if isfield(pa, 'ContFrac1')
        a1= pa.ContFrac1;
        b1= pa.ContFrac1;
    end
    if isfield(pa, 'ContFrac2')
        a2= pa.ContFrac2;
        b2= pa.ContFrac2;
    end
    if isfield(pa, 'ContFrac3')
        a3= pa.ContFrac3;
        b3= pa.ContFrac3;
    end
    if isfield(pa, 'ContFrac4')
        a4= pa.ContFrac4;
        b4= pa.ContFrac4;
    end
end

if pa.assort_flag == 1
    if isfield(pa, 'ContAssort1')
        a1 = pa.ContAssort1;
        %b1 = 1;
        b1= 1/pa.ContAssort1;
    end
    if isfield(pa, 'ContAssort2')
        a2 = pa.ContAssort2;
        %b2 = 1;
        b2= 1/pa.ContAssort2;
    end
    if isfield(pa, 'ContAssort3')
        a3 = pa.ContAssort3;
        %b3 = 1;
        b3= 1/pa.ContAssort3;
    end
    if isfield(pa, 'ContAssort4')
        a4 = pa.ContAssort4;
        %b4 = 1;
        b4= 1/pa.ContAssort4;
    end
end

M_ori = [
  7.96   2.11    2.10    0.35; 
  2.09   6.56    3.97    2.00;
  2.41   4.68    5.74    2.74;
  0.16   1.04    1.43    2.22;   
 ];



Cont = [
 a1 b2 b3 b4
 b1 a2 b3 b4
 b1 b2 a3 b4
 b1 b2 b3 a4
];

M = M_ori.*Cont;
if pa.age_mix_flag == 0
  M = ones(pa.maxa, pa.maxa);
end
%if pa.age_flag == 0 %no age specific boosting but still has age mixing
%effect
%  M = ones(pa.maxa, pa.maxa);
%end

if pa.model == 7
    M = ones(pa.maxa, pa.maxa);
end
%M_ori = [
%  7.96   2.11    2.10    0.35; 
%  2.09   6.56    3.97    2.00;
%  2.41   4.68    5.74    2.74;
%  0.16   1.04    1.43    2.22;   
% ];

if pa.model == 7.1
totalCnt = sum(M);
cnt1= totalCnt(1);
cnt2= totalCnt(2);
cnt3= totalCnt(3);
cnt4= totalCnt(4);

b4=M(4,4)*a4/(cnt4-M(4,4)); %max.a4 = 2.293
b3=M(3,3)*a3/(cnt3-M(3,3)); %max.a3 = 1.307
b2=M(2,2)*a2/(cnt2-M(2,2)); %max.a2 = 1.194
b1=M(1,1)*a1/(cnt1-M(1,1)); %max.a1 = 0.527 

%M = [
%  7.96*(1+a1)   2.11*(1-b2)    2.10*(1-b3)    0.35*(1-b4); 
%  2.09*(1-b1)   6.56*(1+a2)    3.97*(1-b3)    2.00*(1-b4);
%  2.41*(1-b1)   4.68*(1-b2)    5.74*(1+a3)    2.74*(1-b4);
%  0.16*(1-b1)   1.04*(1-b2)    1.43*(1-b3)    2.22*(1+a4);   
%];
end