function [ infecteds CInfecteds CInfecteds_age] = gen_total_infecteds( y, times, pars, strain, agegroup)
    Slu = pars.arrSlu(agegroup,:,:,:);
    Ilu = pars.arrIlu(strain,agegroup,:,:,:);
    CIlu = pars.arrCIlu(strain,agegroup,:,:,:);
    dims_Slu = prod(size(pars.arrSlu(strain,:,:,:))); %[a i j x]
    dims_Ilu = prod(size(pars.arrIlu(strain,:,:,:,:))); %[X a i j k]
    
    %infecteds = sum(y(:,dims_Slu+1:dims_Slu+dims_Ilu),2); % total infected number, regardless which strains
    infecteds = sum(y(times,Ilu),2); % total infected number, regardless which strains
    CInfecteds = sum(y(times,CIlu),2);
    
    %age structured cumulative incidence
    for a = 1:length(CIlu(1,:,1));
    CInfecteds_age(:,a) = sum(y(times,CIlu(1,a,:)),2);
    end
    CInfecteds_age = CInfecteds_age';
end