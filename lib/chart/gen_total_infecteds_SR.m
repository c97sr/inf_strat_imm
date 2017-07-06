function [ infecteds CInfecteds CINaiveInfecteds CIImmInfecteds CINaiveInfectedsFull CIImmInfectedsFull CITreInfectedsFull CInfecteds_age] = gen_total_infecteds_SR( y, times, pars, strain, agegroup)
    Slu = pars.arrSlu(agegroup,:,:,:);
    Ilu = pars.arrIlu(strain,agegroup,:,:,:);
    CIlu = pars.arrCIlu(strain,agegroup,:,:,:);
    CINaivelu = CIlu + 40;
    CITrelu =  CIlu + 40 + 40;
    dims_Slu = prod(size(pars.arrSlu(strain,:,:,:))); %[a i j x]
    dims_Ilu = prod(size(pars.arrIlu(strain,:,:,:,:))); %[X a i j k]
    
    %infecteds = sum(y(:,dims_Slu+1:dims_Slu+dims_Ilu),2); % total infected number, regardless which strains
    infecteds = sum(y(times,Ilu),2); % total infected number, regardless which strains
    CInfecteds = sum(y(times,CIlu),2);
    for a = 1:length(CIlu(1,:,1));
    CInfecteds_age(:,a) = sum(y(times,CIlu(1,a,:)),2);
    %CInfecteds_age(:,a) = sum(y(times,CIlu(1,a,:)),2);
    %CInfecteds_age(:,a) = sum(y(times,CIlu(1,a,:)),2);
    %CInfecteds_age(:,a) = sum(y(times,CIlu(1,a,:)),2);
    end
    CInfecteds_age = CInfecteds_age';
    AbCut = 3; %0(1:5) 1(1:10) 2(1:20) 3(1:40) 
    CINaiveInfecteds = sum(y(times,CINaivelu(:,:,AbCut+1:end)),2); % Titre>=1:40 after infection boosting
    CIImmInfecteds = sum(y(times,CIlu(:,:,1:AbCut)),2);        % Titre<1:40 before boosting
    
    a = 1;
    CINaiveInfectedsFull = y(times,CINaivelu(:,a,1:end));
    CIImmInfectedsFull = y(times,CIlu(:,a,1:end));
    CITreInfectedsFull = y(times,CITrelu(:,a,1:end));
    if length(agegroup)>1
        for a=2:length(CINaivelu(1,:,1))
            CINaiveInfectedsFull = CINaiveInfectedsFull + y(times,CINaivelu(:,a,1:end));
            CIImmInfectedsFull = CIImmInfectedsFull + y(times,CIlu(:,a,1:end));
            CITreInfectedsFull = CITreInfectedsFull + y(times,CITrelu(:,a,1:end));
        end
    end
end