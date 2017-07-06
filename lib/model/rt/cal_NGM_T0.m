function [NGM] = cal_NGM_T0(Ab,pars)
% Calclate NGM for T0
  
    Antibody = Ab;    
    %setup initial conditions
    if pars.maxi == 2 % only 2 titres
        [yini age_arr s0_imm] = make_ics_naive2titres( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, Ab.age);
    else
        [yini age_arr s0_imm] = make_ics_naive( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, Ab.age);
    end
    %[yini age_arr s0_imm] = make_ics_fullnaive( pars, pars.arrSlu, pars.arrIlu);
    it = 1;
        maxX = pars.maxX; %number of strains
        maxa = pars.maxa; %number of age groups
        maxi = pars.maxi; %Ab level to strain A
        maxj = pars.maxj; %Ab level to strain B
        maxk = pars.maxk; %Ab level to strain C
        beta = pars.beta;
        arrg = pars.arrg; %boosting array
        arrf = pars.arrf; %infectivity array
        arrh = pars.arrh; %susceptibility array
        alpha = pars.alpha; %immunity decay rate
        Tg = pars.Tg; %infectious period
        %lookup tables
        arrSlu = pars.arrSlu; 
        arrIlu = pars.arrIlu; 
        arrCIlu = pars.arrCIlu;
        
        %First calculate the force of infection
        matFOI_byS = getMOF_byS(pars,yini); %matrix of lambda

        %Second setup the variables for immune boosting
        %[InBoost OutBoost ImmDecayS ImmDecayI ] = getImmBoost(pars,yini); %Immune boosting and decay
        OutBoostProb = zeros(maxX,maxa,maxi,maxj,maxk);
        for a=1:maxa
            for i=1:maxi
                for j=1:maxj
                    for k=1:maxk
                        % Immune boosting
                        % InBoost(a,i,j,k) = 0;
                        for X=1:maxX
                            %OutBoost(X,a,i,j,k) = 0;
                            for l=1:maxi
                                for m=1:maxj
                                    for n=1:maxk
                                        OutBoostProb(X,a,i,j,k) = OutBoostProb(X,a,i,j,k) + arrg(X,a,i,j,k,l,m,n);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        %Setup matrix T(matT) and matrix Sigma(matS)
        matT = zeros(maxa*maxi, maxa*maxi);
        matS = zeros(maxa*maxi, maxa*maxi);
        X=1;
        j=1;
        k=1;

        %b infects a
        for a=1:maxa
            for b=1:maxa
                for i=1:maxi 
                    for l=1:maxi
                        %matT(bl,ai);
                        m=j;
                        n=k;   
                        idx = (a-1)*maxi+i; 
                        idy = (b-1)*maxi+l;
                        matFOI_byS(X,a,b,l,m,n);
                        S = yini(arrSlu(a,i,j,k));
                        if i==1 % immune but undetectable individuals
							S = S - pars.s0_imm(a);
                        end
                        matT(idx,idy) = matFOI_byS(X,a,b,l,m,n)*arrh(X,a,i,j,k)*(S); %S(a,i,j,k) susceptible individuals chalenged by I(X,b,l,m,n) 
                    end
                end
            end
        end


        X = 1;
        j = 1;
        k = 1;
        m = j;
        n = k;
        for a=1:maxa
            for i=1:maxi
                %matS(bl,ai);
                b=a;
                l=i;
                %should sum all the l
                %matS((a-1)*maxi+i,(b-1)*maxi+l) = -(1./(Tg))*pars.arrg(X,a,i,j,k,l,m,n)
                matS((a-1)*maxi+i,(b-1)*maxi+l) = -(1./(Tg))*OutBoostProb(X,b,l,m,n);
            end
        end
        NGM(it).A = -matT*inv(matS);
end