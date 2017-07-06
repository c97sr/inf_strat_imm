function [ titremat infectedtitremat] = gen_strain_titres_all( y, times, agestates, pars )
%return titre matrix for S + R without infected population
  nots = times;
  nostates = pars.maxi;
  titremat = zeros(length(nots), nostates);
  infectedtitremat = zeros(length(nots), nostates);
  %X = 1; %only plot for the first strain
  for i=1:length(nots)
      t = times(i);
      for l=1:nostates
          a = agestates;
          for m = 1:pars.maxj
              for n = 1:pars.maxk
                titremat(i,l) = titremat(i,l) + y(t,pars.arrSlu(a,l,m,n)) + y(t,pars.arrRlu(a,l,m,n)) + y(t,pars.arrIlu(1,a,l,m,n)); %every one go to test %
                %titremat(i,l) = titremat(i,l) + y(t,pars.arrSlu(a,l,m,n)) + y(t,pars.arrRlu(a,l,m,n)); % assuming no infected go to test serology
                infectedtitremat(i,1) = infectedtitremat(i,1) + y(t,pars.arrIlu(1,a,l,m,n));
              end
          end
      end
  end
end