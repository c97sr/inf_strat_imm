function [ matY] = observe_matrix(maxtitres)

matY = zeros(maxtitres,maxtitres);
EA = 0.02;
err = 0.01;

%iterate true titre from 1 to 10
for s = 1:maxtitres
    for k = 1:maxtitres
      
      if s <=maxtitres && s >= 2
        if k == s
              p = 1-EA-err;
              if s==maxtitres
                  p = 1-EA./2-err
              end
        elseif abs(k-s)==1 
              p = EA./2; 
              
        else 
                p = err/7;
                if s==10
                  p = err/8;
                end
        end
        matY(k,s) = p;
      end
      if s==1
        if k == s
              p = 1-err;
        else 
              p = err/(maxtitres-1);
        end
        matY(k,s) = p;
      end
    end
end
end