function [ matY] = observe_matrix(maxtitres, err)

matY = zeros(maxtitres,maxtitres);

%iterate true titre from 1 to 10
for s = 1:maxtitres
    for k = 1:maxtitres
      
        if k == s
              p = 1-err;
        else 
              p = err/(maxtitres-1);
        end
        matY(k,s) = p;
        
    end
end
end