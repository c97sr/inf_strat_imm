function popweights = get_popweights_4_hk_2013(pa,pop)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    popweights = zeros(1,pa.maxa);
    if pa.maxa~=5
      pop = pop(1:pa.maxa);
      popweights = pop./sum(pop);
    else
      %pop = [3009303, 8731641, 19889635, 12707659, 8453916]; 
      popweights = pop./sum(pop);
    end
end

