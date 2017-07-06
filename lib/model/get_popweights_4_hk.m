function popweights = get_popweights_4_hk(pa,pop)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    popweights = zeros(1,pa.maxa);
    if pa.maxa~=5
      pop = pop(1:pa.maxa);
      popweights = pop./sum(pop);
    else
      popweights = pop./sum(pop);
    end
end

