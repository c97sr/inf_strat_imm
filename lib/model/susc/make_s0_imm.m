function [s0_imm ] = make_s0_imm( pa )
%make_g_simple Summary of this function goes here
% Function to return the most simple immune state boost
% Written by Sean Yuan (hyuan@imperial.ac.uk) 


%s0_imm = pa.initS(:,1)'.*pa.age_arr*pa.PUAb;
if isfield(pa,'age_arr') == 0
    s0_imm = zeros(1,pa.maxa);
else
    s0_imm = pa.age_arr*pa.PUAb;
end

end
