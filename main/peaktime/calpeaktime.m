%calculate the peak time
%http://www.me.rochester.edu/courses/ME406/webexamp5/sir1.pdf

N = 7000000;
i0 = 10/N;
s0 = (N-i0)/N;
%beta = 0.5;
gamma = 1/3.3;
%R0 = beta./gamma;
R0list = 1.22:-0.01:1.18;
for r=1:length(R0list)
R0 = R0list(r);
    %speak = gamma./beta;
speak = 1/R0;
beta = R0*gamma;

s = speak;
ipeak = i0 + 1./R0.*log(s/s0) - (s-s0);
i01 = i0 + 1./R0.*log(s0/s0) - (s0-s0);
fun = @(x) -1./(beta.*x.*(i0 + (1./R0).*log(x/s0) - (x-s0)));
%fun = @(x) exp(-x.^2).*log(x).^2;

q(r) = integral(fun,s0,speak)
end
disp ([q]);