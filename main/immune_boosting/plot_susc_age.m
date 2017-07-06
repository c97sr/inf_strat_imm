
%asigma = 0.845;
asigma = 0;
%a = [2.844; 2.844-1.96*asigma; 2.844+1.96*asigma]
%mua = 2.844; %From ALL
%mub = 1.299; 
%mua = 3.385; %From Hobson
%mub = 2.102;
mua = 3;
%mub = 29.84;
mub = 2.102;
a = [2.15; 3.40; 2.80; 5.76; 3];
%bsigma = 0.376;
bsigma = 0;

%b = [1.299; 1.299-1.96*bsigma; 1.299+1.96*bsigma]
b = mub*ones(1,5);
%b = [mub; mub-1.96*bsigma; mub+1.96*bsigma]

T = 0:2560;
T2 = log2(T)-log2(5);
%figure;
hold on;

T2pts = 1:9;
%Tpts = 2.^(T2pts + log2(5));
Tpts = 2.^(T2pts);

for i = 1:5
	%pi = 1 - 1./(1+exp(b(i)*(log(T)-a(i))))
        %pi = 1 - 1./(1+exp(b(i)*(log2(T)*log(2)-a(i))))
        pi = 1 - 1./(1+exp(b(i)*(T2-a(i))))
	if i == 5
        	plot(T2,1-pi,'k:','LineWidth', 1.5);
	else
		plot(T2,1-pi,'LineWidth', 2);
	    %plot(T2pts,1./(1+exp(b(i)*(log2(Tpts)-a(i)))),'o'); %Plot dot based on 1:10, 1b:20,...,etc.
    end
end
set(gca,'xlim',[0 T2pts(end)]);
%set(gca,'XTick',[0:1:9]');
set(gca,'XTickLabel',{'<1:10';'1:10';'1:20';'1:40';'1:80';'1:160';'1:320';'1:640';'1:1280';'1:2560'});


legend('<20','20-39','40-64','\geq65')
xlabel('Titres');
ylabel('Susceptibility (\rho)');
