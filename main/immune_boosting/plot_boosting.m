function [ ] = plot_boosting(ab_rate)
%Plot antibody boosting probability
pa.AbB = ab_rate;
max_titre = 9;
for i=1:max_titre
    Boosting(i).abl = poisspdf([0:max_titre-i],pa.AbB)/sum(poisspdf([0:max_titre-i],pa.AbB));
end
x = 1:max_titre;
y = 1:max_titre; 
X = linspace(min(x),max(x)+1,length(x)+1); % return the index of x on 2D grid
Y = linspace(min(y),max(y)+1,length(y)+1); % return the index of y on 2D grid
%y = [0.1 0.5 3 3.2 3.5 6 7 8 9];
[X,Y] = meshgrid([x max(x)+1],[y max(y)+1]);
boosting = zeros(max(x),max(y));

for j=1:max_titre
boosting(j:max(y),j) = flipud(Boosting(j).abl);
end
%add a row
boosting(max(y)+1,:) = boosting(max(y),:);
boosting(:,max(x)+1) = boosting(:,max(x));

h = surf(X,Y,boosting);
view(0,90);
set(gca,'xlim',[1 11]);
set(gca,'XTickLabel',[0:10]);
set(gca,'ylim',[1 11]);
set(gca,'YTickLabel',[0:10]);
end