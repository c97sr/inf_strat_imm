function [ ] = figureS3C_boosting(ab_rate)
% The probability distribution of antibody boosting
pa.AbB = ab_rate;
max_titre = 10;
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
set(gca,'XTickLabel',{''});
set(gca,'ylim',[1 11]);
set(gca,'YTickLabel',{''});
titres = [
          '<1:10 ';
          '1:10  ';
          '1:20  ';
          '1:40  ';
          '1:80  ';
          '1:160 ';
          '1:320 ';
          '1:640 ';
          '1:1280';
          '1:2560';
          ];

      ax = axis;    % Current axis limits
      axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
      Yl = ax(3:4);  % Y-axis limits
      % Place the text labels
      Xt = [1.5:1:10.5];
      t = text(Xt,Yl(1)*ones(1,length(Xt)),titres(1:1:10,:));
      set(t,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45);
      Xl = ax(3:4)-0.3;
      t2 = text(Xl(1)*ones(1,length(Xt)),Xt+0.5,titres(1:1:10,:));
      set(t2,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45);
      
      
      for i = 1:length(t)
        ext(i,:) = get(t(i),'Extent');
      end
      
       mTextBox = uicontrol('style','text');
       parentColor = [1 1 1];
       set(mTextBox,'String','Titres (before infection)', 'Position', [20 2  820 20],'FontSize',12,'backgroundcolor',parentColor);
end