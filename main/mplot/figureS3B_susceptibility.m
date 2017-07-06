figure;
h1 = plot_susc_2pars(2);
h2 = plot_susc_2pars(3);
h3 = plot_susc_2pars(5);
legend([h1 h2 h3],{'1:20','1:40','1:160'})
ylabel('Susceptibility (\rho)');

set(gca,'xlim',[0 9]);
set(gca,'XTickLabel',{''});
titres = [
          '1:5   ';
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
      Xt = [0:1:9];
      Xl = ax(3:4)+0.5;
      t = text(Xt,Yl(1)*ones(1,length(Xt)),titres(1:1:10,:));
      set(t,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45);
      %Xl = ax(3:4)+0.5;
      %t2 = text(Xl(1)*ones(1,length(Xt)),Xt+0.5,titres(1:1:10,:));
      %set(t2,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45);
      
      
      for i = 1:length(t)
        ext(i,:) = get(t(i),'Extent');
      end
      
       mTextBox = uicontrol('style','text');
       parentColor = [1 1 1];
       set(mTextBox,'String','Titres', 'Position', [24 2  820 24],'FontSize',14,'backgroundcolor',parentColor);