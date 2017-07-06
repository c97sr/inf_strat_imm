function [ ] = figureS1_sampling( )
% Weekly sampling distribution of baseline and follow-up 432
recruitments
FigL = figure;
set(FigL, 'Position', [100, 500, 820, 550]);

%DateString = {'01/01/2009'};
%formatIn = 'mm/dd/yyyy';
%datenum(DateString,formatIn)

setISL;

sample_numdays = Antibody.numdays;
h = histogram(sample_numdays,130,'FaceAlpha',0.75);
h.FaceColor = [0 0 1];
set(gca,'xlim',[120 423]);
set(gca,'XTick',[120:(365-60-1)/10:120+365-60]);
months = [
          'May';
          'Jun';
          'Jul';
          'Aug';
          'Sep';
          'Oct';
          'Nov';
          'Dec';
          'Jan';
          'Feb';
          '   '
          ];
     
      
      %% Set Text labels 
      ax = axis;    % Current axis limits
      axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
      Yl = ax(3:4);  % Y-axis limits
      % Place the text labels
      Xt = [120:(365-60-1)/10:120+365-60]+15;
      t = text(Xt,Yl(1)*ones(1,length(Xt)),months(1:1:11,:));
      set(t,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45);

      % Remove the default labels
      set(gca,'XTickLabel','')
      
      ylabel('#Recruited samples (weekly)')
      
             set(gca, 'FontSize', 11.5);
       
       mTextBox = uicontrol('style','text');
       parentColor = [1 1 1];
       set(mTextBox,'String','Month', 'Position', [20 2   820 20],'FontSize',12,'backgroundcolor',parentColor);
       
       
       
hold on;
plotObservedIncidence();
end

