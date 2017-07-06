function [] = figure5_Rt()
%Summary of this function goes here
% Comparison of the effective reproductive numbers between the titre
% and threshold models.
%1 Feb 2015
% Written by Sean Yuan (hyuan@imperial.ac.uk)

%chech here for filling the area between the curves,
%http://stackoverflow.com/questions/6245626/matlab-filling-in-the-area-between-two-sets-of-data-lines-in-one-figure
%http://stackoverflow.com/questions/19910982/matlab-fill-area-between-lines

H = figure;
set(H, 'Position', [500, 500, 820, 540]);

%%the first model
out_dir1 = ['out/p0e05/m1/ph1n1/20151024/'];
load([out_dir1 'Rt_m1.mat']);
%out_dir1 = ['out/imm/m4/ph1n1/20160304/'];
%load([out_dir1 'Rt_m4.mat']);
Rt_lower = quantile(Rt,0.025);
Rt_upper = quantile(Rt,0.975);
Rt_mean = mean(Rt)
plot(Rt_mean,'b', 'LineWidth', 2);
hold on;
baseLine = 0.2;
x = 1:360;
h1 = fill([x fliplr(x)], [Rt_upper fliplr(Rt_lower)], 'b','EdgeColor','none');
alpha(.3);


%%the second model
%out_dir2 = ['out/p0e05/m2/ph1n1/20151024/'];
%load([out_dir2 'Rt_m2.mat']);
out_dir2 = ['out/classic/m6/ph1n1/20160307/'];
load([out_dir2 'Rt_m6.mat']);
Rt_lower = quantile(Rt,0.025);
Rt_upper = quantile(Rt,0.975);
Rt_mean = mean(Rt)
plot(Rt_mean,'r', 'LineWidth', 2);
baseLine = 0.2;
x = 1:360;
h2 = fill([x fliplr(x)], [Rt_upper fliplr(Rt_lower)], 'r','EdgeColor','none');
alpha(.3);

%% Still working on this
%%the third model
%out_dir3 = ['out/m3101/ph1n1/20150722/'];
%load([out_dir3 'Rt_m3.101.mat']);
%Rt_lower = quantile(Rt,0.025);
%Rt_upper = quantile(Rt,0.975);
%Rt_mean = mean(Rt)
%plot(Rt_mean,'k');
%baseLine = 0.2;
%x = 1:360;
%h1 = fill([x fliplr(x)], [Rt_upper fliplr(Rt_lower)], 'r','EdgeColor','none');
%alpha(.25);



%% add the ticks 
legend([h1 h2],{'$R_{B}(t)$','$R_{C}(t)$'},'Interpreter','latex');  % Only the blue and green lines appear
%legend({'$\bar{R}_{B}$','$R_{B}$ CI','$\bar{R}_{SIR}$','$R_{SIR}$ CI','$\bar{R}_{SIR}$','$R_{SIR}$ CI'},'Interpreter','latex')

lastsamplingday = 305;
xlim([1 lastsamplingday]);
ylim([0.6 1.4]);
set(gca,'XTick',[1:30.5:lastsamplingday]);
%set('YTickLabel',);

    Xl = [1 lastsamplingday];
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
     
      ax = axis;    % Current axis limits
      axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
      Yl = ax(3:4);  % Y-axis limits
      % Place the text labels
      Xt = [1:(365-60-1)/10:365-60]+15;
      t = text(Xt,Yl(1)*ones(1,length(Xt)),months(1:1:11,:));
      set(t,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45);
      
      % Remove the default labels
      set(gca,'XTickLabel','')
      mTextBox = uicontrol('style','text');
      parentColor = [1 1 1];
      set(mTextBox,'String','Month', 'Position', [20 8   800 24],'FontSize',14,'backgroundcolor',parentColor)
      ylabel('Effective reproductive number','FontSize',14);
end

