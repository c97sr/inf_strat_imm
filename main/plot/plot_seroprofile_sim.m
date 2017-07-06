function [c1 H] = plot_seroprofile_sim( )
% Summary of the function figure1a
% Plot serological data at T1 and T2
% Written by Sean Yuan (hyuan@imperial.ac.uk) 

FigL = figure;
set(FigL, 'Position', [100, 500, 820, 550]);

%% Setup enviromental variabls
global proj Antibody;
init_collect = 1;
second_collect = 2;
third_collect = 3;
k = 2;
lastsamplingday = 150;

%retrieve parameters
pars = InitParameters(); 

%setup initial condition
[yini age_arr] = make_ics_naive( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu);
%load('sim_sample_set_model1.mat');
load('sim_sample_set_model1_ab');

%[yini_k1 age_arr_k1] = make_ics_fromtitres( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, Antibody.K(init_collect).Abl+1);
[yini_k1 age_arr_k1] = make_ics_fromtitres( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, Antibody_sim.K(init_collect).Abl+1);
%[yini_k2 age_arr_k2] = make_ics_fromtitres( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, Antibody.K(second_collect).Abl+1);
[yini_k2 age_arr_k2] = make_ics_fromtitres( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, Antibody_sim.K(second_collect).Abl+1);

%setep simulation time
T0 = pars.OutbreakStartingDay;
meanKdays(1) = mean(Antibody.K(1).numdays - T0);
meanKdays(2) = mean(Antibody.K(2).numdays - T0);
meanKdays(3) = mean(Antibody.K(3).numdays - T0);
sample_time_K1 = round(meanKdays(1));
sample_time_K2 = round(meanKdays(2));
firstsamplingday = 0;
times = firstsamplingday:1:lastsamplingday;
 
%plot disease dynamics
%H = figure;
%set(H, 'Position', [500, 500, 820, 540]);
Yout_k1 = retrieve_histogram(yini_k1, pars, firstsamplingday, 1, 1:pars.maxa);
Yout_k2 = retrieve_histogram(yini_k2, pars, firstsamplingday, 1, 1:pars.maxa);
Xout = [Yout_k1' Yout_k2'];

%calculate percentage of major serovonversion
c = Xout(:,2)-Xout(:,1);
c1 = c/sum(c(2:end));

titres = 0:1:pars.maxi; % x axis
count_titres = Xout(1:pars.maxi,:); % y axis
count_titres(1,:) = count_titres(1,:)/10; % normalize titres level 0
for i=pars.maxi:-1:2
  count_titres(i+1,:) = count_titres(i,:); 
end
count_titres(2,:) = 0;

[hax,hbar1,hbar2] = plotyy(titres+0.3-0.15,count_titres(:,2),titres-0.15,count_titres(:,1),'bar','bar');
set(hax(1),'YTickLabel',[0:20:100]);
set(hax(2),'YLim',[0 0.1]);
set(hax(2),'YTickLabel',[0:2:10]);
set(hax(1),'XLim',[-0.5 10.5]);
set(hax(1),'XTickLabel','');
set(hax(2),'XLim',[-0.5 10.5]);
set(hax(1),'XTick',[]);
set(hax(2),'XTick',[0:10]);
set(hax(2), 'Ticklength', [0 0]);
set(hax(2),'XTickLabel',['n';' ';'1';'2';'3';'4';'5';'6';'7';'8';'9']);
set(hbar1,'FaceColor', [0.6 0.8 0.9], 'EdgeColor', ['none'], 'barWidth', 0.3);
set(hbar2,'FaceColor', [0.2 0.3 0.8], 'EdgeColor', ['none'], 'barWidth', 0.3);
%set(hbar2,'FaceColor', [0.1 0.1 0.8], 'EdgeColor', ['none'], 'barWidth', 0.3);
axes(hax(1)); ylabel('Percentage for undetectable titre','FontSize',12);
axes(hax(2)); ylabel('Percentage for detectable titres','FontSize',12);
xlabel('Titres');
legend('Baseline','Follow-up');

%% Set X axis labels
Titres = [
          '<1:10 ';
          '      ';
          '1:10  ';
          '1:20  ';
          '1:40  ';
          '1:80  ';
          '1:160 ';
          '1:320 ';
          '1:640 ';
          '1:1280';
          '1:2560'
          ];
      
%% Set Text labels 
ax = axis;    % Current axis limits
axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
Yl = ax(3:4);  % Y-axis limits
% Place the text labels
Xl = [1:11];
Xt = [0.3:10.3];
t = text(Xt,Yl(1)*ones(1,length(Xt)),Titres(1:1:11,:));
set(t,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45,'FontSize',9);
set(gca,'XTickLabel','')
%xlabel('Titres');

%parentColor = get(get(h, 'parent'), 'color');
%set(h,'foregroundcolor', [0 0 0], ...
%      'backgroundcolor', parentColor);
  
mTextBox = uicontrol('style','text');
%parentColor = get(get(mTextBox, 'parent'), 'color');
parentColor = [1 1 1];
set(mTextBox,'String','Titres', 'Position', [20 2   820 20],'FontSize',12,'backgroundcolor',parentColor);
end