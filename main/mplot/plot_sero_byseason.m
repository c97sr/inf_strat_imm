function [c1 H] = plot_sero_byseason(agegroup )
% The serological profiles of naive (left y-axis) and immune
% population
% Written by Sean Yuan (hyuan@imperial.ac.uk) 


%% Setup enviromental variabls
global params proj Antibody;

no_season = length(Antibody.K);
startIndex = regexp(params.startdate,'/');
start_year = params.startdate(startIndex(end)+1:end);
prompt = ['virus ' params.strain ' is selected. \n' ... 
          'Samples were recrutied from years ' num2str(str2num(start_year)) ' to ' num2str(str2num(start_year)+no_season-1) '\n'...
          'Enter which season[2-' num2str(no_season) '] the serological profile you want to view: \n'];
in_x = input(prompt);
second_collect = in_x; %strain we want to extract
init_collect = second_collect - 1;

FigL = figure;
set(FigL, 'Position', [100, 500, 820, 550]);

%init_collect = 4;
%second_collect = 5;
%third_collect = 3;
k = 2;
lastsamplingday = 150;

%retrieve parameters
pars = InitParameters(); 
pars.OutbreakStartingDay = 270;
pars.Antibody = Antibody;

%setup initial condition
%[yini age_arr] = make_ics_naive_old( pars);
if exist('agegroup') 
    id_age = find(Antibody.K(init_collect).age >= pars.ages(agegroup,1) & Antibody.K(init_collect).age < pars.ages(agegroup,2));
    ab_init = Antibody.K(init_collect).Abl(id_age)+1;
    for i = 1:pars.maxi
      Yout_k1(i) = sum(ab_init==i);
    end
    Yout_k1 = Yout_k1./sum(Yout_k1);
    
    id_age2 = find(Antibody.K(second_collect).age >= pars.ages(agegroup,1) & Antibody.K(second_collect).age < pars.ages(agegroup,2));
    ab_target = Antibody.K(second_collect).Abl(id_age2)+1;
    for i = 1:pars.maxi
      Yout_k2(i) = sum(ab_target==i);
    end
    Yout_k2 = Yout_k2./sum(Yout_k2);
else
    [yini_k1 age_arr_k1] = make_ics_fromtitres( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, Antibody.K(init_collect).Abl+1);
    [yini_k2 age_arr_k2] = make_ics_fromtitres( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, Antibody.K(second_collect).Abl+1);
end
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
if exist('agegroup') 
    
else
    Yout_k1 = retrieve_histogram(yini_k1, pars, firstsamplingday, 1, 1:pars.maxa);
    Yout_k2 = retrieve_histogram(yini_k2, pars, firstsamplingday, 1, 1:pars.maxa);
end
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
if max(max(count_titres(2:end,:)))<0.1   %max percentage of titre less than 10%
[hax,hbar1,hbar2] = plotyy(titres+0.3-0.15,count_titres(:,2),titres-0.15,count_titres(:,1),'bar','bar');
%set(hax(1),'YTickLabel',[0:100:100]);
set(hax(2),'YLim',[0 0.1]);
set(hax(2),'YTickLabel',(hax(2).YTick)*100);
set(hax(1),'YTickLabel',(hax(1).YTick)*1000);
end
if max(max(count_titres(2:end,:))) >= 0.1   %max percentage of titre between 10-20%
    if max(max(count_titres(2:end,:))) < 0.2
    count_titres(2:end,:) = count_titres(2:end,:)/2; 
    [hax,hbar1,hbar2] = plotyy(titres+0.3-0.15,count_titres(:,2),titres-0.15,count_titres(:,1),'bar','bar');

    set(hax(2),'YLim',[0 0.10001]); %adjust initial titre      
    set(hax(2),'YTick',[0:0.02:0.1]);
    set(hax(2),'YTickLabel',(hax(2).YTick)*100*2);
        
    set(hax(1),'YLim',[0 0.10001]); %adjust target titre
    set(hax(1),'YTick',[0:0.02:0.1]);
    set(hax(1),'YTickLabel',{[(hax(1).YTick)*1000 100]});
    else   %max percentage of titre larger than 20%
            count_titres(2:end,:) = count_titres(2:end,:)/4; 
            [hax,hbar1,hbar2] = plotyy(titres+0.3-0.15,count_titres(:,2),titres-0.15,count_titres(:,1),'bar','bar');
            maxY = max([hax(1).YLim(end) hax(2).YLim(end)]);
            Y1ratio = maxY./hax(1).YLim(end);
            if Y1ratio>1
                count_titres(1,2) = count_titres(1,2)*Y1ratio;
            end
            Y2ratio = maxY./hax(2).YLim(end);
            if Y2ratio>1
                count_titres(1,1) = count_titres(1,1)*Y2ratio;
            end
            [hax,hbar1,hbar2] = plotyy(titres+0.3-0.15,count_titres(:,2),titres-0.15,count_titres(:,1),'bar','bar');
            set(hax(1),'YLim',[0 maxY]);
            set(hax(2),'YLim',[0 maxY]);
            set(hax(2),'YTickLabel',(hax(2).YTick)*100*4);
            set(hax(1),'YTickLabel',(hax(1).YTick)*1000./Y1ratio./Y2ratio);
    end
end


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

axes(hax(1)); ylabel('Percentage of undetectable titre','FontSize',12);
axes(hax(2)); ylabel('Percentage of detectable titres','FontSize',12);
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
title([num2str(str2num(start_year)+in_x-2) '-' num2str(str2num(start_year)+in_x-1)]);

end