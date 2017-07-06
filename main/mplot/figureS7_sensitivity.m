function [] = figureS7_sensitivity()
%read posterior likelihood
%copy mean(PosteriorSamples.LLH) to excel
%output: dscr_out/m2.xxx/ph1n1/20150709 (0707 0708 use the wrong seed number)
%read output peak

listmode = {'m1.11','m1.12','m1.16','m1','m1.17','m1.13','m1.14'};
days = [1 3 5 10 20 30 100];
listmeanllh = [];
listmaxllh = [];
listpeak = [];
for i=1:length(listmode)
%for i=1:2
    mode = char(listmode(i));
    
    out_dir = ['out/p0e05/' mode '/ph1n1/20151024/'];
    fullFileName = [out_dir 'mcmc_output_' mode '.mat'];
    if exist(fullFileName, 'file')
        %   File exists.  Do stuff....
        load(fullFileName);
    else
        out_dir = ['dscr-out/' mode '/ph1n1/20150920/']; 
        fullFileName = [out_dir 'mcmc_output_' mode '.mat'];
        if exist(fullFileName, 'file')
            %   File exists.  Do stuff....
            load(fullFileName);
        else
                out_dir = ['out/' mode '/ph1n1/20150922/'];
                fullFileName = [out_dir 'mcmc_output_' mode '.mat'];
                load(fullFileName);
        end
    end
   
    meanllh = mean(PosteriorSamples.LLH(1200:end));
    maxllh = max(PosteriorSamples.LLH);
    listmeanllh(i) = meanllh;
    listmaxllh(i) = maxllh;
    textFileName = [out_dir 'mylogs.txt'];
    fid=fopen(textFileName);
    C = textscan(fid, '%s','delimiter', '\n');
    str = C{1,1}(2);
    a = regexp(str,['\d+\.?\d+'],'match');
    peak = a{1,1}(1);
    peaklb = a{1,1}(2);
    peakub = a{1,1}(3);
    listpeak(i) = str2num(char(peak));
    listpeaklb(i) = str2num(char(peaklb));
    listpeakub(i) = str2num(char(peakub));
end

mllh = [days' listmaxllh' listmeanllh' listpeak' listpeaklb' listpeakub'];



FigH = figure;
set(FigH, 'Position', [150, 150, 1020, 360]);
hold on;
%subplot(1,2,1);
plot(mllh(:,1),mllh(:,4),'.-')
xlim([0 100]);
%ylabel('Incidence peak', 'FontSize',14);
%ylim([243 365]);
ylim([212 365]);
set(gca,'YTick',[212:30.5:365]);
set(gca,'YTickLabel',{'Aug','Sep','Oct','Nov','Dec','Jan'});
ax = gca;
%ax.XTick = days;
%ax.XTickLabel = {'-12','-10','-8','-6','-4','-2','0','2','4','6','8','10','12'};


listmode = {'m2.11','m2.12','m2.16','m2','m2.17','m2.13','m2.14'};
days = [1 3 5 10 20 30 100];
listmeanllh = [];
listmaxllh = [];
listpeak = [];
for i=1:length(listmode)
%for i=1:2
    mode = char(listmode(i));
    
    out_dir = ['out/p0e05/' mode '/ph1n1/20151024/'];
    fullFileName = [out_dir 'mcmc_output_' mode '.mat'];
    if exist(fullFileName, 'file')
        %   File exists.  Do stuff....
        load(fullFileName);
    else
        out_dir = ['dscr-out/' mode '/ph1n1/20150920/']; 
        fullFileName = [out_dir 'mcmc_output_' mode '.mat'];
        if exist(fullFileName, 'file')
            %   File exists.  Do stuff....
            load(fullFileName);
        else
                out_dir = ['out/' mode '/ph1n1/20150922/'];
                fullFileName = [out_dir 'mcmc_output_' mode '.mat'];
                load(fullFileName);
        end
    end

    
    
   
    meanllh = mean(PosteriorSamples.LLH(1200:end));
    maxllh = max(PosteriorSamples.LLH);
    listmeanllh(i) = meanllh;
    listmaxllh(i) = maxllh;
    textFileName = [out_dir 'mylogs.txt'];
    fid=fopen(textFileName);
    C = textscan(fid, '%s','delimiter', '\n');
    str = C{1,1}(2);
    a = regexp(str,['\d+\.?\d+'],'match');
    peak = a{1,1}(1);
    peaklb = a{1,1}(2);
    peakub = a{1,1}(3);
    listpeak(i) = str2num(char(peak));
    listpeaklb(i) = str2num(char(peaklb));
    listpeakub(i) = str2num(char(peakub));
end

mllh2p = [days' listmaxllh' listmeanllh' listpeak'];

plot(mllh2p(:,1),mllh2p(:,4),'.-')
%xlim([days]);
ylim([212 365]);
set(gca,'YTick',[212:30.5:365]);
set(gca,'YTickLabel',{'Aug','Sep','Oct','Nov','Dec','Jan'});
ylabel('Mean incidence peak', 'FontSize',14);
ax = gca;
%ax.XTick = days;
%ax.XTickLabel = {'-12','-10','-8','-6','-4','-2','0','2','4','6','8','10','12'};

mTextBox = uicontrol('style','text');
parentColor = [1 1 1];
xlabel('Number of Initial Infecteds T_0 (log)', 'FontSize',14);
%set(mTextBox,'String','Starting day (#weeks after 1st May)', 'Position', [290 15   500 30],'FontSize',12,'backgroundcolor',parentColor);

end
