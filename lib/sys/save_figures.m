function [  ] = save_figures(PosteriorSamples,par,BurnIn,sampleno, out_dir, outfile )
%% Save the figure and close the figure windows
setISL
disp('save posterior');
[FigH] = mcmc_posterior_hist(PosteriorSamples);
print('-dpng','-r0',[out_dir 'mcmc_posterior_histogram_m' num2str(par.model) '.png'])
savefig(FigH,[out_dir 'mcmc_posterior_histogram_m' num2str(par.model) '.fig']);
close(FigH);
if (exist('sampleno'))
    SampleNo = sampleno;
else
    SampleNo = 20;    
end
setISL;
disp('save figure2');

%figure2: plot disease dynamics and serological dynamics
[FigH FigL log peak peak_lb peak_ub] = figure2(PosteriorSamples,par,BurnIn,SampleNo); %use 400 for final
savefig(FigH,[out_dir 'fig2a.fig']);
close(FigH);
savefig(FigL,[out_dir 'fig2b.fig']);
close(FigL);

output.out_dir = out_dir;
output.peak = peak;
output.peak_lb = peak_lb;
output.peak_ub = peak_ub;


fileID = fopen([out_dir 'mylogs.txt'],'w');
str = sprintf(log)
fprintf(fileID, '%s', str);
fclose(fileID);

if par.maxi > 2
setISL;
disp('save figure1');
FigH = figure1(PosteriorSamples,par,BurnIn,SampleNo);
savefig(FigH,[out_dir 'fig1.fig']);
close(FigH);

setISL;
[S1 T1]= table1(PosteriorSamples, BurnIn);
writetable(T1,[out_dir 'Table1.csv'],'WriteRowNames',false);

setISL;
[Table2Left Table2Right p] = table2( PosteriorSamples, par,BurnIn,SampleNo);
writetable(Table2Left,[out_dir 'Table2a.csv']);
writetable(Table2Right,[out_dir 'Table2b.csv']);
end

setISL;
[FigH seroprev_diff] = figureS6( PosteriorSamples,par,BurnIn,SampleNo);
savefig(FigH,[out_dir 'figS6.fig']);
close(FigH);

save([out_dir outfile '_peak.mat'],'output','seroprev_diff');

setISL;
disp('save figure5');
[FigH Rt] = figure5(PosteriorSamples,par,BurnIn,SampleNo);
savefig(FigH,[out_dir 'fig5.fig']);
close(FigH);
save([out_dir 'Rt_m' num2str(par.model) '.mat'], 'Rt');

setISL; 
csvwrite(['temp/posterior.csv'],table2array(PosteriorSamples));
Tess = calEffectiveSampleSize; % call R functions
writetable(Tess,[out_dir 'TableESS.csv'],'WriteRowNames',false);
end