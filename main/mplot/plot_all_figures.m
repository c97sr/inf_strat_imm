%main script for plotting figures
%out_dir = 'out/m4.101/ph1n1/20150721/';
%load([out_dir 'mcmc_output_m4.101.mat']);
%out_dir = 'out/m4.41/ph1n1/20150806/';
%load([out_dir 'mcmc_output_m4.41.mat']);
%out_dir = 'out/m3.101/ph1n1/20150801/';
%load([out_dir 'mcmc_output_m3.101.mat']);
out_dir = 'out/m2.118/ph1n1/20150726/';
load([out_dir 'mcmc_output_m2.118.mat']);

BurnIn = 2000;
SampleNo = 20;
[c1 H] = figure1( );
savefig(H,[out_dir 'figure1.fig']);
close(H);

[FigH FigL log peak peak_lb peak_ub] = figure2(PosteriorSamples,par,BurnIn,SampleNo); %use 200 for final
savefig(FigH,[out_dir 'figure2a.fig']);
close(FigH);
savefig(FigL,[out_dir 'figure2b.fig']);
close(FigL);

fileID = fopen([out_dir 'mylogs.txt'],'w');
str = sprintf(log)
fprintf(fileID, '%s', str);
fclose(fileID);

if par.maxi > 2
setISL;
disp('save figure3');
FigH = figure3(PosteriorSamples,par,BurnIn,SampleNo);
savefig(FigH,[out_dir 'figure3.fig']);
close(FigH);

setISL;
FigH = figure4( PosteriorSamples,par,BurnIn,SampleNo);
savefig(FigH,[out_dir 'figure4.fig']);
close(FigH);
end

setISL;
disp('save figure5');
[FigH Rt] = figure5(PosteriorSamples,par,BurnIn,SampleNo);
savefig(FigH,[out_dir 'figure5.fig']);
close(FigH);
save([out_dir 'Rt_m' num2str(par.model) '.mat'], 'Rt');
