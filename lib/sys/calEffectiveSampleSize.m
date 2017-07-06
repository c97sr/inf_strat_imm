function [T] = calEffectiveSampleSize(out_dir)
% Summary of the function calEffectiveSampleSize
  system('set PATH=%PATH%;C:\Program Files\R\R-3.0.0\bin')
  Rcommand = ['"C:\Program Files"\R\R-3.0.0\bin\Rscript ' pwd '\R\test.R'];
  %Rcommand = ['"C:\Program Files"\R\R-3.0.0\bin\R --save ' pwd '\R\testESS.R'];
  %dos('"C:\Program Files"\R\R-3.0.0\bin\Rscript R\test.R')
  dos(Rcommand);
  outputfile = 'temp/MyData1.csv';
  out=textread(outputfile, '%s', 'whitespace',',');
  %out=textread(outputfile, '%s', 'delimiter',',');
  essvalue = {};
  %RowName = PosteriorSamples.Properties.VariableNames(1:end-1);
  %PostBar = post_bar(1:end-1);
  %T = table(RowName',PostBar',post_conf','RowNames',RowName)
  %T.Properties.VariableNames = {'Parameters' 'Mean' 'CI'};
  for i=3:length(out)-1
    if rem(i,2)==0 
        essvalue(i/2-1) = out(i);
    end
  end
  T = table(essvalue');
end