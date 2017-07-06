function [ PosteriorSamples ] = getOutput( out_dir, outfile, burnIn )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
PosteriorSamples = table;
outfilename = [out_dir outfile '.mat'];
if exist(outfilename, 'file') == 2
   dat = load(outfilename);
   pos_prev = dat.PosteriorSamples;
   PosteriorSamples = pos_prev(burnIn+1:end,:);
   newfileid = 2;
   for i=2:10
       outfilename = [out_dir outfile '(' num2str(i) ').mat'];
       if exist(outfilename) == 2
           pos_prev = PosteriorSamples;
           newfileid = i;
           outfilename = [out_dir outfile '(' num2str(newfileid) ').mat'];
           dat1 = load(outfilename);
           pos_new = dat1.PosteriorSamples;
           PosteriorSamples = vertcat(pos_prev, pos_new(burnIn+1:end,:));
       else
           break;
       end
   end
end
save([out_dir outfile '_final.mat'],'PosteriorSamples');
end

