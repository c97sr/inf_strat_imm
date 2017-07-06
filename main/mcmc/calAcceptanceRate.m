llh = PosteriorSamples.LLH;
llh_1 = llh(1:end-1);
llh_2 = llh(2:end);
diff = llh_2-llh_1;
acc = 0;
length(idx)/length(llh)
