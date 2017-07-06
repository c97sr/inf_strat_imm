function [ ] = figureS2_1( )
% Summary of the function figure1b
% Plot seroconversion
% Written by Sean Yuan (hyuan@imperial.ac.uk) 

%global proj Antibody;
setISL;
%setup initial condition
titres_pair = [TitresTablePaired.sr_index TitresTablePaired.Age_1 TitresTablePaired.Levels_T1 TitresTablePaired.Levels_T2];
idx0 = find((titres_pair(:,3)<=1))
idx1 = find((titres_pair(:,3)>=2))

% Plot seroconversion by age
boosting_naive = titres_pair(idx0,4) - titres_pair(idx0,3);
boosting_naive_hist = histc(boosting_naive,1:9);
boosting_naive_hist = boosting_naive_hist/sum(boosting_naive_hist);

boosting_antibody = titres_pair(idx1,4) - titres_pair(idx1,3);
boosting_antibody_hist = histc(boosting_antibody,1:9);
boosting_antibody_hist = boosting_antibody_hist/sum(boosting_antibody_hist);
boosting_hist = [boosting_naive_hist boosting_antibody_hist];
H = figure;
set(H, 'Position', [500, 500, 820, 540]);

bar(boosting_hist,0.6);
%ylabel('Proportion');
%xlabel('Antibody boosting for naive individuals');

%bar(boosting_antibody_hist,0.5);
%set(gca,'XTickLabel',{'Total';'3-19';'20-39';'40-64';'65+'});
%ylabel('Proportion');
%xlabel('Antibody boosting for immune individuals');




%plot disease dynamics
%FigHandle = figure;
%set(FigHandle, 'Position', [100, 100, 1500, 360]);



end