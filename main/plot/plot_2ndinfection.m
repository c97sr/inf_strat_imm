function [ ] = plot_2infection( )
% Summary of the function figure1b
% Plot seroconversion
% Written by Sean Yuan (hyuan@imperial.ac.uk) 

%global proj Antibody;
setISL;
%setup initial condition
titres_pair = [TitresTablePaired.sr_index TitresTablePaired.Age_1 TitresTablePaired.Levels_T1 TitresTablePaired.Levels_T2];
total = length(titres_pair(:,1));
total_converted = length(find((titres_pair(:,4)-titres_pair(:,3))>=2));
seroconvert = total_converted/total;
idx = find(titres_pair(:,3)>=1);
imm_converted = find((titres_pair(idx,4)-titres_pair(idx,3))>=2);
imm_seroconvert = length(imm_converted)/total;

% Plot seroconversion by age
titres_pair_age(1).titres_pair = titres_pair(find(titres_pair(:,2)>0 & titres_pair(:,2)<20),:);
titres_pair_age(2).titres_pair = titres_pair(find(titres_pair(:,2)>=20 & titres_pair(:,2)<40),:);
titres_pair_age(3).titres_pair = titres_pair(find(titres_pair(:,2)>=40 & titres_pair(:,2)<65),:);
titres_pair_age(4).titres_pair = titres_pair(find(titres_pair(:,2)>=65 & titres_pair(:,2)<100),:);
for a=1:4
    titres_pair = titres_pair_age(a).titres_pair;
    total = length(titres_pair(:,1));
    total_converted = length(find((titres_pair(:,4)-titres_pair(:,3))>=2));
    seroconvert(end+1) = total_converted/total; 
    idx = find(titres_pair(:,3)>=1);
    imm_converted = find((titres_pair(idx,4)-titres_pair(idx,3))>=2);
    imm_seroconvert(end+1) = length(imm_converted)/total;
end

H = figure;
set(H, 'Position', [500, 500, 820, 540]);
bar(seroconvert,0.45,'r');
hold on;
Hb = bar(seroconvert - imm_seroconvert, 0.45, 'b');
%set(Hb,'FaceColor',[1,1,1]*0.5);
set(gca,'XTickLabel',{'Total';'3-19';'20-39';'40-64';'65+'});
ylabel('Proportion');
xlabel('Age groups');
%set(gca,'XTickLabel',{'a'; 'b'; 'c'; 'd'})
% initial immunity prevalence from HK data
%subplot(4,4,(a-1)*4+4);
%plot_infecteds_distribution_deter(yini, pars, 1, 1, agegroup);
%ylim([0 0.2]);
%title(['HK sera(T' num2str(init_collect) ')']);
%end %--end of age loop



%plot disease dynamics
%FigHandle = figure;
%set(FigHandle, 'Position', [100, 100, 1500, 360]);



end