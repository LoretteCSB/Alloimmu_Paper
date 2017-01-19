clear
close all
load('/Users/Lorette/Documents/MATLAB/Alloimmu/0_data/Data_MGH_BWH_paper_partI_20160816')
clear sta* freq*  keep precision hosp race

% Version: In this version I calculate the standard deviation


%% Now we restrict the analysis to evanescent antibodies
%
r_DB=r_DB(r_DB.hosp_acq==1&r_DB.nb_scfrom1>1&r_DB.evanes==1,:);
r_mult          = r_mult(r_mult.evanes==1,:);
r_seq           = r_seq(r_seq.evanes==1,:);
r_simult        = r_simult(r_simult.evanes==1,:);
r_seq_and_simult= r_seq_and_simult(r_seq_and_simult.evanes==1,:);
%}
%%

[n,list]=histcounts(categorical(r_DB.empi));
% nb of patients with at least 2 ab
freq_total_ab=tabulate(n(n>1));
freq_total_ab(:,3)=freq_total_ab(:,3)./100;
freq_total_ab(:,2)=[];
%% Most divergent persitence 
%  filter evanes==1 must be specifed before
[  persis_mult,stat_delta_mult ] =...
    CompareWithinVariabiliy( r_mult );

[  persis_seq,stat_delta_seq ] =...
    CompareWithinVariabiliy( r_seq);
[  persis_simult ,stat_delta_simult] =...
    CompareWithinVariabiliy( r_simult );
[  persis_both,stat_delta_both ] = ...
    CompareWithinVariabiliy( r_seq_and_simult );
stat_delta_persis = [stat_delta_seq ;stat_delta_simult;stat_delta_both;stat_delta_mult];
stat_delta_persis.Properties.RowNames={'seq','simult','both','mult'};

%% test: seq and simult are different
[ p_diff_seq_vs_simult_min ] = ranksum(persis_seq.delta_min,persis_simult.delta_min);
[ p_diff_seq_vs_simult_mid ] = ranksum(persis_seq.delta_mid,persis_simult.delta_mid);
[ p_diff_seq_vs_simult_max ] = ranksum(persis_seq.delta_max,persis_simult.delta_max);
[ p_std_seq_vs_simult_min ] = ranksum(persis_seq.std_min,persis_simult.std_min);
[ p_std_seq_vs_simult_mid ] = ranksum(persis_seq.std_mid,persis_simult.std_mid);
[ p_std_seq_vs_simult_max ] = ranksum(persis_seq.std_max,persis_simult.std_max);


%{
% ttest2_withmean( stat_delta_seq,stat_delta_simult,'max' );%2.4166e-07
[ pval_seq_vs_simult_mid ] = %ttest2_withmean( stat_delta_seq,stat_delta_simult,'mid' );%2.2399e-05
[ pval_seq_vs_simult_min ] = %ttest2_withmean( stat_delta_seq,stat_delta_simult,'min' );%0.0047

[ pval_seq_vs_both_max ] = ttest2_withmean( stat_delta_seq,stat_delta_both,'max' );%0.1903
[ pval_seq_vs_both_mid ] = ttest2_withmean( stat_delta_seq,stat_delta_both,'mid' );%0.1699
[ pval_seq_vs_both_min ] = ttest2_withmean( stat_delta_seq,stat_delta_both,'min' );%0.2183

[ pval_simult_vs_both_max ] = ttest2_withmean( stat_delta_simult,stat_delta_both,'max' );% 6.4037e-06
[ pval_simult_vs_both_mid ] = ttest2_withmean( stat_delta_simult,stat_delta_both,'mid' );%5.5151e-05
[ pval_simult_vs_both_min ] = ttest2_withmean( stat_delta_simult,stat_delta_both,'min' );%   0.0016
%}
%clear stat_delta_seq stat_delta_simult stat_delta_both


%% DELTA DAYS EXPECTED if persis is random to observe 2 days
size_simu=500;
[delta_rdn_mid,persis_random_mid] =GetDistanceBetween2Dates(r_DB.persis_dur_mid,freq_total_ab,size_simu);
[delta_rdn_min,persis_random_min] =GetDistanceBetween2Dates(r_DB.persis_dur_min,freq_total_ab,size_simu);
[delta_rdn_max,persis_random_max] =GetDistanceBetween2Dates(r_DB.persis_dur_max,freq_total_ab,size_simu);


persis_random=array2table(horzcat(persis_random_min,persis_random_mid,persis_random_max),...
    'VariableName',{'delta_min','delta_mid','delta_max'});
clear persis_random_*
%% test delta observed are not random
% seq vs random

stat=stat_delta_seq;

type='max';std_rdn=eval(strcat('delta_rdn_',type));
pval_seq_vs_random_max  = ttest1_withmean( stat.(strcat('mean_',type)),stat.(strcat('std_',type)),stat.N,std_rdn.mean);
%0.5306
type='mid';std_rdn=eval(strcat('delta_rdn_',type));
pval_seq_vs_random_mid  = ttest1_withmean( stat.(strcat('mean_',type)),stat.(strcat('std_',type)),stat.N,std_rdn.mean);
%0.0035
type='min';std_rdn=eval(strcat('delta_rdn_',type));
pval_seq_vs_random_min  = ttest1_withmean( stat.(strcat('mean_',type)),stat.(strcat('std_',type)),stat.N,std_rdn.mean);
%10-5

% simult vs random
stat=stat_delta_simult;
type='max';std_rdn=eval(strcat('delta_rdn_',type));
pval_simult_vs_random_max  = ttest1_withmean( stat.(strcat('mean_',type)),stat.(strcat('std_',type)),stat.N,std_rdn.mean);
%0
type='mid';std_rdn=eval(strcat('delta_rdn_',type));
pval_simult_vs_random_mid  = ttest1_withmean( stat.(strcat('mean_',type)),stat.(strcat('std_',type)),stat.N,std_rdn.mean);
%0
type='min';std_rdn=eval(strcat('delta_rdn_',type));
pval_simult_vs_random_min  = ttest1_withmean( stat.(strcat('mean_',type)),stat.(strcat('std_',type)),stat.N,std_rdn.mean);
%0.1979
%
stat=stat_delta_both;
type='max';std_rdn=eval(strcat('delta_rdn_',type));
pval_both_vs_random_max  = ttest1_withmean( stat.(strcat('mean_',type)),stat.(strcat('std_',type)),stat.N,std_rdn.mean);
%0.53
type='mid';std_rdn=eval(strcat('delta_rdn_',type));
pval_both_vs_random_mid  = ttest1_withmean( stat.(strcat('mean_',type)),stat.(strcat('std_',type)),stat.N,std_rdn.mean);
%0.0035
type='min';std_rdn=eval(strcat('delta_rdn_',type));
pval_both_vs_random_min  = ttest1_withmean( stat.(strcat('mean_',type)),stat.(strcat('std_',type)),stat.N,std_rdn.mean);
%10-5


stat=stat_delta_mult;

type='max';std_rdn=eval(strcat('delta_rdn_',type));
pval_mult_vs_random_max  = ttest1_withmean( stat.(strcat('mean_',type)),stat.(strcat('std_',type)),stat.N,std_rdn.mean);
%0.5306
type='mid';std_rdn=eval(strcat('delta_rdn_',type));
pval_mult_vs_random_mid  = ttest1_withmean( stat.(strcat('mean_',type)),stat.(strcat('std_',type)),stat.N,std_rdn.mean);
%0.0035
type='min';std_rdn=eval(strcat('delta_rdn_',type));
pval_mult_vs_random_min  = ttest1_withmean( stat.(strcat('mean_',type)),stat.(strcat('std_',type)),stat.N,std_rdn.mean);
%10-5

%% compare within std in persistence
% Only keep patients with at least 2 evanescent ab
r_seq=r_seq(r_seq.evanes==1,:);
[n,list]=histcounts(categorical(r_seq.empi));
r_seq(ismember(categorical(r_seq.empi),list(n==1)),:)=[];

r_simult=r_simult(r_simult.evanes==1,:);
[n,list]=histcounts(categorical(r_simult.empi));
r_simult(ismember(categorical(r_simult.empi),list(n==1)),:)=[];

%analysis can be performed for different definitions of persistence
persis='persis_dur_mid';

%for each patient caluclate the std in persistence duration
list_pt=unique(r_seq.empi);
std_seq=NaN(length(list_pt),1);
for pt=1:length(list_pt)
    std_seq(pt,1)=std(r_seq.(persis)(r_seq.empi==list_pt(pt)));
end

list_pt=unique(r_simult.empi);
std_simult=NaN(length(list_pt),1);
for pt=1:length(list_pt)
    std_simult(pt,1)=std(r_simult.(persis)(r_simult.empi==list_pt(pt)));
end

[h,p] = ttest2(std_seq,std_simult);% not signif
display('mean std for seq')
mean(std_seq)
display('mean std for simult')
mean(std_simult)
p = ranksum(std_seq,std_simult)% signif



%% DELTA DAYS EXPECTED if persis is random to observe 2 days
% I did not use that in the paper
for s=1:100
    simu_size=height(persis_seq);
    [delta_rdn_mid,persi_simu_mid] =GetDistanceBetween2Dates(r_DB.persis_dur_mid,freq_total_ab,size_simu);
    [delta_rdn_min,persi_simu_min] =GetDistanceBetween2Dates(r_DB.persis_dur_min,freq_total_ab,size_simu);
    [delta_rdn_max,persi_simu_max] =GetDistanceBetween2Dates(r_DB.persis_dur_max,freq_total_ab,size_simu);
    p_seq_vs_rdn_min(s,1)=ranksum(persis_seq.delta_min,persi_simu_min);
    p_seq_vs_rdn_mid(s,1)=ranksum(persis_seq.delta_mid,persi_simu_mid);
    p_seq_vs_rdn_max(s,1)=ranksum(persis_seq.delta_max,persi_simu_max);
end

for s=1:100
    simu_size=height(persis_simult);
    [delta_rdn_mid,persi_simu_mid] =GetDistanceBetween2Dates(r_DB.persis_dur_mid,freq_total_ab,size_simu);
    [delta_rdn_min,persi_simu_min] =GetDistanceBetween2Dates(r_DB.persis_dur_min,freq_total_ab,size_simu);
    [delta_rdn_max,persi_simu_max] =GetDistanceBetween2Dates(r_DB.persis_dur_max,freq_total_ab,size_simu);
    p_simult_vs_rdn_min(s,1)=ranksum(persis_simult.delta_min,persi_simu_min);
    p_simult_vs_rdn_mid(s,1)=ranksum(persis_simult.delta_mid,persi_simu_mid);
    p_simult_vs_rdn_max(s,1)=ranksum(persis_simult.delta_max,persi_simu_max);
end

clear stat_delta_seq stat_delta_simult stat_delta_both stat_delta_mult

%% Il me faut une table pour faire box plot avec la diff max 




%% Regression
% all patient with multiple evanescent ab
% Persis = effect pt + effectsimult +effectAb
%{
r_mult= r_mult(r_mult.evanes==1,:);
r_mult.id = nominal(r_mul.empi);
Persistence = Effet Pt + EffectSimult + Effect AbSpecificity
% comapre std among evanescent antibodies
%}
save('/Users/Lorette/Documents/MATLAB/Alloimmu/0_data/Results_MGH_BWH_paper_partIII_20160816')


