function [ result ] = Compare2Gr(r_DB1,r_dem1,r_DB2,r_dem2,suff1,suff2 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
rowname={'nb_pt','nb_ab','nb_ab_per_pt','age','pct_female','nb_transf_before1',...
    'FU_from1','nb_scfrom1','pct_evanes','nb_evanes_2sc',...
    'pct_recov','time_to1','time_to1_men' 'time_to1_women','nb_transf_before1_men',...
    'nb_transf_before1_women','dur_to_1st_neg_scr','dur_mid'};

stat_desc1=GetStatDes(r_dem1,r_DB1);
stat_desc2=GetStatDes(r_dem2,r_DB2);
% add p-values
pval=array2table(NaN(height(stat_desc2),1));
pval.Properties.RowNames=stat_desc1.Properties.RowNames;
pval.Properties.VariableNames={'pvalue'};

t1=r_DB1;t2=r_DB2;
r1=unique(t1(:,{'empi','total_ab'}));r2=unique(t2(:,{'empi','total_ab'}));

[h,p]=ttest2(t1.nb_transf_before1,t2.nb_transf_before1);
p = ranksum(t1.nb_transf_before1,t2.nb_transf_before1);
pval('nb_transf_before1','pvalue')={p};

t1=r_DB1;t2=r_DB2;
r1=unique(t1(:,{'empi','total_ab'}));r2=unique(t2(:,{'empi','total_ab'}));

[h,p]=ttest2(t1.age,t2.age);
pval('age','pvalue')={p};
[h,p]=ttest2(t1.gender01,t2.gender01);
pval('pct_female','pvalue')={p};
[h,p]=ttest2(t1.FU_from1,t2.FU_from1);
p = ranksum(t1.FU_from1,t2.FU_from1);
pval('FU_from1','pvalue')={p};

[h,p]=ttest2(t1.nb_scfrom1,t2.nb_scfrom1);
p = ranksum(t1.nb_scfrom1,t2.nb_scfrom1);
pval('nb_scfrom1','pvalue')={p};

%[h,p]=ttest2(t1.time_to1,t2.time_to1);
p = ranksum(t1.time_to1,t2.time_to1);

pval('time_to1','pvalue')={p};


%evanes
t1=t1(t1.nb_scfrom1>1,:);t2=t2(t2.nb_scfrom1>1,:);
[h,p]=ttest2(t1.evanes,t2.evanes);
pval('pct_evanes','pvalue')={p};

%[h,p]=ttest2(t1.persis_dur_max(t1.evanes==1),t2.persis_dur_max(t2.evanes==1));
p = ranksum(t1.persis_dur_max(t1.evanes==1),t2.persis_dur_max(t2.evanes==1));
pval('dur_to_1st_neg_scr','pvalue')={p};

%[h,p]=ttest2(t1.persis_dur_mid(t1.evanes==1),t2.persis_dur_mid(t2.evanes==1));
p = ranksum(t1.persis_dur_mid(t1.evanes==1),t2.persis_dur_mid(t2.evanes==1));

pval('dur_mid','pvalue')={p};


%recov
t1=t1(t1.nb_scfrom_evanes>1,:);t2=t2(t2.nb_scfrom_evanes>1,:);
[h,p]=ttest2(t1.recov,t2.recov);
pval('pct_recov','pvalue')={p};


clear r_DB1 r_DB2
clear h p t1 r1 t2 r2
stat_desc1.Properties.VariableNames=strcat(stat_desc1.Properties.VariableNames,suff1);
stat_desc2.Properties.VariableNames=strcat(stat_desc2.Properties.VariableNames,suff2);
stat_des=horzcat(stat_desc1,stat_desc2);
result=join(stat_des,pval,'Keys','RowNames');

end





