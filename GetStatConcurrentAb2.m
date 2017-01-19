function [stats,Freq_conc]=GetStatConcurrentAb2(r_DB,precision,month)
%get stat for concurrent antibodies
%antibody are matching if:
% - all are persistent and have the same persistence duration
% - all are evanescent and have the same persistence duration
% if some antibodies are persistent / evanescent ==> then we do not know

% r_dur: sur chaque ligne: empi, date1 et une colonne par ab avec la persis 


r_DB.date_1=round(r_DB.date_1,0);
%r_DB(r_DB.sh0lg1==999999,:)=[];%remove ab with unknown category

% Reorganize data
r=r_DB(:,{'empi','ab_name','date_1','persis_dur_mid'});
r_dur=unstack(r,'persis_dur_mid','ab_name');% empi / date1 and then the various ab in colum
ix_col_ab=width(r_dur);
r_dur.nb_simult_ab=sum(~isnan(table2array(r_dur(:,3:ix_col_ab))'))';

r=r_DB(:,{'empi','ab_name','date_1','evanes'});
r_evanes=unstack(r,'evanes','ab_name');
ix_col_ab=width(r_evanes);
r_evanes.nb_simult_ab=sum(~isnan(table2array(r_evanes(:,3:ix_col_ab))'))';
r_evanes.nb_evanes=nansum(table2array(r_evanes(:,3:ix_col_ab))')';

r_evanes.same_fate01=(r_evanes.nb_simult_ab==r_evanes.nb_evanes)+(r_evanes.nb_evanes==0);%1 if same fate / 0 otherwise
r_evanes.not_all_persis01=(r_evanes.nb_evanes~=0);
r_evanes.all_persis01=(r_evanes.nb_evanes==0);
r_evanes.all_evanes01=(r_evanes.nb_simult_ab==r_evanes.nb_evanes);
r_evanes=r_evanes(:,{'empi','date_1','nb_evanes','same_fate01','not_all_persis01','all_persis01','all_evanes01'});

r_2orMore_ab=join(r_dur,r_evanes);
r_2orMore_ab=r_2orMore_ab(r_2orMore_ab.nb_simult_ab>1,:);


r_2orMore_ab.min_persis =nanmin(table2array(r_2orMore_ab(:,3:ix_col_ab))')';
r_2orMore_ab.max_persis =nanmax(table2array(r_2orMore_ab(:,3:ix_col_ab))')';
r_2orMore_ab.diff=(max(table2array(r_2orMore_ab(:,3:ix_col_ab))')-min(table2array(r_2orMore_ab(:,3:ix_col_ab))'))';
r_2orMore_ab.strict_conc=(r_2orMore_ab.diff==0).*r_2orMore_ab.same_fate01;

r_2orMore_ab.conc_within_precision=(r_2orMore_ab.diff<=precision).*r_2orMore_ab.same_fate01;
r_2orMore_ab.discordant=(r_2orMore_ab.diff>=precision).*(r_2orMore_ab.same_fate01==0);

%% LOOK AT MOST FREQUENT PAIR AB DEVELOPED SIMULTANEOUSLY
%by antibody
clear freq_ab_concordant
ab_name=r_2orMore_ab.Properties.VariableNames(3:ix_col_ab);
for ab=3:ix_col_ab
    %freq_ab_concordant(ab-2,1)=r_2orMore_ab.Properties.VariableNames(ab);
    freq_ab_concordant(ab-2,1)=sum(~isnan(r_2orMore_ab.(ab)));
end 
freq_ab_concordant=table(freq_ab_concordant,'VariableName',{'count'});

freq_ab_concordant.name=ab_name';
freq_ab_concordant.pct=freq_ab_concordant.count./length(unique(r_2orMore_ab.empi));
freq_ab_concordant=sortrows(freq_ab_concordant,'count','descend');


%by pair
ab_name=r_2orMore_ab.Properties.VariableNames(3:ix_col_ab);
table_combi=table;row=0;
for pt=1:height(r_2orMore_ab)
    pt_ab=ab_name(find(~isnan(table2array(r_2orMore_ab(pt,3:ix_col_ab)))));
    n_ab=length(pt_ab);
    switch n_ab
    case 2 
    row=row+1; table_combi.ab(row)=strcat(pt_ab(1),'-',pt_ab(2));
    case 3
    row=row+1; table_combi.ab(row)=strcat(pt_ab(1),'-',pt_ab(2));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(1),'-',pt_ab(3));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(2),'-',pt_ab(3));
    case 4
    row=row+1; table_combi.ab(row)=strcat(pt_ab(1),'-',pt_ab(2));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(1),'-',pt_ab(3));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(1),'-',pt_ab(4));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(2),'-',pt_ab(3));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(2),'-',pt_ab(4));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(3),'-',pt_ab(4));
    case 5
    row=row+1; table_combi.ab(row)=strcat(pt_ab(1),'-',pt_ab(2));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(1),'-',pt_ab(3));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(1),'-',pt_ab(4));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(1),'-',pt_ab(5));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(2),'-',pt_ab(3));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(2),'-',pt_ab(4));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(2),'-',pt_ab(5));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(3),'-',pt_ab(4));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(3),'-',pt_ab(5));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(4),'-',pt_ab(5));
    case 6
    row=row+1; table_combi.ab(row)=strcat(pt_ab(1),'-',pt_ab(2));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(1),'-',pt_ab(3));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(1),'-',pt_ab(4));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(1),'-',pt_ab(5));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(1),'-',pt_ab(6));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(2),'-',pt_ab(3));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(2),'-',pt_ab(4));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(2),'-',pt_ab(5));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(2),'-',pt_ab(6));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(3),'-',pt_ab(4));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(3),'-',pt_ab(5));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(3),'-',pt_ab(6));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(4),'-',pt_ab(5));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(4),'-',pt_ab(6));
    row=row+1; table_combi.ab(row)=strcat(pt_ab(5),'-',pt_ab(6));
end
end
freq_pair_concordant=tabulate(table_combi.ab);
 
freq_pair_concordant=cell2table(freq_pair_concordant);
freq_pair_concordant=sortrows(freq_pair_concordant,2,'descend');

%%
ix=(r_2orMore_ab.not_all_persis01==1);%not all persistent

% FIGURES
figure
scatter(r_2orMore_ab.min_persis(ix)/month,r_2orMore_ab.max_persis(ix)/month)
hold on
plot([0, max(r_2orMore_ab.max_persis)/month],[0,max(r_2orMore_ab.max_persis)/month],'k')
plot([0, max(r_2orMore_ab.max_persis)/month],[precision/month,max(r_2orMore_ab.max_persis)/month+precision/month],'k-.')
xlim([0,36])
ylim([0,36])
xlabel({'Min'},'FontSize',21);
ylabel({'Max'},'FontSize',21);

%figure
%scatter(r_2orMore_ab.max_persis,r_2orMore_ab.diff)
%title('max persistence vs difference between the 2 antibodies evanescence date')


% stat
stats = grpstats(r_2orMore_ab,'nb_simult_ab',{'mean','std','median','iqr'},...
    'DataVars','diff');

name=r_2orMore_ab.Properties.VariableNames;

list_nb_ab=stats.nb_simult_ab;
for i=1:length(list_nb_ab);
    
    t_n_simul_ab(i,1)=list_nb_ab(i);
    t_n_pt(i,1)=length(unique(r_2orMore_ab.empi(r_2orMore_ab.nb_simult_ab==list_nb_ab(i),:)));%height(r_2orMore_ab(r_2orMore_ab.nb_simult_ab==list_nb_ab(i),:));
    t_n_ab(i,1)=sum(r_2orMore_ab.nb_simult_ab(r_2orMore_ab.nb_simult_ab==list_nb_ab(i),:));
    t_nb_event(i,1)=height(r_2orMore_ab(r_2orMore_ab.nb_simult_ab==list_nb_ab(i),:));
    
    t_n_all_evanes01(i,1)=sum(r_2orMore_ab.all_evanes01(r_2orMore_ab.nb_simult_ab==list_nb_ab(i),:));
    t_n_all_evanes01_same_date(i,1)=sum(r_2orMore_ab.all_evanes01(r_2orMore_ab.nb_simult_ab==list_nb_ab(i)&r_2orMore_ab.strict_conc==1,:));
    t_n_all_evanes01_diff_date(i,1)=sum(r_2orMore_ab.all_evanes01(r_2orMore_ab.nb_simult_ab==list_nb_ab(i)&r_2orMore_ab.strict_conc==0,:));

    t_n_all_persis01(i,1)=sum(r_2orMore_ab.all_persis01(r_2orMore_ab.nb_simult_ab==list_nb_ab(i),:));
    t_n_diff_fate(i,1)=sum(r_2orMore_ab.same_fate01==0&r_2orMore_ab.nb_simult_ab==list_nb_ab(i));
   
    t_n_not_all_persis01(i,1)=sum(r_2orMore_ab.not_all_persis01(r_2orMore_ab.nb_simult_ab==list_nb_ab(i),:));
    t_n_strict_concord(i,1)=sum(r_2orMore_ab.strict_conc(r_2orMore_ab.nb_simult_ab==list_nb_ab(i)&ix));
    t_n_concord_within_precision(i,1)=sum(r_2orMore_ab.conc_within_precision(r_2orMore_ab.nb_simult_ab==list_nb_ab(i)&ix));
    t_n_discord(i,1)=sum(r_2orMore_ab.discordant(r_2orMore_ab.nb_simult_ab==list_nb_ab(i)));
    t_frac_strict_concord(i,1)=t_n_strict_concord(i,1)/t_n_not_all_persis01(i,1);
    t_frac_conc_within_precision(i,1)=sum(r_2orMore_ab.conc_within_precision(r_2orMore_ab.nb_simult_ab==list_nb_ab(i)&ix))/t_n_not_all_persis01(i,1);
    %t_frac_strict_concord(i,1)=sum(r_2orMore_ab.strict_conc(r_2orMore_ab.nb_simult_ab==list_nb_ab(i)))/t_n_pt(i,1);
    %t_frac_conc_within_precision(i,1)=sum(r_2orMore_ab.conc_within_precision(r_2orMore_ab.nb_simult_ab==list_nb_ab(i)))/t_n_pt(i,1);
    t_frac_strict_disconcord(i,1)=sum(r_2orMore_ab.discordant(r_2orMore_ab.nb_simult_ab==list_nb_ab(i)))/t_n_pt(i,1);
    
    
end
i=i+1;

t_n_simul_ab(i,1)=mean(r_2orMore_ab.nb_simult_ab);
t_n_pt(i,1)=length(unique(r_2orMore_ab.empi));%height(r_2orMore_ab);
t_n_ab(i,1)=sum(r_2orMore_ab.nb_simult_ab);
t_n_all_persis01(i,1)=sum(r_2orMore_ab.all_persis01);
t_n_not_all_persis01(i,1)=sum(r_2orMore_ab.not_all_persis01);
t_n_strict_concord(i,1)=sum(r_2orMore_ab.strict_conc(ix));
t_n_concord_within_precision(i,1)=sum(r_2orMore_ab.strict_conc(ix));
t_n_discord(i,1)=sum(r_2orMore_ab.discordant);
t_frac_strict_concord(i,1)=t_n_strict_concord(i,1)/t_n_not_all_persis01(i,1);
t_frac_conc_within_precision(i,1)=sum(r_2orMore_ab.conc_within_precision(ix))/t_n_not_all_persis01(i,1);
t_frac_strict_disconcord(i,1)=sum(r_2orMore_ab.discordant)/t_n_pt(i,1);
t_n_all_evanes01(i,1)=sum(t_n_all_evanes01);
t_n_all_evanes01_same_date(i,1)=sum(t_n_all_evanes01_same_date);
 t_n_all_evanes01_diff_date(i,1)=sum(r_2orMore_ab.all_evanes01(r_2orMore_ab.strict_conc==0,:));

t_nb_event(i,1)=height(r_2orMore_ab);%sum(t_nb_event);
 t_n_diff_fate(i,1)=sum( t_n_diff_fate);


var_name={ 'n_simul_ab','n_pt','n_ab','n_event',...
    'n_all_persis01','n_diff_fate',...
    't_n_all_evanes01_same_date','t_n_all_evanes01_diff_date','n_all_evanes01'};
Freq_conc=table(t_n_simul_ab,t_n_pt,t_n_ab,t_nb_event, ...
t_n_all_persis01,t_n_diff_fate,...
t_n_all_evanes01_same_date, t_n_all_evanes01_diff_date, t_n_all_evanes01,'VariableName',var_name)


                     end