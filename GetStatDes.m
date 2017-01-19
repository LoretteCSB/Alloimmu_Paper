function [ result,race ] =GetStatDes(r_dem,r_DB)
%   Detailed explanation goes here

months=365.25/12;
r_dem.race(strcmp(r_dem.race, 'African American-AFRICAN' ))={'Black'};

race=tabulate(r_dem.race);

varname={'mean','std','median'};
rowname={'nb_pt','nb_ab','nb_ab_per_pt','age','pct_female','nb_transf_before1',...
    'FU_from1','nb_scfrom1','pct_evanes','nb_evanes_2sc',...
    'pct_recov','time_to1','time_to1_men' 'time_to1_women',...
    'nb_transf_before1_men','nb_transf_before1_women','dur_to_1st_neg_scr','dur_mid'};


result(1,1)={length(unique(r_dem.empi))};%'nb_pt'
result(2,1)={height(unique(r_DB))};%nb_ab
%ab_per_pt=height(unique(r_DB))/length(unique(r_dem.empi))
%result(3,1)={ab_per_pt};%nb_ab_per_pt
r=unique(r_DB(:,{'empi','total_ab'}));
result(3,1:3)=getstat(r.total_ab);%nb_ab_per_pt
result(4,1:3)=getstat(r_DB.age);%age
result(5,1:3)={sum(r_dem.gender01)/height(r_dem)};%pct_female
result(6,1:3)=getstat(r_DB.nb_transf_before1(r_DB.nb_transf_before1<999999));%nb_transf_before1
result(7,1:3)=getstat(r_DB.FU_from1/months);%FU_from1
result(8,1:3)=getstat(r_DB.nb_scfrom1);%nb_scfrom1
%evanesc
r_DB_2sc=r_DB(r_DB.nb_scfrom1>1,:);
result(9,1)={mean(r_DB_2sc.evanes)};% pct evanesc

%redetection
r_DB_evan2sc=r_DB(r_DB.nb_scfrom_evanes>=2,:);

result(10,1)={height(unique(r_DB_evan2sc))};%nb of ab with at least 2 screens (neg screen+antiher one)
result(11,1)={mean(r_DB_evan2sc.recov)};%pct_recov

% time to 1
result(12,1:3)=getstat(r_DB.time_to1/months);
result(13,1:3)=getstat(r_DB.time_to1(r_DB.gender01==0)/months);
result(14,1:3)=getstat(r_DB.time_to1(r_DB.gender01==1)/months);
result(15,1:3)=getstat(r_DB.nb_transf_before1(r_DB.gender01==0));
result(16,1:3)=getstat(r_DB.nb_transf_before1(r_DB.gender01==1));

result(17,1:3)=getstat(r_DB.persis_dur_max(r_DB.evanes==1)/months);
result(18,1:3)=getstat(r_DB.persis_dur_mid(r_DB.evanes==1)/months);

%result(17,1:3)=getstat(r_DB.persis_dur_max(r_DB.nb_scfrom_evanes>1)/months);


result=cell2table(result,'RowNames',rowname);
result.Properties.VariableNames=varname;



end

function stat = getstat(data)
    stat={nanmean(data),nanstd(data),nanmedian(data)};
end

