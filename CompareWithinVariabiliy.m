function [  persis,stat ] = CompareWithinVariabiliy( r_DB )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%r=r_DB(r_DB.evanes==1,{'empi','persis_dur_min','persis_dur_mid','persis_dur_max'});

r=r_DB(:,{'empi','persis_dur_min','persis_dur_mid','persis_dur_max'});

% remove patients with only 1 ab
[n,list]=histcounts(categorical(r.empi));
r(ismember(categorical(r.empi),list(n==1)),:)=[];
%calculate distance 
list= unique(r.empi);
persis=array2table(NaN(length(list),8),'VariableName',...
    {'empi','n_ab','delta_min','delta_mid','delta_max','std_min','std_mid','std_max'});
for pt=1:length(list)
    persis.empi(pt)=list(pt);
    persis.n_ab(pt)=height(r(r.empi==list(pt),:));
    persis.delta_min(pt)=max(r.persis_dur_min(r.empi==list(pt)))-min(r.persis_dur_min(r.empi==list(pt)));
    persis.delta_mid(pt)=max(r.persis_dur_mid(r.empi==list(pt)))-min(r.persis_dur_mid(r.empi==list(pt)));
    persis.delta_max(pt)=max(r.persis_dur_max(r.empi==list(pt)))-min(r.persis_dur_max(r.empi==list(pt)));
    persis.std_min(pt)=std(r.persis_dur_min(r.empi==list(pt)));
    persis.std_mid(pt)=std(r.persis_dur_mid(r.empi==list(pt)));
    persis.std_max(pt)=std(r.persis_dur_max(r.empi==list(pt)));
end

stat=[length(list),sum(persis.n_ab),mean(persis.delta_min),std(persis.delta_min),median(persis.delta_min),...
    mean(persis.delta_mid),std(persis.delta_mid),median(persis.delta_mid),...
    mean(persis.delta_max),std(persis.delta_max),median(persis.delta_max)    ];

stat=array2table(stat,'VariableName',{'N','n_ab','mean_min','std_min','median_min'...
    'mean_mid','std_mid','median_mid',...
    'mean_max','std_max','median_max'});



end

