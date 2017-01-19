clear
close all

clear; close all;
cutoff=90;
month=365.25/12;
%%
% the first time only: generate code
%dir_name=strcat('/Users/Lorette/Documents/MATLAB/Alloimmu/');
%[ r_DB,r_dem] = CombineMGHBWHData(dir_name);
%save('/Users/Lorette/Documents/MATLAB/Alloimmu/0_data/Data_MGH_BWH_paper20160816','r_*')

load('/Users/Lorette/Documents/MATLAB/Alloimmu/0_data/Data_MGH_BWH_paper20160816','r_*')

%load('/Users/Lorette/Documents/MATLAB/Alloimmu/0_data/Data_MGH_BWH_paper','r_*')
r_DB.dod=[]; %until I inverstigate difference I delete it 

r_DB.date_1=round(r_DB.date_1,0);
filename='Data_MGH_BWH_paper_partI_20160816';
%% filter gender
%{
r_DB=r_DB(ismember(r_DB.gender01,gender),:);

if gender==0
    filename='Data_MGH_BWH_paper_partI_20160816_men';
end
if gender==1
    filename='Data_MGH_BWH_paper_partI_20160816_women';
end
%}
%% add nb ab per pt after cleaning data
r_DB.total_ab_corrected=NaN(height(r_DB),1);
list_pt=unique(r_DB.empi);
for pt=1:length(list_pt)
    r_DB.total_ab_corrected(r_DB.empi==list_pt(pt))=length(unique(r_DB.ab(r_DB.empi==list_pt(pt))));
end
r_DB.total_ab_before_clean=r_DB.total_ab;
r_DB.total_ab=r_DB.total_ab_corrected;r_DB.total_ab_corrected=[];

clear pt list_pt

list_preexisting_ab=unique(r_DB.empi(r_DB.hosp_acq==0));

%% characteristic pop
%
% j'enleve cette partie quand j'etudie les genders
[stat_desc_All,race]= GetStatDes(r_dem,r_DB);%GetDescriptiveStat(r_dem,r_DB);

% description db : ne pas utiliser le % d'ab evanescent
r_single=r_DB(r_DB.total_ab==1,:);r_dem_single=r_dem(ismember(r_dem.empi,r_single.empi),:);
r_mult=r_DB(r_DB.total_ab>1,:);r_dem_mult=r_dem(ismember(r_dem.empi,r_mult.empi),:);
stat_single_vs_mult= Compare2Gr(r_single,r_dem_single,r_mult,r_dem_mult,'single','mult' );


% frequency ab
%{
freq_ab_spec=array2table(tabulate(r_DB.ab_name),'VariableName',{'ab_name','n','pct'});
freq_ab_spec_male=array2table(tabulate(r_DB.ab_name(r_DB.gender01==0)),'VariableName',{'ab_name','n','pct'});
freq_ab_spec_female=array2table(tabulate(r_DB.ab_name(r_DB.gender01==1)),'VariableName',{'ab_name','n','pct'});
freq_ab_spec=outerjoin(freq_ab_spec,freq_ab_spec_male,'Keys','ab_name','RightVariable',{'n','pct'});
freq_ab_spec=outerjoin(freq_ab_spec,freq_ab_spec_female,'Keys','ab_name','RightVariable',{'n','pct'});
freq_ab_spec.Properties.VariableNames={'ab_name','n','pct','n_male','pct_male','n_female','pct_female'};

freq_ab_spec=GetPvalPct(freq_ab_spec);
freq_ab_spec=sortrows(freq_ab_spec,'n','Descend');

clear freq_ab_spec_*
height(r_DB(strcmp(r_DB.ab_name,'D')&r_DB.gender01==0,:));

% most frequent pair
freq_pair = GetMostFreqPair(r_DB); 
ix=find(freq_pair.total_ab<30);freq_pair(ix,:)=[];freq_pair(:,ix)=[];

freq_pair_male = GetMostFreqPair(r_DB(r_DB.gender01==0,:));
ix=find(freq_pair_male.total_ab<20);freq_pair_male(ix,:)=[];freq_pair_male(:,ix)=[];

freq_pair_female = GetMostFreqPair(r_DB(r_DB.gender01==1,:));
ix=find(freq_pair_female.total_ab<30);freq_pair_female(ix,:)=[];freq_pair_female(:,ix)=[];

clear ix
%}
%% ---------------FILTER 1 SCREEN-------------------
%% FROM NOW I FOCUS ON PATIENTS WITH AT LEAST ONE SCREEN
%%
r_DB=r_DB(r_DB.nb_scfrom1>1,:);
r_dem=r_dem(ismember(r_dem.empi,r_DB.empi),:);
stat_desc_All_sc2=GetStatDes(r_dem,r_DB);%GetDescriptiveStat(r_dem,r_DB);

% description db
r_single=r_DB(r_DB.total_ab==1,:);r_dem_single=r_dem(ismember(r_dem.empi,r_single.empi),:);
r_mult=r_DB(r_DB.total_ab>1,:);r_dem_mult=r_dem(ismember(r_dem.empi,r_mult.empi),:);
stat_single_vs_mult_2scr= Compare2Gr(r_single,r_dem_single,r_mult,r_dem_mult,'single','mult' );

%}

%% FROM NOW I FOCUS ON HOSPITAL ACQUIRED ANTIBODIES
%% single vs multiple HA
% que les patients sans pre-existing ab
r_DB=r_DB(r_DB.hosp_acq==1&~ismember(r_DB.empi,list_preexisting_ab),:);%r_DB.empi(r_DB.hosp_acq==0)),:);
% PATCH : je remplace 999999 par 0
r_DB.nb_transf_before1(r_DB.nb_transf_before1==999999)=0;
r_dem=r_dem(ismember(r_dem.empi,r_DB.empi),:);


stat_desc_All_sc2HA=GetStatDes(r_dem,r_DB);%GetDescriptiveStat(r_dem,r_DB);

r_single=r_DB(r_DB.total_ab==1,:);
r_dem_single=r_dem(ismember(r_dem.empi,r_single.empi),:);

%focus on patient with at least 2 newly acquired ab
r_mult=r_DB(r_DB.total_ab>1,:);

[n,list]=histcounts(categorical(r_mult.empi));

r_mult(ismember(categorical(r_mult.empi),list(n==1)),:)=[];
r_dem_mult=r_dem(ismember(r_dem.empi,r_mult.empi),:);

stat_single_vs_mult_2scrHA= Compare2Gr(r_single,r_dem_single,r_mult,r_dem_mult,'single','mult' );




%% sequential/simulataneous

% patient with simultaneous: 2 ab at 1 date
% patient with sequential: 2 ab at 2 different dates
%{
vl={'empi','total_ab','date_1','FU_from1','nb_scfrom1','persis',...
    'persis_dur_mid','persis_dur_min','persis_dur_max','evanes',...
    'recov','hosp_acq','ab_name','gender01','age','ab'};

t=outerjoin(r_mult(:,vl),r_mult(:,vl),'Key','empi','RightVariable',{'total_ab','date_1','FU_from1','persis',...
    'persis_dur_mid','persis_dur_min','persis_dur_max','evanes',...
    'recov','hosp_acq','ab_name','ab'});
%}
vl={'empi','date_1','FU_from1','ab_name','letter'};
r_mult.letter=cellfun(@(x) x(1),r_mult.ab_name);
t=outerjoin(r_mult(:,vl),r_mult,'Key','empi','RightVariable',{'date_1','FU_from1',...
    'ab_name','letter'});
modifiedStr= strrep(t.Properties.VariableNames, '_r_mult_1','1');modifiedStr= strrep(modifiedStr, '_r_mult_2','2');
modifiedStr= strrep(modifiedStr, '_left','1');modifiedStr= strrep(modifiedStr, '_right','2');
modifiedStr= strrep(modifiedStr, '_r_mult','2');
t.Properties.VariableNames=modifiedStr;
clear modif*

t(t.letter1>t.letter2,:)=[];
t(strcmp(t.ab_name1,t.ab_name2),:)=[];
t.diffdatetest=round(t.date_11,0)-round(t.date_12,0);

list_simult= unique(t.empi(abs(t.diffdatetest)==0)); %same date
list_seq= unique(t.empi(abs(t.diffdatetest)>0));
list_seq_simult=unique(intersect(list_seq,list_simult));
list_seq_only=unique(setdiff(list_seq,list_seq_simult));%set(A,B) data in A and not in B
list_simult_only= unique(setdiff(list_simult,list_seq_simult));

%{
t(strcmp(t.ab_name1,t.ab_name2),:)=[];
t.diffdatetest=t.date_11-t.date_12;
t(t.diffdatetest>0,:)=[];

list_simult= unique(t.empi(t.diffdatetest>=-1)); %same date
list_seq= unique(t.empi(t.diffdatetest<-1));
list_seq_simult=unique(intersect(list_seq,list_simult));
list_seq_only=unique(setdiff(list_seq,list_seq_simult));%set(A,B) data in A and not in B
list_simult_only= unique(setdiff(list_simult,list_seq_simult);
%}
r_simult_large=r_mult(ismember(r_mult.empi,list_simult),:);r_dem_ssimult_large=r_dem(ismember(r_dem.empi,list_simult),:);
r_seq=r_mult(ismember(r_mult.empi,list_seq_only),:);r_dem_seq=r_dem(ismember(r_dem.empi,list_seq_only),:);
r_simult=r_mult(ismember(r_mult.empi,list_simult_only),:);r_dem_simult=r_dem(ismember(r_dem.empi,list_simult_only),:);
r_seq_and_simult=r_mult(ismember(r_mult.empi,list_seq_simult),:);r_dem_seq_and_simult=r_dem(ismember(r_dem.empi,list_seq_simult),:);
stat_seq_and_simult_HA2sc=GetStatDes(r_dem_seq_and_simult,r_seq_and_simult);
stat_seq_vs_simult_HA2sc= Compare2Gr(r_seq,r_dem_seq,r_simult,r_dem_simult,'seq','simult' );
clear vl n list* t



save(strcat('/Users/Lorette/Documents/MATLAB/Alloimmu/0_data/',filename))


%freq multiple ab
r=r_DB(r_DB.hosp_acq==1&r_DB.nb_scfrom1>1,:);
[n,list]=histcounts(categorical(r_DB.empi));

freq_total_ab=tabulate(n);
freq_total_ab(:,3)=freq_total_ab(:,3)./100;
freq_total_ab(:,2)=[];


clearvars -except stat_*%stat_single_vs_mult_* stat_seq_vs_simult_*