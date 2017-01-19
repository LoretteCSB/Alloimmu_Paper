clear
close all
load('/Users/Lorette/Documents/MATLAB/Alloimmu/0_data/Data_MGH_BWH_paper_partI_20160816')
%load('/Users/Lorette/Documents/MATLAB/Alloimmu/0_data/Data_MGH_BWH_paper_partI')

clear sta* freq*  keep precision hosp race
addpath('/Users/Lorette/Documents/MATLAB/fileexchange')

% At this stage, r_DB should only contain patients with HA and more than 1 sc,
%but to be sure, I still put the filter
r_DB(r_DB.hosp_acq==0&r_DB.nb_scfrom1==1,:)=[];

%list_persis=unique(r_DB.empi(r_DB.persis==1));
%r_DB(ismember(r_DB.empi,list_persis),:)=[];
%}
%% filter gender
gender=[0 1];
r_single=r_single(ismember(r_single.gender01,gender),:);
r_seq=r_seq(ismember(r_seq.gender01,gender),:);
r_simult=r_simult(ismember(r_simult.gender01,gender),:);
r_mult=r_mult(ismember(r_mult.gender01,gender),:);

filename='Results_MGH_BWH_paper_partII_20160816';
pic='survival_curve_persistence';
if gender==0
    filename='Results_MGH_BWH_paper_partII_men';
    pic='survival_curve_persistence_men';
end
if gender==1
    filename='Results_MGH_BWH_paper_partII_women';
    pic='survival_curve_persistence_women';
end
%% SURVIVAL CURVE FOR ANTIBODY PERSISTENCE
%
[Ssingle,xsingle,Smin,Smax] = ecdf(r_single.persis_dur_mid,'function','survivor','censoring',r_single.persis);%GetSurvivalCurve(HA(ismember(HA.empi,list_pt_single),:),'persis_dur_mid','persis',keep,'Time to evanescence',month,dir_name,hosp,withIC);
[Sseq,xseq,Smin,Smax] = ecdf(r_seq.persis_dur_mid,'function','survivor','censoring',r_seq.persis);%GetSurvivalCurve(HA(ismember(HA.empi,list_pt_sequentially),:),'persis_dur_mid','persis',keep,'Time to evanescence',month,dir_name,hosp,withIC);
[Ssimult,xsimult,Smin,Smax] = ecdf(r_simult.persis_dur_mid,'function','survivor','censoring',r_simult.persis);%GetSurvivalCurve(HA(ismember(HA.empi,list_pt_simult),:),'persis_dur_mid','persis',keep,'Time to evanescence',month,dir_name,hosp,withIC);
[Smult,xmult,Smin,Smax] = ecdf(r_mult.persis_dur_mid,'function','survivor','censoring',r_mult.persis);%GetSurvivalCurve(HA(ismember(HA.empi,list_pt_with_both),:),'persis_dur_mid','persis',keep,'Time to evanescence',month,dir_name,hosp,withIC);

figure
stairs(xsimult/month,Ssimult,'b','LineWidth',1.5)
hold on
stairs(xsingle/month,Ssingle,'k','LineWidth',1.5)
stairs(xseq/month,Sseq,'r','LineWidth',1.5)

%stairs(xmult/month,Smult,'--','LineWidth',1.5)

legend('Simultaneous antibodies','Single antibody','Sequential antibodies')%,'Multiple')
legend('boxoff')
xlabel('Time before evanescence (month)')
ylabel('% ab detected')
set(gca,'FontSize',14)
ylim([0 1])
hold off
saveas(gcf,strcat('/Users/Lorette/Documents/POSTDOC/1-Alloimmunization/report/figures paper/',pic),'epsc')

% log rank
xsingle=table2array(r_single(:,{'persis_dur_mid','persis'}));
xseq=table2array(r_seq(:,{'persis_dur_mid','persis'}));
xsimult=table2array(r_simult(:,{'persis_dur_mid','persis'}));
xmult=table2array(r_mult(:,{'persis_dur_mid','persis'}));

figure
display('*****************************************************')
display('*****  SINGLE VERSUS MULTIPLE **********')
logrank(xsingle,xmult)
display('*****************************************************')
display('*****  SINGLE VERSUS SEQUENTIAL**********')
logrank(xsingle,xseq)
display('*****************************************************')
display('******   SINGLE VERSUS SIMULT **********')
logrank(xsingle,xsimult)
display('*****************************************************')
display('******  SEQUENTIAL VERSUS SIMULTANEOUS **********')
logrank(xseq,xsimult)
% Syntax: 	logrank(x1,x2,alpha,censflag)
% Inputs:
%           X1 and X2 (mandatory)- Nx2 data matrix:
%           (X:,1) = survival time of the i-th subject
%           (X:,2) = censored flag 
%           (0 if not censored; 1 if censored)
%           ALPHA (optional) - significance level (default 0.05) 
%           CENSFLAG (optional) - Censored Plot flag (default 0). If 0
%           censored data will be plotted spreaded on the horizontal
%           segment; if 1 they will be plotted at the given time of
%           censoring.
%}
%% FATE ANTIBODIES 

%%
precision=7;
[stat_simultAb,Freq_simultAb]=GetStatConcurrentAb2(r_DB,precision,month);


list_pt=unique(r_DB.empi);
r_DB= sortrows(r_DB,'date_1','ascend');
r_DB.date1 = r_DB.date_1;

for pt=1:length(list_pt)
    list_date_pt=r_DB.date1(r_DB.empi==list_pt(pt));
    for date_pt=1:length(list_date_pt)
        r_DB.date1(r_DB.empi==list_pt(pt)&r_DB.date_1==list_date_pt(date_pt))=date_pt;
    end
end
r_DB=r_DB(:,{'empi','date1','ab','total_ab','gender01', 'ab_name','evanes','recov'});
list_pt_Ab=r_DB.empi(ismember(r_DB.ab_name,'Cw'));
length(unique(list_pt_Ab))

%% LOOKING AT PAIRS

%'empi','date_1','persis_dur_mid','evanes','evanes_dur_mid','evanes_dur_min','evanes_dur_max','nb_scfrom_evanes','recov','transf_to_recov','unit_to_recov','persis2_dur_mid','persis2_dur_min','persis2_dur_max','persis2','evanes2','FU_from_redect','scr_from_recov','hosp_acq','gr_persis_min','gr_persis_max','persisminmax','ab_name','cold','clinically_relevant','bening_ab','auto','D','NatOccur','other','gender','dob','race','vital_status','dod','gender01','age','minWBC_transf1','maxWBC_transf1','medianWBC_transf1','meanWBC_transf1','precise','ab','total_ab_before_clean'}


r_simult_large.letter=cellfun(@(x) x(1),r_simult_large.ab_name);
vl={'empi','date_1','persis_dur_max','evanes','ab_name','letter'};

t=outerjoin(r_simult_large(:,vl),r_simult_large(:,vl),'Key','empi','RightVariable',...
    {'date_1','persis_dur_max','evanes','ab_name','letter'});
modifiedStr= strrep(t.Properties.VariableNames, '_r_simult_large_1','1');modifiedStr= strrep(modifiedStr, '_r_simult_large','2');
modifiedStr= strrep(modifiedStr, '_left','1');modifiedStr= strrep(modifiedStr, '_right','2');
t.Properties.VariableNames=modifiedStr;
clear modif*

t(t.letter1>t.letter2,:)=[];t(:,{'letter1','letter2'})=[];
t(strcmp(t.ab_name1,t.ab_name2),:)=[];

t.diffdatetest=round(t.date_11,0)-round(t.date_12,0);
t(abs(t.diffdatetest)>0,:)=[];
t.diffpersis=t.persis_dur_max1-t.persis_dur_max2;
t.both_evanes=t.evanes1==1&t.evanes2==1;
t.both_persis=t.evanes1==0&t.evanes2==0;
t.diff_fate=(t.evanes1==0&t.evanes2==1)|(t.evanes1==1&t.evanes2==0);

% number of pair
display(strcat('# simult pairs:',num2str(height(t)))) 
display(strcat('# pairs evanesc same date:',num2str(height(t(t.both_evanes==1&t.diffpersis<1,:))))) 
display(strcat('# pairs persis:',num2str(height(t(t.both_persis,:))))) 
display(strcat('# pairs evanesc 1 to 15d:',num2str(height(t(t.both_evanes==1&t.diffpersis>=1&t.both_evanes==1&t.diffpersis<15,:))))) 
display(strcat('# pairs evanes diff date:',num2str(height(t(t.both_persis==1,:))))) 
display(strcat('# pairs different fate:',num2str(height(t(t.diff_fate==1,:))))) 


list_simult= unique(t.empi(abs(t.diffdatetest)<1));

t=t(ismember(t.empi,list_simult),:);
%}

save(strcat('/Users/Lorette/Documents/MATLAB/Alloimmu/0_data/',filename))
