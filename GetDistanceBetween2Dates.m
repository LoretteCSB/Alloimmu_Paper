function [delta,delta_date] =GetDistanceBetween2Dates(persis_dur,freq_total_ab,size_simu)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

[f,x] = ecdf(persis_dur);
delta_date=[];
for ab=2:length(freq_total_ab)
    yq=rand(round(size_simu*freq_total_ab(ab,2),0),ab);
    xq=interp1(f,x,yq);
    diff_persis=(max(xq')-min(xq'))';
    delta_date=[delta_date;diff_persis];
end
    
delta.mean=mean(delta_date);
delta.std=std(delta_date);
delta.med=median(delta_date);
delta.n=length(persis_dur);

end

