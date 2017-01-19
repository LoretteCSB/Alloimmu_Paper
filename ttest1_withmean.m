function [ pvalue ] = ttest1_withmean( meansample,stdsample,nsample,mu )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
%S=stdsample*sqrt(nsample/(nsample-1));
S=stdsample;
DF=nsample-1;
T=(meansample-mu)/S*sqrt(nsample);
pvalue= (1-tcdf(abs(T),DF))*2;
pvalue= tcdf(T,DF);

end

