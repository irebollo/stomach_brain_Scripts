%{
paper_bayes_MAIN

Compute Bayes Factor for the influence of different variables on
gastric-BOLD coupling strength
1- Framewise displacement, one sample ttest
2- EGG power, frequency, traitAnxiety, BMI. Correlation
3- Gender, two sample ttest

Inputs
Framewise displacement/ coupling strenght r values
    Y:\data4paper\rFWDxCS_subjects.mat
Demographical variables
    Y:\data4paper\demographics.mat

Outputs:
    Bayes factor, output on command line

IR 07/04/2017
%}


subjects = global_subjectList
cfgMain= global_getcfgmain


%% 1 - Framewise displacement (FWD)
% one sample ttest bayes factor

% Prior H1: corresponds to an effect differing from 0 with a p-value of 0.05
nobs = 30; % number of observations
xref = +1.699; % significant t value for one-sided ttest 29 degrees of freedom reference effect size (significant effect)
% source: http://www.ttable.org/uploads/2/1/7/9/21795380/9754276.png?852




disp('FWD:')
load([global_path2root 'data4paper' filesep 'rFWDxCS_subjects.mat'])
% data obtained form script paper_control_corrMovement_CS

xdat = Mov_cs_r_gasnet_fisherZ

% perform ttest
[bf_log10]= my_ttest_bayes(xdat, xref);

% interpret bf
[res_Bayes, bf] = interpret_Bayes(bf_log10)

% intert it to 
bf_inverted_unlog = 1/(10^bf)

%% 2 - demographics
% perform JZS bayes correlation and two sample ttest

% load demograpical variables and CS

load([global_path2root 'data4paper' filesep 'demographics.mat'])
% Cs obtained from networks_couplingStrenght_obtain

% correlation
for iCluster = 1:12
bf_BMI(iCluster)= corrbf(BMIr(iCluster),30)
end

for iCluster = 1:12
bf_anxiety(iCluster)= corrbf(anxietyr(iCluster),30)
end

for iCluster = 1:12
bf_EGGpower(iCluster)= corrbf(powerEGGr(iCluster),30)
end

for iCluster = 1:12
bf_EGGPeak(iCluster)= corrbf(peakEGGr(iCluster),30)
end

%ttests
for iCluster = 1:12
bf_Gender(iCluster)= t2smpbf(stats.tstat(iCluster),15,15)
end

% p value from ttest to compare 
[h pGender ci stats] = ttest(clusterCS(~females,:),clusterCS(females,:))

bar([mean(clusterCS(~females,:));mean(clusterCS(females,:))])
shg





table4exportEffectSizes = [stats.tstat',BMIr,anxietyr,powerEGGr,peakEGGr]

table4exportBF = [bf_Gender;bf_BMI;bf_anxiety;bf_EGGpower;bf_EGGPeak]'

table4exportp = [pGender;BMIp;anxietyp;powerEGGp;peakEGGp]'





























