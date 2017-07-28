function paper_controlPeaksOustideScanner(subjects,cfgMain)

%{
Control to check if the EGG peaks are the same inside and outside
the scanner. First calls prepro_egg_outsideScanner to preprocess EGG outside the scaner. 
Manual corrections in 4 subjects is applied. Then perform a ttest between peaks inside and outside the scanner. 
Output: results of ttest in the command line
%}


for iSubj = 1:length(subjects)
    subj_idx = subjects(iSubj)
prepro_egg_outsideScanner(subj_idx,cfgMain)
end

peaksInsideScanner = global_getEGGpeaks;
peaksOutsideScanner = global_getEGGpeaksOutsideScanner(subjects,cfgMain);

peaksOutsideScanner(6,3)= 0.042 % Corrected manually, automatic peak detected 0,033 which is not EGG
peaksOutsideScanner(13,3)= 0.039 % Corrected manually, automatic peak detected 0,066 which is not EGG
peaksOutsideScanner(23,3)= 0.059 % Corrected manually, automatic peak detected 0,033 which is not a peak
peaksInsideScanner (29,:)= [] % subject 44 outside scanner is unexploitable so we removed it from analysis of peaks inside outside
peaksOutsideScanner (29,:)= []


[h p ci stats] = ttest(peaksInsideScanner(:,3),peaksOutsideScanner(:,3))

figure
plot(peaksInsideScanner(:,3),peaksOutsideScanner(:,3),'ok','markersize',5,'LineWidth',3)
xlabel('peaks InsideScanner','fontsize',15)
ylabel('peaks OutsideScanner','fontsize',15)
xlim([0.03 0.06])
ylim([0.03 0.06])

h=gca;
set(h,'fontsize',15)
grid on
line = refline(1,0)
set(line,'color',[0 0 0],'LineWidth',2)

title(['Inside mean ' num2str(mean(peaksInsideScanner(:,3))) ' oustside ' num2str(mean(peaksOutsideScanner(:,3))) ' p ' num2str(p) ' t ' num2str(stats.tstat)])



