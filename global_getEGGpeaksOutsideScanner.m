function EGGpeaks =  global_getEGGpeaksOutsideScanner(subjects,cfg)

%{
Retrieve the frequency of the EGG peak recorded when subjects where outside
the scanner in order to do the analysis performed in paper_controlPeaksOustideScanner

Output:EGGpeaks First column subj _idx, second column channel number, third column Peak
frequency


IR 28/06/2017
%}

EGGpeaks = zeros(3,length(subjects))';

for iSubj = 1 : length(subjects)
    %load
    subj_idx = subjects(iSubj);
       
    filenameEGG = global_filename(subj_idx,cfg,'EGGoutsideScannerFilename');
    % Y:\Subjects\Subject13\Timeseries\EGGtimeseries\EGGoutsideScannerInfo13.mat
    
    load(filenameEGG)
    
    EGGpeaks(iSubj,1)= subj_idx;
    
    EGGpeaks(iSubj,2)= logEGGpreprocessing.bestChannel; 
    
    EGGpeaks(iSubj,3)= logEGGpreprocessing.mostPowerfullFrequency;
end

end