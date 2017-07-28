% paper_randomization_randomize
%{
This script performs the group level statistics with PLV obtained with
random timeshifts to the EGG timeseries across subjects.
The PLV obtained with all possible timeshifts is obtained with 
randomization of EGG phases using the script timeseries_AllRotations_Regression
The randomization matrix is loaded, which contains for each iteration
(1:1000), which PLV of time-shifted EGG (out of 360), will be assigned to
each of the 30 subjects
It then calls paper_randomization_statsCluster to perform group level
statistics on these surragate data

%}

paper_randomization_createRandomizationMatrix % creates randomization matrix

cfgMain = global_getcfgmain
cfgMain.randomized = 1;
subjects= global_subjectList;

filename_randomizationMatrix = [global_path2root 'scripts4paper' filesep 'files' filesep 'randomizationRotationMatrix'];

load(filename_randomizationMatrix)
cfgMain.randomizationMatrix = randomizationMatrix;

for iRandomization = 1000
randomizationNumber = i;

cfgMain.currentRandomizationIteration = randomizationNumber;
paper_randomization_statsCluster(subjects,cfgMain)

end

paper_randomization_figure