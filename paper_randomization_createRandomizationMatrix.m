function    paper_randomization_createRandomizationMatrix (cfg)

%{
Creates and saves a matrix defining the time-shifts that
will be used for each subject at each iteration for the randomization
control. 
To call only once, and then load the output matrix stored in global_path2root,'Scripts\Files\randomizationMatrix'
output is stored in 
Y:\scripts\files\randomizationRotationMatrix


%}
nRandomizations = 1000;
subjects = global_subjectList;

randomizationMatrix = zeros(nRandomizations,length(subjects)); % Columns subjects, rows = iterations


for iRandomization = 1:nRandomizations
       
    randomizationMatrix (iRandomization,:) = randi([1 360],1,length(subjects));
    
end


filename_randomizationMatrix = [global_path2root 'scripts' filesep 'files' filesep 'randomizationRotationMatrix'];
save(filename_randomizationMatrix,'randomizationMatrix')

end