%%  MULTIPLE SUBJECT ANALYSIS %%
clear;
close all;
expDir = pwd;
dataDir = 'data_master';
experiment = 'exp';
subjectName = '001';
whomst = subjectName;
cd(dataDir)

%% Load all subject data (both a & b trials)
% put all information into a master type structure that is seperated by trial number
load('visualmemory_condition_order')
load('visualmemory_subjectsRan')   
[numFiles, ~] = size(files);
possibleFileNames = cell(size(visualmemory_subjectsRan,2),1);
%preallocate a cell that will load theData structures from each participant
%into one cell
master_subjectData = cell(length(possibleFileNames),2); %files by 2 columns
for i = 1:length(visualmemory_subjectsRan)
    filename = strcat('analyzed_visualmemorymf_exp_',visualmemory_subjectsRan{1,i},'.mat');
    if exist(filename,'file') ~= 0
        load(filename); % Loads HC, test, and Regular trials
        runNumbers = 1:length(theData);
        [fields, numRuns] = size(theData);
    end
    possibleFileNames{i,1} = filename;
    master_subjectData{i,1} = theData;
    master_subjectData{i,2} = subject;
end
