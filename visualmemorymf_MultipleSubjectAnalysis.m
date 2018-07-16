%% MULTIPLE SUBJECT ANALYSIS %% 
% This script is designed to run based off of data calculated in the single
% subject analysis. It is designed to take all the average data collected
% from each subject who has ran (and ideally completeld 4) visualmemorymf
% trials. It will compara this data amongst subjects to produce concise
% results.

%% SETUP %%
clear;
close all;
expDir = '/Users/juliaschwartz/Desktop/visualmemorymf';
dataDir = '/Users/juliaschwartz/Desktop/visualmemorymf/data_master';
cd(expDir)
load('visualmemory_condition_order')
load('visualmemory_subjectsRan')  

cd(dataDir)
files = struct2cell(dir(dataDir))';

[numFiles, ~] = size(files);
possibleFileNames = cell(1:length(visualmemory_subjectsRan));
for i = length(visualmemory_subjectsRan)
    filename = strcat('data_visualmemorymf_exp_',visualmemory_subjectsRan{i},'.mat');
    possibleFileNames{i} = filename;
end

%preallocate a cell that will load theData structures from each participant
%into one cell
master_theData = cell(1:length(possibleFileNames));

% if any of files.name = possibleFilesNames then load the file and put
% into a cell array
for currfilenum = 1:numFiles
    for i = length(possibleFileNames)
        dataFile = files{currfilenum,1};
        if strcmp(dataFile,possibleFileNames{i,1}) == 1
            load(dataFile)
            master_theData{1} = theData;
        end
    end
end

