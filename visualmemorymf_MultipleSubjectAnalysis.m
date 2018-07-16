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

plotVar = 1; %if equal to 0, doesnt plot
printVar = 1; %if equal to 0, doesnt print

%preallocate a cell that will load theData structures from each participant
%into one cell
master_subjectData = cell(1:length(possibleFileNames),2); %files by 2 columns

% if any of files.name = possibleFilesNames then load the file and put
% into a cell array
for currfilenum = 1:numFiles
    for i = length(possibleFileNames)
        dataFile = files{currfilenum,1};
        if strcmp(dataFile,possibleFileNames{i,1}) == 1
            load(dataFile)
            %index into the subjectsRan and compare to p.subject to find
            %out who is who if needed
            master_subjectData{i,1} = theData;
            master_subjectData{i,2} = subject;
        end
    end
end
[subjectsLong,~] = size(master_subjectData);

%% ANALYSIS %%

% Preallocate matrices for subject data %
master.avgContEstimation = nan(length(visualmemory_subjectsRan),master_subjectData{1,1}.p.numContrasts,3);
master.avgDiffLocation = nan(length(visualmemory_subjectsRan),master_subjectData{1,1}.p.numContrasts,3);
master.avgDiffContrast = nan(length(visualmemory_subjectsRan),master_subjectData{1,1}.p.numContrasts,3);

% Loop through the subject data and add to matrices
for subj = 1:subjectsLong
    currentSubject = master_subjectData{subj,2};%There shouldnt be any Nans here if all people have been run 4 times total. For testing purposes take out nans
    for contNum = 1:master_subjectData{1,1}.p.numContrasts
        master.avgContEstimation(subj,contNum,1) = subject.meanEstContPerception(contNum);
        master.avgContEstimation(subj,contNum,2) = subject.meanEstContWorkingMemory(contNum);
        master.avgContEstimation(subj,contNum,3) = subject.meanEstContBaseline(contNum);
        
        master.avgDiffContrast(subj,contNum,1) = subject.meanDiffContPerception(contNum);
        master.avgDiffContrast(subj,contNum,2) = subject.meanDiffContWorkingMemory(contNum);
        master.avgDiffContrast(subj,contNum,3) = subject.meanDiffContBaseline(contNum);
        
        master.avgDiffLocation(subj,contNum,1) = subject.meanEstContPerception(contNum);
        master.avgDiffLocation(subj,contNum,2) = subject.meanEstContWorkingMemory(contNum);
        master.avgDiffLocation(subj,contNum,3) = subject.meanEstContBaseline(contNum);
    end
end

%Take mean down the column of each page of each field of "master" to find
%averages over subjects. Then plot these as done in single analysis.
%Final average over all subjects will be in the last row, subjectsLong+1

master.avgContEstimation(subjectsLong+1,:,:) = mean(master.avgContEstimation,1);
master.avgDiffContrast(subjectsLong+1,:,:) = mean(master.avgDiffContrast,1);
master.avgDiffLocation(subjectsLong+1,:,:) = mean(master.avgDiffLocation,1);

if plotVar ~= 0
    % Total Avg. Estimated Contrast
    
    %%%% THIS IS FROM SINGLE SUBJECT CHANGE TO MULT COMPATABILITY %%%%
        figure(1)
        set(gcf, 'Name', sprintf('Perceived Contrast Versus Center Contrast'));
        loglog(p.centerContrast,mean(master.avgContEstimation(:,:,3),1),'-o')
        hold on
        if isnan(mean(subject.avgEstContrast(:,:,1),1)) == 0
        loglog(p.centerContrast,mean(subject.avgEstContrast(:,:,1),1),'-o')
        end
        if isnan(mean(subject.avgEstContrast(:,:,2),1)) == 0
        hold on
        loglog(p.centerContrast,mean(subject.avgEstContrast(:,:,2),1),'-o')
        end
        hold on 
        loglog([0.1 0.8],[0.1 0.8],'k--')
        xlabel('Center Contrast')
        ylabel('Perceived Contrast')
        xticks([0.1 0.8]); yticks([0.1 0.8]);
        xticklabels({'10','80'});yticklabels({'10','80'});
        legend('Baseline','Perception','Working Memory','Log Scale')
        hold off 




   
    

