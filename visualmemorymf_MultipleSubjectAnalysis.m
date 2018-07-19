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
possibleFileNames = cell(length(visualmemory_subjectsRan),1);
for i = 1:length(visualmemory_subjectsRan)
    filename = strcat('data_visualmemorymf_exp_',visualmemory_subjectsRan{i},'.mat');
    possibleFileNames{i,1} = filename;
end

plotVar = 1; %if equal to 0, doesnt plot
printVar = 1; %if equal to 0, doesnt print

%preallocate a cell that will load theData structures from each participant
%into one cell
master_subjectData = cell(length(possibleFileNames),2); %files by 2 columns

% if any of files.name = possibleFilesNames then load the file and put
% into a cell array
for currfilenum = 1:numFiles
    dataFile = files{currfilenum,1};
    for i = 1:length(possibleFileNames)
        if strcmp(dataFile,possibleFileNames{i,1}) == 1
            load(dataFile)
            %index into the subjectsRan and compare to p.subject to find
            %out who is who if needed
            master_subjectData{i,1} = theData;
            master_subjectData{i,2} = subject;
            fprintf('\n%s and %s are compatible. Loading',dataFile, possibleFileNames{i,1})
        else
            fprintf('\n%s and %s not compatible',dataFile, possibleFileNames{i,1})
        end
    end
end
[subjectsLong,~] = size(master_subjectData);
p = theData(1).p;
centerContrast = (10.^linspace(log10(p.minContrast),log10(p.maxContrast),p.numContrasts));

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
        figure(1)
        set(gcf, 'Name', sprintf('Perceived (Estimated) Contrast Versus Center Contrast'));
        loglog(centerContrast,master.avgContEstimation(subjectsLong+1,:,3),'-o') %Baseline
        hold on
        if isnan(mean(subject.avgEstContrast(:,:,1),1)) == 0
        loglog(centerContrast,master.avgContEstimation(subjectsLong+1,:,1),'-o') %Perception
        end
        if isnan(mean(subject.avgEstContrast(:,:,2),1)) == 0
        hold on
        loglog(centerContrast,master.avgContEstimation(subjectsLong+1,:,2),'-o') %Working Memory
        end
        hold on 
        loglog([0.1 0.8],[0.1 0.8],'k--') %log scale line
        xlabel('Center Contrast')
        ylabel('Perceived Contrast')
        xticks([0.1 0.8]); yticks([0.1 0.8]);
        xticklabels({'10','80'});yticklabels({'10','80'});
        % Legend incorporating nans
        if sum(isnan(master.avgContEstimation(subjectsLong+1,:,1))) ~= 0 && sum(isnan(master.avgContEstimation(subjectsLong+1,:,2))) ~= 0
            legend('Baseline','Log Scale') %perception and working memory are nans
        elseif sum(isnan(master.avgContEstimation(subjectsLong+1,:,1))) ~= 0 && sum(isnan(master.avgContEstimation(subjectsLong+1,:,2))) == 0
            legend('Baseline','Working Memory','Log Scale') % perception is nan
        elseif sum(isnan(master.avgContEstimation(subjectsLong+1,:,1))) == 0 && sum(isnan(master.avgContEstimation(subjectsLong+1,:,2))) ~= 0
            legend('Baseline','Perception','Log Scale') % working memory is nan
        else
            legend('Baseline','Perception','Working Memory','Log Scale')
        end 
        hold off 
        
    % Histogram Plotting
        figure(2)
        set(gcf,'Name',sprintf('Histograms of Estimated Contrast at Each Contrast/Condition'));
        for i = 1:theData(1).p.numContrasts
            subplot(3,theData(1).p.numContrasts,i)
            hist(master.avgContEstimation(1:(end-1),i,1)) %perception
            xlim([0 1])
            subplot(3,theData(1).p.numContrasts,i+5)
            hist(master.avgContEstimation(1:(end-1),i,2)) %Working mem
            xlim([0 1 ])
            subplot(3,theData(1).p.numContrasts,i+10)
            hist(master.avgContEstimation(1:(end-1),i,3)) %baseline
            xlim([0 1])
        end               
end
 
