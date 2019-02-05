%%% vmmf_subject_analysis.m

%% Prepare data and working directories
clear all; close all; clc;

expDir=pwd;
dataDir='data_master';

% Load in data
cd(dataDir)
load('visualmemory_subjectsRan')
subjectProfile=struct('SubjectName',[] ,'Order', [], 'Condition', [], 'TheData',[],'OrganizedData',[],'LocationError', [],'ContrastEstimate',[]);

subjectProfile.SubjectName = cell(1,10);
subjectProfile.Order = nan(1,10);
subjectProfile.TheData = cell(1,10);
subjectProfile.OrganizedData = cell(1,10);
subjectProfile.LocationError = nan(1,10);
subjectProfile.ContrastEstimate = nan(10,5,3);

subjectProfile.SubjectName = cellfun(@str2double,{visualmemory_subjectsRan{1,:}});
for currSubj=1:numel(subjectProfile.SubjectName)
    if exist(['data_vmmf_00' num2str(subjectProfile.SubjectName(currSubj)) '.mat'],'file') ~= 0
        load(['data_vmmf_00' num2str(subjectProfile.SubjectName(currSubj))]);
    end
    subjectProfile.Order(currSubj) = strcmp(visualmemory_subjectsRan(2,currSubj),'a'); % report order
    subjectProfile.Condition(:,currSubj) = theData(1).p.trialSchedule; %subject condition schedule
    subjectProfile.TheData{currSubj} = theData; %subject data
end
cd(expDir)
%% Data Organization
centerContrast = unique(subjectProfile.TheData{1}(1).p.trialEvents(:,3));

% Grab contrast estimates
for currSubj = 1:numel(subjectProfile.SubjectName)
    allData = struct('Perception',{[] [] [] [] []},'WorkingMemory',{[] [] [] [] []},'Baseline',{[] [] [] [] []}); % 5 slots for each contrast
    dataFields = fieldnames(allData);
    locationData = nan(length(subjectProfile.TheData{currSubj}(1).p.trialEvents),numel(subjectProfile.TheData{currSubj}));
    for currRun = 1:numel(subjectProfile.TheData{currSubj})
        
        % Keep track of current order
        if currRun <= 4
            currOrder = subjectProfile.Order(currSubj);
        else
            currOrder = ~subjectProfile.Order(currSubj);
        end
        
        for currField = 1:numel(dataFields) % goes through each field in allData structure (Perception, Working Memory, Baseline)
            for currContrast = 1:numel(centerContrast)
                relevantTrials = subjectProfile.TheData{currSubj}(currRun).p.trialEvents(:,1) == currField & subjectProfile.TheData{currSubj}(currRun).p.trialEvents(:,3) == centerContrast(currContrast);
                allData(currContrast).(dataFields{currField})(currRun,1) = nanmean(subjectProfile.TheData{currSubj}(currRun).data.EstimatedContrast(relevantTrials));
                allData(currContrast).(dataFields{currField})(currRun,2) = currOrder;
            end
        end
        
        locationData(:,currRun) = subjectProfile.TheData{currSubj}(currRun).data.DifferenceLocation;
    end
    
    subjectProfile.OrganizedData{currSubj} = allData;
    subjectProfile.LocationError(currSubj) = mean(mean(locationData));
    
    for currField = 1:numel(dataFields)
        for currContrast = 1:numel(centerContrast)
            subjectProfile.ContrastEstimate(currSubj,currContrast,currField) = nanmean(allData(currContrast).(dataFields{currField})(:,1),1);
        end
    end
end
