%%% vmmf_subject_analysis.m
%%% This script is designed to run analyses on data from the visual memory surround suppression task, at the individual* and group level, as specified.
%%% Subjects estimate location and contrast in a surround suppression
%%% paradigm. Report order is counter-balanced, therefore there is a group A (location -> contrast) and group B (contrast -> location).
%%% Most subjects were called back to perform the order opposite from their
%%% original. Therefore, the first 4 runs are from their original order,
%%% runs 5-8 are the opposite order.
%%% The following script will produce the following:
%%% - center vs. percieved contrast plot for subject or group (all
%%% subjects and group A/B)
%%% - variances in task performance (contrast & location error) for subject
%%% or group (all subjects and group A/B), for all conditions and each
%%% separate condition (1-perception, 2-memory, 3-baseline)

%%% *There is a super subject option that treats all data as one observer.
%% Prepare data and working directories
clear all; close all; clc;

expDir=pwd;
dataDir='data_master';
analysisLVL='subject'; % choose level of analysis ('subject'/'group'/'super')

% Load in relevant subject data
cd(dataDir)
load('visualmemory_subjectsRan')
subjects={visualmemory_subjectsRan{1,:}};
subjectProfiles=cell(numel(subjects),4); % subjProfiles = {subjectName, reportOrder, runSchedule, theData}
for nSubj=1:numel(subjects)
    if exist(['data_vmmf_' subjects{nSubj} '.mat'],'file') ~= 0
        load(['data_vmmf_' subjects{nSubj}]);
    end
    subjectProfiles{nSubj,1} = subjects{nSubj}; %subject name
    subjectProfiles{nSubj,2} = {visualmemory_subjectsRan{2:3,nSubj}}; % report order
    subjectProfiles{nSubj,3} = theData(1).p.trialSchedule; %subject condition schedule
    subjectProfiles{nSubj,4} = theData; %subject data
end
cd(expDir)

%% Organize & analyze data 

pIndx=1; wmIndx=2; blIndx=3;
centerContrasts = unique(subjectProfiles{1,4}(1).p.trialEvents(:,3));
conditions = 1:3;

% analysis matrix [perceived contrasts avg., " " STE, location difference
% avg., " " STE, contrast difference avg., " " STE]
for nSubj = 1:numel(subjects)
    for nRun = 1:numel(subjectProfiles{nSubj,4})
        currOrder = (ceil(nRun/4)); %subjectProfiles{nSubj,2}(ceil(nRun/4));
        currCond = subjectProfiles{nSubj,4}(nRun).p.testCondition;
        currRunData = subjectProfiles{nSubj,4}(nRun).data;
        currRunEvents = subjectProfiles{nSubj,4}(nRun).p.trialEvents;
        for currCond = 1:numel(conditions)
            for currContrast = 1:numel(centerContrasts)
                relevantTrials = subjectProfiles{nSubj,4}(nRun).p.trialEvents(:,1) == currCond & subjectProfiles{nSubj,4}(nRun).p.trialEvents(:,3) == centerContrasts(currContrast);
                
            end
        end
    end
end

%%
switch analysisLVL
    case 'subject'
        
    case 'group'
        
    case 'super'
end