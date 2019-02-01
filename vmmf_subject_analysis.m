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
subjectProfiles=struct('SubjectName',[] , 'Order', [], 'Condition', [], 'TheData',[]); % subjProfiles = {subjectName, reportOrder, runSchedule, theData}
subjects={visualmemory_subjectsRan{1,:}};
subjectProfiles.SubjectName = cellfun(@str2double,subjects);

subjectProfiles.Order = nan(1,10);
subjectProfiles.TheData = cell(1,10);
for currSubj=1:numel(subjects)
    if exist(['data_vmmf_' subjects{currSubj} '.mat'],'file') ~= 0
        load(['data_vmmf_' subjects{currSubj}]);
    end
    subjectProfiles.Order(currSubj) = strcmp(visualmemory_subjectsRan(2,currSubj),'a'); % report order
    subjectProfiles.Condition(:,currSubj) = theData(1).p.trialSchedule; %subject condition schedule
    subjectProfiles.TheData{currSubj} = theData; %subject data
end
cd(expDir)
%% Data Organization

perceptionIndx=1; wmIndx=2; baselineIndx=3;
centerContrasts = unique(subjectProfiles{1,4}(1).p.trialEvents(:,3));
conditions = unique(subjectProfiles{1,4}(1).p.testCondition);

% setup super subject profile
superProfile = cell(sum(numRuns),numel(centerContrasts),numel(conditions),2); % [numRuns, contrasts, conditions, orders
superRunCount = 0;

for currSubj = 1:numel(subjects)
    condTrials = nan(4,5,7,2); %hard coded for now [runs, contrasts, data fields, orders]
    baselineTrials = nan(4,5,7,2);
    condCount = [0 0]; %counter for how many times the condition has been pulled
    for currRun = 1:numel(subjectProfiles{currSubj,4})
        superRunCount = superRunCount + 1; %keep track of run number for super subject
        % keep track of the run number within the order
        if currRun > 4
            trueRunNum = currRun - 4;
        else
            trueRunNum = currRun;
        end
        % Check which order to index
        orderIndx = (ceil(currRun/4))+1;
        if strcmp(visualmemory_subjectsRan{orderIndx,currSubj},'a')
            currOrder = 1;
        elseif strcmp(visualmemory_subjectsRan{orderIndx,currSubj},'b')
            currOrder = 2;
        end
        % Keep track of condition counter, #of times condition has appeared
        % in the order
        if sum(condCount == 2) >= 1
            condCount(condCount == 2) = 0;
        end
        currCond = subjectProfiles{currSubj,4}(currRun).p.testCondition;
        condCount(currCond) = condCount(currCond)+1; %counter for how many times the condition has been pulled
        currRunData = subjectProfiles{currSubj,4}(currRun).data;
        currRunEvents = subjectProfiles{currSubj,4}(currRun).p.trialEvents;
        dataFields = fieldnames(subjectProfiles{currSubj,4}(currRun).data);
        for currContrast = 1:numel(centerContrasts)
            condTrialsIndx = subjectProfiles{currSubj,4}(currRun).p.trialEvents(:,1) == currCond & subjectProfiles{currSubj,4}(currRun).p.trialEvents(:,3) == centerContrasts(currContrast);
            baselineTrialsIndx = subjectProfiles{currSubj,4}(currRun).p.trialEvents(:,1) == baselineIndx & subjectProfiles{currSubj,4}(currRun).p.trialEvents(:,3) == centerContrasts(currContrast);
            for currField=1:numel(dataFields)
                currDataField = subjectProfiles{currSubj,4}(currRun).data.(dataFields{currField});
                condTrials(trueRunNum,currContrast,currField,currOrder) = currDataField(condTrialsIndx);
                baselineTrials(trueRunNum,currContrast,currField,currOrder) = currDataField(baselineTrialsIndx);
            end
        end
    end
    
    % Data Analysis
    analysis.cond_avgs = nan(2,5,7,2); % hard-coded for now [runs per cond, contrasts, fields, orders]
    analysis.cond_error = nan(2,5,7,2);
    analysis.base_avgs = nan(2,5,7,2); % hard-coded for now [runs per cond, contrasts, fields, orders]
    analysis.base_error = nan(2,5,7,2);
    
    currTrialSchedule = subjectProfiles{currSubj,3};
%     for currOrder = 1:2
%         for currRun = 1:numel(currTrialSchedule)
%             currIndx = [currRun,currContrast,currField,currOrder];
%             currCond = currTrialSchedule(currCond);
%             for currContrast = 1:numel(centerContrasts)
%                 for currField = 1:numel(dataFields)
%                     currData = condTrials{currRun,currContrast,currField,currOrder};
%                     %                 analysis.cond_avgs(currIndx) = mean(condTrials{});
%                     %                 analysis.cond_error(currIndx) = std()/sqrt(numel(subjectProfiles{currSubj,4}));
%                     %
%                     %                 analysis.base_avgs(currIndx) = mean();
%                     %                 analysis.base_error(currIndx) = std()/sqrt(numel(subjectProfiles{currSubj,4}));
%                 end
%             end
%         end
%     end
end

%% Data Analysis

