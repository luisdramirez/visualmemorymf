%%  SINGLE SUBJECT ANALYSIS %%
clear;
close all;
expDir = pwd;
dataDir = 'data_master';
experiment = 'exp';
subjectName = '009';
whomst = subjectName;
cd(dataDir)

%% Load all subject data (both a & b trials)
% put all information into a master type structure that is seperated by trial number
load('visualmemory_condition_order')
load('visualmemory_subjectsRan')   

%Load run data
if exist(['data_visualmemorymf_' experiment '_' subjectName '.mat'],'file') ~= 0
    load(['data_visualmemorymf_' experiment '_' subjectName '.mat']);
    load('visualmemory_subjectsRan'); load('visualmemory_condition_order');
    visualmemory_condition_order = visualmemory_condition_order_real;
    runNumbers = 1:length(theData);
    [fields, numRuns] = size(theData);
else
    error('data file does not exist')
end
subjectCondSchedule = theData(1).p.trialSchedule;

%% Split trials
% Index into the condition order structure to find the persons data, and
% what their beginnning order is. 
% Trials 1-4 of the master structure can be indexed out into either a
% "structureA" or "structureB" depending on the beginning condition order
% Following this, the remaining trials (5-8) would be put into the opposite
% order matrix

if length(theData) <= 4
    runsCompleted = subjectCondSchedule(1:length(theData));
else
    runsCompleted = [subjectCondSchedule subjectCondSchedule(1:length(theData)-4)];
end

for currRun = 1:numel(theData)
    if currRun <= 4
        order1.order = visualmemory_subjectsRan{2,str2double(subjectName)};
        order1.p(currRun) = theData(currRun).p;
        order1.t(currRun) = theData(currRun).t;
        order1.EstimatedLocation(currRun,:) = theData(currRun).data.EstimatedLocation;
        order1.DifferenceLocation(currRun,:) = theData(currRun).data.DifferenceLocation;
        order1.ResponseTime_location(currRun,:) = theData(currRun).data.ResponseTime_location;
        order1.EstimatedContrast(currRun,:) = theData(currRun).data.EstimatedContrast;
        order1.DifferenceContrast(currRun,:) = theData(currRun).data.DifferenceContrast;
        order1.ResponseTime_Contrast(currRun,:) = theData(currRun).data.ResponseTime_Contrast;
        order1.responseTime(currRun,:) = theData(currRun).data.responseTime;
    else
        currRun = currRun - length(subjectCondSchedule);
        order2.order = visualmemory_subjectsRan{3,str2double(subjectName)};
        order2.p(currRun) = theData(currRun).p;
        order2.t(currRun) = theData(currRun).t;
        order2.EstimatedLocation(currRun,:) = theData(currRun).data.EstimatedLocation;
        order2.DifferenceLocation(currRun,:) = theData(currRun).data.DifferenceLocation;
        order2.ResponseTime_location(currRun,:) = theData(currRun).data.ResponseTime_location;
        order2.EstimatedContrast(currRun,:) = theData(currRun).data.EstimatedContrast;
        order2.DifferenceContrast(currRun,:) = theData(currRun).data.DifferenceContrast;
        order2.ResponseTime_Contrast(currRun,:) = theData(currRun).data.ResponseTime_Contrast;
        order2.responseTime(currRun,:) = theData(currRun).data.responseTime;
    end
end
    
%% Loop through the matrices / or display split figures (condition A and condition B)
% This would display graphs we wished to keep, and seperate each by
% the order conditon

perceptionIndex = 1;
workingmemIndex = 2;
baselineIndex = 3;

wmMat = nan(order1.p(1).numTrials/2,sum(runsCompleted==workingmemIndex));
pcMat = nan(order1.p(1).numTrials/2,sum(runsCompleted==perceptionIndex));

