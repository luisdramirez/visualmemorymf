%%  MASTER ANALYSIS %%
% This script is meant to loop through all trials with either hardcoded
% test runs, test runs, or actual experiment runs. It will show represent
% meaningful data within one single run, or if there are multiple runs for
% one single subject/experiment matching, it will loop through to collect
% meaningful results between data sets.

%% SETUP %%
% Preliminary data loading and setup %
clear;
close all;
expDir = '/Users/juliaschwartz/Desktop/visualmemorymf';
dataDir = 'data_master';
allP.experiment = 'test';
allP.subject = 'JS';
cd(dataDir)

%Load run data
if exist(['data_visualmemorymf_' allP.experiment '_' allP.subject '.mat'],'file') ~= 0
    load(['data_visualmemorymf_' allP.experiment '_' allP.subject '.mat']); % Loads HC, test, and Regular trials
    cd(expDir)
    load('visualmemory_subjectsRan'); load('visualmemory_condition_order');
    runNumbers = 1:length(theData);
    [fields, nTrials] = size(theData);
else
    error('data file does not exist')
end

%Plotting and Printing Settings
plotVar = 0; %Set equal to 1 to display plots, equal to 0 to not display plots.
printVar = 0; %Set equal to 0 to not print information, equal to 1 to print information.

% Pre-allocate Data cells.
allP = cell(1,length(runNumbers)); % Parameters
allT = cell(1,length(runNumbers)); % Experiment Timing Information
allData = cell(1,length(runNumbers)); % Subject-entered Data
stats = cell(1,length(runNumbers));

% Assign cells based off of data size.
for nRun = 1:nTrials
    allP{nRun} = theData(nRun).p;
    allT{nRun} = theData(nRun).t;
    allData{nRun} = theData(nRun).data;
    stats{nRun} = theData(nRun).stats;
end

%Finding subject name and indexing to condition order
condIndex = find(strcmp(visualmemory_subjectsRan,allP{1,1}.subject));
if condIndex > 24
    condIndex = condIndex - 24; %The condition order resets after 24, this matches to the reset.
end
subjectCondSchedule = visualmemory_condition_order(condIndex,:); %Gives the current condition schedule, indexes to the row we are on, columns 1-4 represent the condition for each run.


%Save out Relevant Information
subject.avgEstContrast = nan(nTrials,theData(1).p.numContrasts,3);
subject.avgDiffContrast = nan(nTrials,theData(1).p.numContrasts,3);
subject.avgDiffLoc = nan(nTrials,theData(1).p.numContrasts,3);

%% ANALYSIS LOOP %%
% Will include 2 nested for loops.
    % Main for loop: Goes through all runs completed by a single subject.
        % Inner Condition Loop, uses the subject condition schedule, loops
        % through each column to access the current condition.
        % 1 = PERCEPTION
        % 2 = WORKING MEMORY
        % 3 = BASELINE
 
        
% Shortens the condition schedule to only go through trials already ran.
 if nTrials < 4
     subjectCondSchedule = subjectCondSchedule(1:nTrials);
 end
 
 % Main for loop.
 for nRun = 1:nTrials
     thisRunsCond = subjectCondSchedule(nRun);
    if thisRunsCond == 1
        variableCondition = 'Perception';
        baselineCondition = 'Baseline';
    elseif thisRunsCond == 2
        variableCondition = 'Working Memory';
        baselineCondition = 'Baseline';
    else 
        error('Condition numbers/scheduling are set up incorrectly.')
    end
    p = allP{nRun};
    t = allT{nRun};
    data = allData{nRun};
    data = cell2mat(struct2cell(data));
    data = data';
    [dataTrials,dataParams] = size(data);
    p.centerContrast = [10.^linspace(log10(p.minContrast),log10(p.maxContrast),p.numContrasts)];
    
    % Add on Trial Number to end of p.trialEvents (sorting purposes)
    trialNum = (1:length(p.trialEvents))';
    p.trialEvents = [p.trialEvents trialNum];
    
    % Seperate Variable and Baseline Condition
    baseline.TE = p.trialEvents((p.trialEvents(:,1)==3),:);
    baseline.Data = zeros(length(baseline.TE),dataParams);
    for j = 1:length(baseline.TE)
        baseline.Data(j,:) = data(baseline.TE(j,6),:);
    end
    if strcmp(variableCondition,'Perception') == 1
        perception.TE = p.trialEvents((p.trialEvents(:,1)==thisRunsCond),:);
        perception.Data = zeros(length(perception.TE),dataParams);
        for j = 1:length(perception.TE)
            perception.Data(j,:) = data(perception.TE(j,6),:);
        end
        perception.Data.contrast(i) = zeros(length(perception.TE)/p.numContrasts,dataParams);
        for i = 1:p.numContrasts
            perception.TE.contrast(i) = perception.TE((perception.TE(:,3) == p.centerContrast(1)),:);
            perception.Data.contrast(i) = zeros(length(perception.TE)/p.numContrasts,dataParams);
        end
             
    elseif strcmp(variableCondition,'Working Memory') == 1
        workingmem.TE = p.trialEvents((p.trialEvents(:,1)==thisRunsCond),:);
        workingmem.Data = zeros(length(workingmem.TE),dataParams);
        for j = 1:length(workingmem.TE)
            workingmem.Data(j,:) = data(workingmem.TE(j,6),:);
        end
    end

    % NEED TO KNOW NUMBER OF CONTRASTS BEFORE, UNLESS YOU DO A DIFFERENT
    % PAGE OF TRIAL EVENTS AND DATA FOR EACH NEW CONTRAST!
    
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
%     % Seperate trials based off of contrast (within their condition)
%     cond1cont1TE = cond1TE((cond1TE(:,3) == p.centerContrast(1)),:);
%     cond1cont2TE = cond1TE((cond1TE(:,3) == p.centerContrast(2)),:);
%     cond1cont3TE = cond1TE((cond1TE(:,3) == p.centerContrast(3)),:);
%     cond1cont4TE = cond1TE((cond1TE(:,3) == p.centerContrast(4)),:);
%     cond1cont5TE = cond1TE((cond1TE(:,3) == p.centerContrast(5)),:);
%     cond1orgTE = [cond1cont1TE; cond1cont2TE; cond1cont3TE; cond1cont4TE; cond1cont5TE];
%     
%     cond2cont1TE = cond2TE((cond2TE(:,3) == p.centerContrast(1)),:);
%     cond2cont2TE = cond2TE((cond2TE(:,3) == p.centerContrast(2)),:);
%     cond2cont3TE = cond2TE((cond2TE(:,3) == p.centerContrast(3)),:);
%     cond2cont4TE = cond2TE((cond2TE(:,3) == p.centerContrast(4)),:);
%     cond2cont5TE = cond2TE((cond2TE(:,3) == p.centerContrast(5)),:);
%     cond2orgTE = [cond2cont1TE; cond2cont2TE; cond2cont3TE; cond2cont4TE; cond2cont5TE]; 
%     
%      % Seperate data based off of contrast
%     cond1cont1Data = zeros(length(cond1cont1TE),dataParams);
%     cond1cont2Data = zeros(length(cond1cont2TE),dataParams);
%     cond1cont3Data = zeros(length(cond1cont3TE),dataParams);
%     cond1cont4Data = zeros(length(cond1cont4TE),dataParams);
%     cond1cont5Data = zeros(length(cond1cont5TE),dataParams);
%     
%     cond2cont1Data = zeros(length(cond2cont1TE),dataParams);
%     cond2cont2Data = zeros(length(cond2cont2TE),dataParams);
%     cond2cont3Data = zeros(length(cond2cont3TE),dataParams);
%     cond2cont4Data = zeros(length(cond2cont4TE),dataParams);
%     cond2cont5Data = zeros(length(cond2cont5TE),dataParams);
%     for j = 1:length(cond1cont1TE)
%         cond1cont1Data(j,:) = data(cond1cont1TE(j,6),:);
%         cond1cont2Data(j,:) = data(cond1cont2TE(j,6),:);
%         cond1cont3Data(j,:) = data(cond1cont3TE(j,6),:);
%         cond1cont4Data(j,:) = data(cond1cont4TE(j,6),:);
%         cond1cont5Data(j,:) = data(cond1cont5TE(j,6),:);
%         
%         cond2cont1Data(j,:) = data(cond2cont1TE(j,6),:);
%         cond2cont2Data(j,:) = data(cond2cont2TE(j,6),:);
%         cond2cont3Data(j,:) = data(cond2cont3TE(j,6),:);
%         cond2cont4Data(j,:) = data(cond2cont4TE(j,6),:);
%         cond2cont5Data(j,:) = data(cond2cont5TE(j,6),:);
%     end
%      cond1orgData = [cond1cont1Data; cond1cont2Data; cond1cont3Data; cond1cont4Data; cond1cont5Data];
%      cond2orgData = [cond2cont1Data; cond2cont2Data; cond2cont3Data; cond2cont4Data; cond2cont5Data];
%     
%     
% 
% %     theData(i).p.cond1Data = cond1Data; 
% %     theData(i).p.cond2Data = cond2Data;

    
    
 end
        

    

