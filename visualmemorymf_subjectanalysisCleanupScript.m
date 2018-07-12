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
if sum(strcmp(allP{1,1}.experiment,{'test','test_HC'})) == 1
    subjectCondSchedule = [1 1 1 1]; % Fixed to perception for test trials.
else
condIndex = find(strcmp(visualmemory_subjectsRan,allP{1,1}.subject));
    if condIndex > 24
        condIndex = condIndex - 24; %The condition order resets after 24, this matches to the reset.
    end
subjectCondSchedule = visualmemory_condition_order(condIndex,:); %Gives the current condition schedule, indexes to the row we are on, columns 1-4 represent the condition for each run.
end

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
    
    % Seperate Baseline Condition, loop for contrast seperation
    baseline.TE = p.trialEvents((p.trialEvents(:,1)==3),:);
    baseline.Data = zeros(length(baseline.TE),dataParams);
    % Loop that adds baseline data to a seperate matrix
    for j = 1:length(baseline.TE)
        baseline.Data(j,:) = data(baseline.TE(j,6),:);
    end
    for i = 1:p.numContrasts
        baseline.contTE(:,:,i) = baseline.TE(baseline.TE(:,3) == p.centerContrast(i),:);
        baseline.contData(:,:,i) = baseline.Data((baseline.TE(:,3) == p.centerContrast(i)),:);
    end
    %Organizes the contrast seperation trial events and data based on
    %contrast
    %%%%%%% COME BACK TO THIS, ITS CONCATENATING
    %%%%%%% INCORRECTLY also needs to be added to each
    %%%%%%% condition
    for i = 1:p.numContrasts
        if isfield(baseline,'orgTE') == 1
            baseline.orgTE = vertcat(baseline.orgTE,baseline.contTE(:,:,i));
        else
            baseline.orgTE = baseline.contTE(:,:,i);
        end
    end
    
    % Perception versus working memory loop
    if strcmp(variableCondition,'Perception') == 1
        perception.TE = p.trialEvents((p.trialEvents(:,1)==thisRunsCond),:);
        perception.Data = zeros(length(perception.TE),dataParams);
        % Loop that adds perception data to a seperate matrix
        for j = 1:length(perception.TE)
            perception.Data(j,:) = data(perception.TE(j,6),:);
        end
        % Loop that adds TE/data to a 3D matrix, pages seperate by contrast
        for i = 1:p.numContrasts
            perception.contTE(:,:,i) = perception.TE(perception.TE(:,3) == p.centerContrast(i),:);
            perception.contData(:,:,i) = perception.Data((perception.TE(:,3) == p.centerContrast(i)),:);
        end
        
        % MEANS
        % Estimated Contrast -- %% WHY IS THIS ALSO COMING OUT TO 400!!!!!,
        % fix whatever bug this is also
        for i = 1:p.numContrasts
            if exist('percepMeanVec','var') == 1
                percepMeanVec = [percepMeanVec mean(perception.contData(:,4,i)) * ones(1,length(perception.contData(:,4,i)))];
            else
                percepMeanVec = mean(perception.contData(:,4,i)) * ones(1,length(perception.contData(:,4,i)));
            end
        end
%         
%         for i = 1:p.numContrasts
%             if i = 1
%             elseif any(i = [2:5]) == 1
%             end
%         end
        % Concatenation loop to organize the contrast seperated data - COME
        % BACK TO THIS!!!!
          
        
%          % Estimated Contrast Means
%         cond1cont1mean = ones([1 length(cond1cont1Data)])*mean(cond1cont1Data(:,4));
%         cond1cont2mean = ones([1 length(cond1cont2Data)])*mean(cond1cont2Data(:,4));
%         cond1cont3mean = ones([1 length(cond1cont3Data)])*mean(cond1cont3Data(:,4));
%         cond1cont4mean = ones([1 length(cond1cont4Data)])*mean(cond1cont4Data(:,4));
%         cond1cont5mean = ones([1 length(cond1cont5Data)])*mean(cond1cont5Data(:,4));
%         cond1meanVec = [cond1cont1mean  cond1cont2mean cond1cont3mean cond1cont4mean cond1cont5mean];
        
        
        
        
        
    elseif strcmp(variableCondition,'Working Memory') == 1
        workingmem.TE = p.trialEvents((p.trialEvents(:,1)==thisRunsCond),:);
        workingmem.Data = zeros(length(workingmem.TE),dataParams);
        % Loop that adds wm  data to a seperate matrix
        for j = 1:length(workingmem.TE)
            workingmem.Data(j,:) = data(workingmem.TE(j,6),:);
        end
        % Loop that adds contrast TE/data to a 3D matrix, pages sep by contrast
        for i = 1:p.numContrasts
            workingmem.contTE(:,:,i) = workingmem.TE(workingmem.TE(:,3) == p.centerContrast(i),:);
            workingmem.contData(:,:,i) = workingmem.Data((workingmem.TE(:,3) == p.centerContrast(i)),:);
        end 
    end
 end
% for i = 1:p.numContrasts
%     if exist('C','var') == 1
%         C = vertcat(C,baseline.contTE(:,:,i))
%     else
%         C = (vertcat(baseline.contTE(:,:,i)))
%     end
% end