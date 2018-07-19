%%  SINGLE SUBJECT ANALYSIS %%
% This script is meant to loop through all trials with either hardcoded
% test runs, test runs, or actual experiment runs for one single subject.
% It will show or print meaningful data for all runs completed for one
% subject.

%% SETUP %%
% Preliminary data loading and setup %
clear;
close all;
expDir = '/Users/juliaschwartz/Desktop/visualmemorymf';
dataDir = 'data_master';
allP.experiment = 'exp';
allP.subject = 'SL';
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
plotVar = 1; %Set equal to 1 to display plots, equal to 0 to not display plots.
printVar = 1; %Set equal to 0 to not print information, equal to 1 to print information.

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
end

%Finding subject name and indexing to condition order
if sum(strcmp(allP{1,1}.experiment,{'test','test_HC'})) == 1
    subjectCondSchedule = [1 1 1 1]; % Fixed to perception for test trials.
else
    condIndex = find(strcmp(visualmemory_subjectsRan,allP{1,1}.subject));
    if condIndex > 24
        condIndex = condIndex - 24; %The condition order resets after 24, this matches to the reset.
    end
    subjectCondSchedule = visualmemory_condition_order(condIndex,:); %Gives the current condition schedule, indexes to the row we are on, columns 1-4 represent the condition for each run
end

%Save out Relevant Information, 3 for bl+percep+wm
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
 if nTrials == 1
     subjectCondSchedule = allP{1,1}.trialSchedule(1); 
 elseif nTrials < 4
     subjectCondSchedule = subjectCondSchedule(1:nTrials);
 end
 
%% MAIN FOR LOOP: NUMBER OF TRIALS %%

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
    
    p.centerContrast = (10.^linspace(log10(p.minContrast),log10(p.maxContrast),p.numContrasts));
    
    % Add on Trial Number to end of p.trialEvents/data
    trialNum = (1:length(p.trialEvents))';
    p.trialEvents = [p.trialEvents trialNum];
    [dataTrials,dataParams] = size(data);
    
    %% BASELINE CONDITION %%
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
    baseline.orgTE = sortrows(baseline.TE,[3 6]);
    sort = baseline.TE(:,6);
    contsort = baseline.TE(:,3);
    baseline.Data = [baseline.Data contsort sort];
    baseline.orgData = sortrows(baseline.Data, [8 9]);
    % remove contrast and order number from data
    baseline.Data = baseline.Data(:,(1:7));
    baseline.orgData = baseline.orgData(:,(1:7));
    
    % MEANS
    % Estimated Contrast/Contrast Diff/ Location Diff
    
    baseline.EstContMeanVec = ones(1,p.numContrasts);
    baseline.ContDiffMeanVec = ones(1,p.numContrasts);
    baseline.LocDiffMeanVec = ones(1,p.numContrasts);
    for i = 1:p.numContrasts
        % BASELINE IS PAGE 3
        subject.avgEstContrast(nRun,i,3) = mean(baseline.contData(:,4,i));
        subject.avgDiffContrast(nRun,i,3) = abs(mean(baseline.contData(:,5,i)));
        subject.avgDiffLoc(nRun,i,3) = mean(baseline.contData(:,2,i));
        baseline.EstContMeanVec(i) = mean(baseline.contData(:,4,i));
        baseline.ContDiffMeanVec(i) = abs(mean(baseline.contData(:,5,i)));
        baseline.LocDiffMeanVec(i) = mean(baseline.contData(:,2,i));
    end
    
     %% PERCEPTION VS WORKING MEMORY LOOP %%
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
        %Organizes the contrast seperation trial events and data based on
        %contrast
        perception.orgTE = sortrows(perception.TE,[3 6]);
        Psort = perception.TE(:,6);
        Pcontsort = perception.TE(:,3);
        perception.Data = [perception.Data Pcontsort Psort];
        perception.orgData = sortrows(perception.Data, [8 9]);
        % remove contrast and order number from data
        perception.Data = perception.Data(:,(1:7));
        perception.orgData = perception.orgData(:,(1:7));
    
        % MEANS
        % Estimated Contrast/Contrast Diff/ Location Diff
        perception.EstContMeanVec = ones(1,p.numContrasts);
        perception.ContDiffMeanVec = ones(1,p.numContrasts);
        perception.LocDiffMeanVec = ones(1,p.numContrasts);
        for i = 1:p.numContrasts
            % PERCEPTION IS PAGE 1
            subject.avgEstContrast(nRun,i,1) = mean(perception.contData(:,4,i));
            subject.avgDiffContrast(nRun,i,1) = abs(mean(perception.contData(:,5,i)));
            subject.avgDiffLoc(nRun,i,1) = mean(perception.contData(:,2,i));
            perception.EstContMeanVec(i) = mean(perception.contData(:,4,i));
            perception.ContDiffMeanVec(i) = abs(mean(perception.contData(:,5,i)));
            perception.LocDiffMeanVec(i) = mean(perception.contData(:,2,i));
        end
        
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
        %Organizes the contrast seperation trial events and data based on
        %contrast
        workingmem.orgTE = sortrows(workingmem.TE,[3 6]);
        WMsort = workingmem.TE(:,6);
        WMcontsort = workingmem.TE(:,3);
        workingmem.Data = [workingmem.Data WMcontsort WMsort];
        workingmem.orgData = sortrows(workingmem.Data, [8 9]);
        % remove contrast and order number from data
        workingmem.Data = workingmem.Data(:,(1:7));
        workingmem.orgData = workingmem.orgData(:,(1:7));
        
        % MEANS
        % Estimated Contrast/Contrast Diff/ Location Diff
        workingmem.EstContMeanVec = ones(1,p.numContrasts);
        workingmem.ContDiffMeanVec = ones(1,p.numContrasts);
        workingmem.LocDiffMeanVec = ones(1,p.numContrasts);
        for i = 1:p.numContrasts
            % WORKING MEMORY IS PAGE 2
            subject.avgEstContrast(nRun,i,2) = mean(workingmem.contData(:,4,i));
            subject.avgDiffContrast(nRun,i,2) = abs(mean(workingmem.contData(:,5,i)));
            subject.avgDiffLoc(nRun,i,2) = mean(workingmem.contData(:,2,i));
            workingmem.EstContMeanVec(i) = mean(workingmem.contData(:,4,i));
            workingmem.ContDiffMeanVec(i) = abs(mean(workingmem.contData(:,5,i)));
            workingmem.LocDiffMeanVec(i) = mean(workingmem.contData(:,2,i));
        end
    end
    
    %% Plotting (conditional plot variable must not equal 0 to display) %%
    % The plots within this for loop display location and contrast visuals
    % for each trial.
    if plotVar ~= 0 
        % ESTIMATED CONTRAST PLOTTING %
        figure(nRun)
        set(gcf, 'Name', sprintf('Estimated Contrast Statistics over %i Contrasts for %s versus %s Trials',p.numContrasts,baselineCondition,variableCondition));
        
        % BASELINE %
        % Contrast overview plot %
        subplot(2,p.numContrasts+1,1)
        plot(baseline.orgData(:,4));
        hold on
        plot(baseline.orgTE(:,3),'Linewidth',2);
        ylim([0 1])
        title('BASELINE')
        hold on
        plot(1:length(baseline.Data),repelem(subject.avgEstContrast(nRun,:,3),20),'Linewidth',2);
        legend('Estimated Contrast','Actual Contrast','Avg. Estimated Contrast');
        hold off
        
        % Histogram Plots for each contrast %
        for i = 1:p.numContrasts
            subplot(2,p.numContrasts+1,i+1)
            hist(baseline.contData(:,4,i))
            xlim([0 1])
            ylim([0 8])
            hold on
            line([baseline.contTE(1,3,i) baseline.contTE(1,3,i)],ylim,'Linewidth',1.75,'Color','r')
            hold on
            line([subject.avgEstContrast(nRun,i,3) subject.avgEstContrast(nRun,i,3)],ylim,'Linewidth',1.75,'Color','g');
            hold off
            title(sprintf('Histogram for %.2f Trials',baseline.contTE(1,3,i)))
            legend('Est. Contrast Bins','Actual Contrast','Avg. Est Contrast')
        end
        
        % VARIABLE CONDITION CONTRAST PLOTTING %
        if sum(strcmp(variableCondition,'Perception')) == 1
            % PERCEPTION %
            % Contrast overview plot %
            subplot(2,p.numContrasts+1,p.numContrasts+2)
            plot(perception.orgData(:,4));
            hold on
            plot(perception.orgTE(:,3),'Linewidth',2);
            ylim([0 1])
            title('PERCEPTION')
            hold on
            plot(1:length(perception.Data),repelem(subject.avgEstContrast(nRun,:,1),20),'Linewidth',2);
            legend('Estimated Contrast','Actual Contrast','Avg. Estimated Contrast');
            hold off
            
            % Histogram Plots for each contrast %
            for i = 1:p.numContrasts
                subplot(2,p.numContrasts+1,p.numContrasts+2+i)
                hist(perception.contData(:,4,i))
                xlim([0 1])
                ylim([0 8])
                hold on
                line([perception.contTE(1,3,i) perception.contTE(1,3,i)],ylim,'Linewidth',1.75,'Color','r')
                hold on
                line([subject.avgEstContrast(nRun,i,1) subject.avgEstContrast(nRun,i,1)],ylim,'Linewidth',1.75,'Color','g');
                hold off
                title(sprintf('Histogram for %.2f Trials',perception.contTE(1,3,i)))
                legend('Est. Contrast Bins','Actual Contrast','Avg. Est Contrast')
            end
  
        elseif sum(strcmp(variableCondition,'Working Memory')) == 1
            % WORKING MEMORY %
            % Contrast overview plot %
            subplot(2,p.numContrasts+1,p.numContrasts+2)
            plot(workingmem.orgData(:,4));
            hold on
            plot(workingmem.orgTE(:,3),'Linewidth',2);
            ylim([0 1])
            title('WORKING MEMORY')
            hold on
            plot(1:length(workingmem.Data),repelem(subject.avgEstContrast(nRun,:,2),20),'Linewidth',2);
            legend('Estimated Contrast','Actual Contrast','Avg. Estimated Contrast');
            hold off
            
            % Histogram Plots for each contrast %
            for i = 1:p.numContrasts
                subplot(2,p.numContrasts+1,p.numContrasts+2+i)
                hist(workingmem.contData(:,4,i))
                xlim([0 1])
                ylim([0 8])
                hold on
                line([workingmem.contTE(1,3,i) workingmem.contTE(1,3,i)],ylim,'Linewidth',1.75,'Color','r')
                hold on
                line([subject.avgEstContrast(nRun,i,2) subject.avgEstContrast(nRun,i,2)],ylim,'Linewidth',1.75,'Color','g');
                hold off
                title(sprintf('Histogram for %.2f Trials',workingmem.contTE(1,3,i)))
                legend('Est. Contrast Bins','Actual Contrast','Avg. Est Contrast')
            end
        end
        
        % LOCATION DIFFERENCE PLOTTING %
        figure(nRun+nTrials)
        set(gcf, 'Name', sprintf('Estimated Location Statistics over %i Contrasts for %s versus %s Trials',p.numContrasts,baselineCondition,variableCondition));
        
        % BASELINE %
        % Location overview plot %
        subplot(2,p.numContrasts+1,1)
        plot(baseline.orgData(:,2));       
        title('BASELINE')
        hold on
        plot(1:length(baseline.Data),repelem(subject.avgDiffLoc(nRun,:,3),20),'Linewidth',2);
        legend('Location Difference','Avg. Loc. Diff. per Contrast');
        hold off
        
        % Histogram Plots for each contrast %
        for i = 1:p.numContrasts
            subplot(2,p.numContrasts+1,i+1)
            hist(baseline.contData(:,2,i))
            hold on
            line([subject.avgDiffLoc(nRun,i,3) subject.avgDiffLoc(nRun,i,3)],ylim,'Linewidth',1.75,'Color','r');
            hold off
            title(sprintf('Histogram for %.2f Trials',baseline.contTE(1,3,i)))
            legend('Location Difference','Avg. Loc. Diff. Per Contrast')
        end
        
        % VARIABLE CONDITION LOCATION DIFFERENCE %
        if sum(strcmp(variableCondition,'Perception')) == 1
            % PERCEPTION %
            % Location overview plot %
            subplot(2,p.numContrasts+1,p.numContrasts+2)
            plot(perception.orgData(:,2));
            hold on
            title('PERCEPTION')
            plot(1:length(perception.Data),repelem(subject.avgDiffLoc(nRun,:,1),20),'Linewidth',2);
            legend('Location Difference','Avg. Location Difference');
            hold off
            
            % Histogram Plots for each contrast %
            for i = 1:p.numContrasts
                subplot(2,p.numContrasts+1,p.numContrasts+2+i)
                hist(perception.contData(:,2,i))
                hold on
                line([subject.avgDiffLoc(nRun,i,1) subject.avgDiffLoc(nRun,i,1)],ylim,'Linewidth',1.75,'Color','r');
                hold off
                title(sprintf('Histogram for %.2f Trials',perception.contTE(1,3,i)))
                legend('Location Diff. Bins','Avg. Location Diff.')
            end
  
        elseif sum(strcmp(variableCondition,'Working Memory')) == 1
            % WORKING MEMORY %
            % Location Difference overview plot %
            subplot(2,p.numContrasts+1,p.numContrasts+2)
            plot(workingmem.orgData(:,2));
            title('WORKING MEMORY')
            hold on
            plot(1:length(workingmem.Data),repelem(subject.avgDiffLoc(nRun,:,2),20),'Linewidth',2);
            legend('Estimated Contrast','Actual Contrast','Avg. Estimated Contrast');
            hold off
            
            % Histogram Plots of Location Difference for each Contrast %
            for i = 1:p.numContrasts
                subplot(2,p.numContrasts+1,p.numContrasts+2+i)
                hist(workingmem.contData(:,4,i),10)
                xlim([0 1])
                ylim([0 8])
                hold on
                line([subject.avgDiffLoc(nRun,i,2) subject.avgDiffLoc(nRun,i,2)],ylim,'Linewidth',1.75,'Color','g')
                hold off
                title(sprintf('Histogram for %.2f Trials',workingmem.contTE(1,3,i)))
                legend('Est. Contrast Bins','Actual Contrast','Avg. Est Contrast')
            end
        end
    end    

    %% Printing Data (conditional print variable must not equal 0 to display) %%
    if printVar ~= 0
        fprintf('\n\nTRIAL:%i\n',nRun);
        for i = 1:p.numContrasts
        fprintf('  BASELINE:The mean contrast estimation was %.4f at the %.3f contrast level. The abs. difference is %.3f',subject.avgEstContrast(nRun,i,3),p.centerContrast(i),abs(subject.avgEstContrast(nRun,i,3)-p.centerContrast(i)));
            if sum(strcmp(variableCondition,'Perception')) == 1
                fprintf('\nPERCEPTION:The mean contrast estimation was %.4f at the %.3f contrast level. The abs. difference is %.3f\n\n',subject.avgEstContrast(nRun,i,1),p.centerContrast(i),abs(subject.avgEstContrast(nRun,i,1)-p.centerContrast(i)));
            elseif sum(strcmp(variableCondition,'Working Memory')) == 1
                fprintf('\nWORKING MEMORY:The mean contrast estimation was %.4f at the %.3f contrast level. The abs. difference is %.3f\n\n',subject.avgEstContrast(nRun,i,2),p.centerContrast(i),abs(subject.avgEstContrast(nRun,i,2)-p.centerContrast(i)));
            end
        fprintf('  BASELINE:The average location difference was %0.4f at the %0.3f contrast level.\n',subject.avgDiffLoc(nRun,i,3),p.centerContrast(i))
            if sum(strcmp(variableCondition,'Perception')) == 1
                fprintf('PERCEPTION:The average location difference was %0.4f at the %0.3f contrast level.\n\n',subject.avgDiffLoc(nRun,i,1),p.centerContrast(i))
            elseif sum(strcmp(variableCondition,'Working Memory')) == 1
                fprintf('WORKING MEMORY:The average location difference was %0.4f at the %0.3f contrast level.\n\n',subject.avgDiffLoc(nRun,i,2),p.centerContrast(i))
            end
        end
    end
 end
  %% Plotting (conditional plot variable must not equal 0 to display) %%
    % The plots outside of the main for loop display visuals that compare
    % average information from each run.
    if plotVar ~= 0 
        
        % avg Estimated Contrast
        figure(nTrials*2 + 1)
        set(gcf, 'Name', sprintf('Perceived Contrast Versus Center Contrast'));
        loglog(p.centerContrast,mean(subject.avgEstContrast(:,:,3),1),'-o')
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
         % Legend incorporating nans
        if sum(isnan(mean(subject.avgEstContrast(:,:,1),1))) ~= 0 && sum(isnan(mean(subject.avgEstContrast(:,:,2),1))) ~= 0
            legend('Baseline','Log Scale')
        elseif sum(isnan(mean(subject.avgEstContrast(:,:,1),1))) ~= 0 && sum(isnan(mean(subject.avgEstContrast(:,:,2),1))) == 0
            legend('Baseline','Working Memory','Log Scale')
        elseif sum(isnan(mean(subject.avgEstContrast(:,:,1),1))) == 0 && sum(isnan(mean(subject.avgEstContrast(:,:,2),1))) ~= 0
            legend('Baseline','Perception','Log Scale')
        else
            legend('Baseline','Perception','Working Memory','Log Scale')
        end
        hold off 
        
        
        % avg Contrast Difference
        figure(nTrials*2 + 2)
        set(gcf, 'Name', ('Contrast Difference versus Center Contrast'))
        hold on
        loglog(p.centerContrast,mean(subject.avgDiffContrast(:,:,3),1),'-o')
        hold on
        if isnan(mean(subject.avgDiffContrast(:,:,2),1)) == 0
        loglog(p.centerContrast,mean(subject.avgDiffContrast(:,:,1),1),'-o') 
        end
        hold on
        if isnan(mean(subject.avgDiffContrast(:,:,2),1)) == 0
        loglog(p.centerContrast,mean(subject.avgDiffContrast(:,:,2),1),'-o') 
        end
        xlabel('Center Contrast')
        ylabel('Difference in Contrast')
         % Legend incorporating nans
        if sum(isnan(mean(subject.avgDiffContrast(:,:,1),1))) ~= 0 && sum(isnan(mean(subject.avgDiffContrast(:,:,2),1))) ~= 0
            legend('Baseline')
        elseif sum(isnan(mean(subject.avgDiffContrast(:,:,1),1))) ~= 0 && sum(isnan(mean(subject.avgDiffContrast(:,:,2),1))) == 0
            legend('Baseline','Working Memory')
        elseif sum(isnan(mean(subject.avgDiffContrast(:,:,1),1))) == 0 && sum(isnan(mean(subject.avgDiffContrast(:,:,2),1))) ~= 0
            legend('Baseline','Perception')
        else
            legend('Baseline','Perception','Working Memory')
        end
        hold off

        % avg Location Difference
        figure(nTrials*2 + 3)
        set(gcf, 'Name',('Location Difference versus Center Contrast'));
        hold on
        loglog(p.centerContrast,mean(subject.avgDiffLoc(:,:,3),1),'-o') 
        hold on
        if isnan(mean(subject.avgDiffLoc(:,:,1),1)) == 0
        loglog(p.centerContrast,mean(subject.avgDiffLoc(:,:,1),1),'-o') 
        end
        hold on
        if isnan(mean(subject.avgDiffLoc(:,:,2),1)) == 0
        loglog(p.centerContrast,mean(subject.avgDiffLoc(:,:,2),1),'-o')
        end
        hold off
        xlabel('Center Contrast')
        ylabel('Difference in Location')
         % Legend incorporating nans
        if sum(isnan(mean(subject.avgDiffLoc(:,:,1),1))) ~= 0 && sum(isnan(mean(subject.avgDiffLoc(:,:,2),1))) ~= 0
            legend('Baseline')
        elseif sum(isnan(mean(subject.avgDiffLoc(:,:,1),1))) ~= 0 && sum(isnan(mean(subject.avgDiffLoc(:,:,2),1))) == 0
            legend('Baseline','Working Memory')
        elseif sum(isnan(mean(subject.avgDiffLoc(:,:,1),1))) == 0 && sum(isnan(mean(subject.avgDiffLoc(:,:,2),1))) ~= 0
            legend('Baseline','Perception')
        else
            legend('Baseline','Perception','Working Memory')
        end
        hold off 
    end
    
%% SAVE OUT AVERAGES (AVG OVER TRIALS RAN FOR A SINGLE SUBJECT) %%
subject.meanEstContPerception = mean(subject.avgEstContrast(:,:,1),1);
subject.meanEstContWorkingMemory = mean(subject.avgEstContrast(:,:,2),1);
subject.meanEstContBaseline = mean(subject.avgEstContrast(:,:,3),1);

subject.meanDiffContPerception = mean(subject.avgDiffContrast(:,:,1),1);
subject.meanDiffContWorkingMemory = mean(subject.avgDiffContrast(:,:,2),1);
subject.meanDiffContBaseline = mean(subject.avgDiffContrast(:,:,3),1);

subject.meanDiffLocPerception = mean(subject.avgDiffLoc(:,:,1),1);
subject.meanDiffLocWorkingMemory = mean(subject.avgDiffLoc(:,:,2),1);
subject.meanDiffLocBaseline = mean(subject.avgDiffLoc(:,:,3),1);
    
%% SAVE SUBJECT STRUCTURE %%
cd(dataDir)
save(['data_visualmemorymf_' p.experiment '_' p.subject '.mat'], 'subject','theData')

      