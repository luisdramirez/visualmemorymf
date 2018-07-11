%%  MASTER ANALYSIS %%
% This script is meant to loop through all trials with either hardcoded
% test runs, test runs, or actual experiment runs. It will show represent
% meaningful data within one single run, or if there are multiple runs for
% one single subject/experiemtn matching, it will loop through to collect
% meaningful results between data sets.

% Preliminary data loading and setup %
clear all;
close all;
expDir = 'visualmemorymf';
dataDir = 'data_master';
p.experiment = 'test';
p.subject = 'JS';
cd(dataDir)

if exist(['data_visualmemorymf_' p.experiment '_' p.subject '.mat'],'file') ~= 0
    load(['data_visualmemorymf_' p.experiment '_' p.subject '.mat']); %works for HC, test, and regular trials
else
    error('data file does not exist')
end

%% LOOP FOR ALL TRIALS POSSIBLE %%

[fields, runsCompleted] = size(theData);

%Save out Relevant Information
subject.avgEstContrast = nan(runsCompleted,theData(1).p.numContrasts,3);
subject.avgDiffContrast = nan(runsCompleted,theData(1).p.numContrasts,3);
subject.avgDiffLoc = nan(runsCompleted,theData(1).p.numContrasts,3);

for i = 1:runsCompleted
    theDataCurrent = theData(i); %index to trial number
    p = theDataCurrent.p; data = theDataCurrent.data; t = theDataCurrent.t;
    data = cell2mat(struct2cell(data));
    data = data';
    [trials,params] = size(p.trialEvents);
    if p.numContrasts == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %Test_HC and 1 Contrast Runs%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if any(p.trialEvents(:,1) == 1) % Perception trials, condition # = 1
            cond1Results = p.trialEvents(p.trialEvents(:,1)==1,:);
            cond1Data = data(p.trialEvents(:,1)==1,:);
            cond1Name = 'Perception';
        end
        if any(p.trialEvents(:,1) == 2) % Working Memory trials, condition # = 2
            if isempty(cond1Data) == 1  &&  isempty(firstcondResults) == 1 && ~strcmp(cond1Name,'Perception')
                cond1Results = p.trialEvents(p.trialEvents(:,1)==2,:);
                cond1Data = data(p.trialEvents(:,1)==2,:);
                cond1Name = 'Working Memory';
            else
                cond2Results = p.trialEvents(p.trialEvents(:,1)==2,:);
                cond2Data = data(p.trialEvents(:,1)==2,:);
                cond2Name = 'Working Memory';
            end
        end
        if any(p.trialEvents(:,1) == 3)  % Baseline trials, condition = 3
            if isempty(cond1Data) == 1  &&  isempty(firstcondResults) == 1 && ~strcmp(cond1Name,'Perception') && ~strcmp(cond1Name,'Working Memory')
                cond1Results = p.trialEvents(p.trialEvents(:,1)==3,:);
                cond1Data = data(p.trialEvents(:,1)==3,:);
                cond1Name = 'Baseline';
            else
                cond2Results = p.trialEvents(p.trialEvents(:,1)==3,:);
                cond2Data = data(p.trialEvents(:,1)==3,:);
                cond2Name = 'Baseline';
            end
        end
        %%  First Condition Estimated Contrast vs. Actual Contrast  %%

        % Compare Perception and Baseline on same y-limit
        minContrast = min([min(cond1Data(:,4)) min(cond2Data(:,4))]);
        maxContrast = max([max(cond1Data(:,4)) max(cond2Data(:,4))]);
        %Lower bound
        if minContrast - 0.1 <= 0
            yMin = 0;
        else 
            yMin = minContrast - 0.1;
        end
        %upper bound
        if maxContrast + 0.1 >= 1
            yMax = 1;
        else
            yMax = maxContrast + 0.1;
        end

        figure(i)
        set(gcf, 'Name', sprintf('%s Versus %s Visual Analysis: Contrast',cond1Name,cond2Name));
        subplot(2,3,1)
        plot(cond1Data(:,4)); %Estimated Contrast
        hold on
        fit = polyfit(1:length(cond1Data(:,4)),cond1Data(:,4)',1);
        plot(1:length(cond1Results(:,3)),cond1Results(:,3));
        %%%% LOG FIT %%%% ?!?
        plot(polyval(fit,1:length(cond1Data(:,4)))); %Line of best fit for the Estimated Contrast
        xlim([1 length(cond1Data(:,4))]);
        ylim([yMin yMax]);
        title(sprintf('%s Contrast Compared to Actual Contrast',cond1Name))
        legend('Subject Contrast Estimation','Actual Contrast', 'Average Contrast Estimation')
        xlabel('Trial Number')
        ylabel('Contrast')

        %%  First Condition Estimated Contrast vs. Actual Contrast  %%
        subplot(2,3,4)
        plot(cond2Data(:,4)); %Estimated Contrast
        hold on
        plot(1:length(cond2Results(:,3)),cond2Results(:,3));
        fit = polyfit(1:length(cond2Data(:,4)),cond2Data(:,4)',1);
        plot(polyval(fit,1:length(cond2Data(:,4)))); %Line of best fit for the Estimated Contrast
        xlim([1 length(cond2Data(:,4))]);
        ylim([yMin yMax]);
        title(sprintf('%s Contrast Estimation Compared to Actual Contrast',cond2Name));
        legend('Subject Contrast Estimation','Actual Contrast', 'Average Contrast Estimation')
        xlabel('Trial Number')
        ylabel('Contrast')

      %% Histogram of First Condition Responses %%
        subplot(2,3,2)
        [cond1bins1,edges1] = histcounts(cond1Data(:,4),20);
        [cond1bins2,edges2] = histcounts(cond2Data(:,4),20);
        histmax = max([max(cond1bins1) max(cond1bins2)]);
        pHist = histogram(cond1Data(:,4),20);
        hold on
        xlim([0 1]);
        ylim([0 histmax]);
        title(sprintf('Histogram of %s Contrast Estimations',cond1Name));
        line([p.centerContrast p.centerContrast],ylim, 'Linewidth',1,'Color','r');
        lgd = legend('Bins','Actual Contrast');
        lgd.Location = 'northwest';
        xlabel('Contrast')
        ylabel('Frequency of response')

     %% Histogram of Second Condition Responses %%
        subplot(2,3,5)
        secondCondHist = histogram(cond2Data(:,4),20);
        xlim([0 1]);
        ylim([0 histmax+1]);
        title(sprintf('Histogram of %s Contrast Estimations',cond2Name));
        hold on
        line([p.centerContrast p.centerContrast],ylim, 'Linewidth',1,'Color','r');
        lgd = legend('Bins','Actual Contrast');
        lgd.Location = 'northwest';
        xlabel('Contrast')
        ylabel('Frequency of response')

    %% Misc. Statistical Information %%

        %Means & Standard Deviations of estimated contrasts
        cond1Std = std(cond1Data(:,4));
        cond1Mean = mean2(cond1Data(:,4));
        cond2Std =std(cond2Data(:,4));
        cond2Mean = mean2(cond2Data(:,4));
        fprintf('\nThe average response for %s trials was a contrast of %.4f, with a standard deviation of %.4f.\nThe absolute difference between the average and actual contrast is %f\n\n',cond1Name,cond1Mean, cond1Std,abs(cond1Mean-p.centerContrast));
        fprintf('The average response for %s trials was a contrast of %.4f, with a standard deviation of %.4f.\nThe absolute difference between the average and actual contrast is %f\n\n',cond2Name,cond2Mean, cond2Std, abs(cond2Mean-p.centerContrast));

        %Percent Error of Contrast Estimations
        cond1PE = zeros(size(1:length(cond1Data(:,4))));
        cond2PE = zeros(size(1:length(cond2Data(:,4))));
        for j = 1:length(cond1Data(:,4))
            cond1PE(j) = (abs((cond1Data(j,4)-cond1Results(j,3))/cond1Results(j,3)))*100;
            cond2PE(j) = (abs((cond2Data(j,4)-cond2Results(j,3))/cond2Results(j,3)))*100;
        end
        maxPE = max([max(cond1PE) max(cond2PE)]);
        if maxPE + 5 >= 100
            yLimPE = 100;
        else
            yLimPE = (maxPE + 5);
        end

        % Paired Sample T-test
        [h,p] = ttest(cond1Data(:,5),cond2Data(:,5));
        if h == 0
            fprintf('The reported contrasts for %s and %s were not statistically significant\n\n',cond1Name,cond2Name);
        elseif h == 1
            fprintf('The reported contrasts for %s and %s were statistically significant from each other, with a p value of %.3f\n\n',cond1Name,cond2Name,p);
        end
    %% Percent Error Plot First Condition %%
        subplot(2,3,3)
        plot(cond1PE);
        xlim([1 length(cond1Data)])
        ylim([0 yLimPE])
        title(sprintf('%s Contrast Percent Error',cond1Name))
        hold on
        plot(repmat(mean2(cond1PE),1,length(cond1Data)))
        legend('Percent Error of Each Estimation','Average Percent Error')
        xlabel('Trial Number')
        ylabel('Percent Error')
    %% Percent Error Plot Second Condition %%
        subplot(2,3,6)
        plot(cond2PE);
        xlim([1 length(cond2Data)])
        ylim([0 yLimPE])
        title(sprintf('%s Contrast Percent Error',cond2Name))
        hold on
        plot(repmat(mean2(cond2PE),1,length(cond2Data)))
        legend('Percent Error of Each Estimation','Average Percent Error')
        xlabel('Trial Number')
        ylabel('Percent Error')
        if cond1Mean > cond2Mean
            fprintf('%s had a higher estimated contrast than %s',cond1Name,cond2Name);
            fprintf('with a percent difference of %.3f%\n\n',(cond1Mean-cond2Mean)*100);
        elseif cond1Mean < cond2Mean
            fprintf('%s had a higher estimated contrast than %s',cond2Name,cond1Name);
            fprintf('with a percent difference of %.3f%\n\n',abs(cond1Mean-cond2Mean)*100);
        end

    
    elseif p.numContrasts == 5 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %Test & 5 Contrast Runs%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Will always be baseline + either working memory or perception...
        %1 - PERCEPTION
        %2 - WORKING MEMORY
        %3 - BASELINE
        
    %Add trial numbers to the end of trialEvents for organizational
    %purposes
    trialNum = (1:length(p.trialEvents))';
    p.trialEvents = [p.trialEvents trialNum];
    thisRunConds = sortrows(unique(p.trialEvents(:,1))); % First index is wm or perception (1 or 2), second is baseline (3)
    
    % Get condition names
    if thisRunConds(1) == 1
        condition1 = 'Perception';
    elseif thisRunConds(1) == 2
        condition1 = 'Working Memory';
    else 
        error('Condition names not set up correctly.')
    end
    if thisRunConds(2) == 3
        condition2 = 'Baseline';
    else
        error('Condition names not set up correctly.')
    end
    
    % Seperate trials based off of condition
    cond1TE = p.trialEvents((p.trialEvents(:,1)==thisRunConds(1)),:);
    cond2TE = p.trialEvents((p.trialEvents(:,1)==thisRunConds(2)),:);
    [dataTrials,dataParams] = size(data);
    cond1Data = zeros(length(cond1TE),dataParams);
    cond2Data = zeros(length(cond2TE),dataParams);
    for j = 1:length(cond1TE)
        cond1Data(j,:) = data(cond1TE(j,6),:);
        cond2Data(j,:) = data(cond2TE(j,6),:);
    end
    theData(i).p.cond1Data = cond1Data; %THIS ISNT WORKING FOR SOME REASON COME BACK
    theData(i).p.cond2Data = cond2Data;
    p.centerContrast = [10.^linspace(log10(p.minContrast),log10(p.maxContrast),p.numContrasts)]; %change if contrasts change
    % Seperate trials based off of contrast (within their condition)
    cond1cont1TE = cond1TE((cond1TE(:,3) == p.centerContrast(1)),:);
    cond1cont2TE = cond1TE((cond1TE(:,3) == p.centerContrast(2)),:);
    cond1cont3TE = cond1TE((cond1TE(:,3) == p.centerContrast(3)),:);
    cond1cont4TE = cond1TE((cond1TE(:,3) == p.centerContrast(4)),:);
    cond1cont5TE = cond1TE((cond1TE(:,3) == p.centerContrast(5)),:);
    cond1orgTE = [cond1cont1TE; cond1cont2TE; cond1cont3TE; cond1cont4TE; cond1cont5TE];
    
    cond2cont1TE = cond2TE((cond2TE(:,3) == p.centerContrast(1)),:);
    cond2cont2TE = cond2TE((cond2TE(:,3) == p.centerContrast(2)),:);
    cond2cont3TE = cond2TE((cond2TE(:,3) == p.centerContrast(3)),:);
    cond2cont4TE = cond2TE((cond2TE(:,3) == p.centerContrast(4)),:);
    cond2cont5TE = cond2TE((cond2TE(:,3) == p.centerContrast(5)),:);
    cond2orgTE = [cond2cont1TE; cond2cont2TE; cond2cont3TE; cond2cont4TE; cond2cont5TE]; 
    
    % Seperate data based off of contrast
    cond1cont1Data = zeros(length(cond1cont1TE),dataParams);
    cond1cont2Data = zeros(length(cond1cont2TE),dataParams);
    cond1cont3Data = zeros(length(cond1cont3TE),dataParams);
    cond1cont4Data = zeros(length(cond1cont4TE),dataParams);
    cond1cont5Data = zeros(length(cond1cont5TE),dataParams);
    
    cond2cont1Data = zeros(length(cond2cont1TE),dataParams);
    cond2cont2Data = zeros(length(cond2cont2TE),dataParams);
    cond2cont3Data = zeros(length(cond2cont3TE),dataParams);
    cond2cont4Data = zeros(length(cond2cont4TE),dataParams);
    cond2cont5Data = zeros(length(cond2cont5TE),dataParams);
    for j = 1:length(cond1cont1TE)
        cond1cont1Data(j,:) = data(cond1cont1TE(j,6),:);
        cond1cont2Data(j,:) = data(cond1cont2TE(j,6),:);
        cond1cont3Data(j,:) = data(cond1cont3TE(j,6),:);
        cond1cont4Data(j,:) = data(cond1cont4TE(j,6),:);
        cond1cont5Data(j,:) = data(cond1cont5TE(j,6),:);
        
        cond2cont1Data(j,:) = data(cond2cont1TE(j,6),:);
        cond2cont2Data(j,:) = data(cond2cont2TE(j,6),:);
        cond2cont3Data(j,:) = data(cond2cont3TE(j,6),:);
        cond2cont4Data(j,:) = data(cond2cont4TE(j,6),:);
        cond2cont5Data(j,:) = data(cond2cont5TE(j,6),:);
    end
     cond1orgData = [cond1cont1Data; cond1cont2Data; cond1cont3Data; cond1cont4Data; cond1cont5Data];
     cond2orgData = [cond2cont1Data; cond2cont2Data; cond2cont3Data; cond2cont4Data; cond2cont5Data];
    
    % Estimated Contrast Means
    cond1cont1mean = ones([1 length(cond1cont1Data)])*mean(cond1cont1Data(:,4));
    cond1cont2mean = ones([1 length(cond1cont2Data)])*mean(cond1cont2Data(:,4));
    cond1cont3mean = ones([1 length(cond1cont3Data)])*mean(cond1cont3Data(:,4));
    cond1cont4mean = ones([1 length(cond1cont4Data)])*mean(cond1cont4Data(:,4));
    cond1cont5mean = ones([1 length(cond1cont5Data)])*mean(cond1cont5Data(:,4));
    cond1meanVec = [cond1cont1mean  cond1cont2mean cond1cont3mean cond1cont4mean cond1cont5mean];
    
    cond2cont1mean = ones([1 length(cond2cont1Data)])*mean(cond2cont1Data(:,4));
    cond2cont2mean = ones([1 length(cond2cont2Data)])*mean(cond2cont2Data(:,4));
    cond2cont3mean = ones([1 length(cond2cont3Data)])*mean(cond2cont3Data(:,4));
    cond2cont4mean = ones([1 length(cond2cont4Data)])*mean(cond2cont4Data(:,4));
    cond2cont5mean = ones([1 length(cond2cont5Data)])*mean(cond2cont5Data(:,4));
    cond2meanVec = [cond2cont1mean  cond2cont2mean cond2cont3mean cond2cont4mean cond2cont5mean];
    
    %Estimated Location Means
   
    cond1cont1LocMean = ones([1 length(cond1cont1Data)])*mean(cond1cont1Data(:,2));
    cond1cont2LocMean = ones([1 length(cond1cont2Data)])*mean(cond1cont2Data(:,2));
    cond1cont3LocMean = ones([1 length(cond1cont3Data)])*mean(cond1cont3Data(:,2));
    cond1cont4LocMean = ones([1 length(cond1cont4Data)])*mean(cond1cont4Data(:,2));
    cond1cont5LocMean = ones([1 length(cond1cont5Data)])*mean(cond1cont5Data(:,2));
    cond1LocMeanVec = [cond1cont1LocMean  cond1cont2LocMean cond1cont3LocMean cond1cont4LocMean cond1cont5LocMean];
    cond1LocMeanVec = unique(cond1LocMeanVec);
    
    cond2cont1LocMean = ones([1 length(cond2cont1Data)])*mean(cond2cont1Data(:,2));
    cond2cont2LocMean = ones([1 length(cond2cont2Data)])*mean(cond2cont2Data(:,2));
    cond2cont3LocMean = ones([1 length(cond2cont3Data)])*mean(cond2cont3Data(:,2));
    cond2cont4LocMean = ones([1 length(cond2cont4Data)])*mean(cond2cont4Data(:,2));
    cond2cont5LocMean = ones([1 length(cond2cont5Data)])*mean(cond2cont5Data(:,2));
    cond2LocMeanVec = [cond2cont1LocMean  cond2cont2LocMean cond2cont3LocMean cond2cont4LocMean cond2cont5LocMean];
    cond2LocMeanVec = unique(cond2LocMeanVec);
    
    subject.avgDiffLoc(i,:,2) = cond1LocMeanVec; %page 2 is perception, page 3 is wm
    subject.avgDiffLoc(i,:,1) = cond2LocMeanVec; %page 1 is baseline

    % Contrast Difference Means
    cond1cont1ContDiffMean = ones([1 length(cond1cont1Data)])*mean(cond1cont1Data(:,5));
    cond1cont2ContDiffMean = ones([1 length(cond1cont2Data)])*mean(cond1cont2Data(:,5));
    cond1cont3ContDiffMean = ones([1 length(cond1cont3Data)])*mean(cond1cont3Data(:,5));
    cond1cont4ContDiffMean = ones([1 length(cond1cont4Data)])*mean(cond1cont4Data(:,5));
    cond1cont5ContDiffMean = ones([1 length(cond1cont5Data)])*mean(cond1cont5Data(:,5));
    cond1ContDiffMeanVec = abs(unique([cond1cont1ContDiffMean  cond1cont2ContDiffMean cond1cont3ContDiffMean cond1cont4ContDiffMean cond1cont5ContDiffMean]));
    
    cond2cont1ContDiffMean = ones([1 length(cond2cont1Data)])*mean(cond2cont1Data(:,5));
    cond2cont2ContDiffMean = ones([1 length(cond2cont2Data)])*mean(cond2cont2Data(:,5));
    cond2cont3ContDiffMean = ones([1 length(cond2cont3Data)])*mean(cond2cont3Data(:,5));
    cond2cont4ContDiffMean = ones([1 length(cond2cont4Data)])*mean(cond2cont4Data(:,5));
    cond2cont5ContDiffMean = ones([1 length(cond2cont5Data)])*mean(cond2cont5Data(:,5));
    cond2ContDiffMeanVec = abs(unique([cond2cont1ContDiffMean  cond2cont2ContDiffMean cond2cont3ContDiffMean cond2cont4ContDiffMean cond2cont5ContDiffMean]));
    
    subject.avgContDiff(i,:,2) = cond1ContDiffMeanVec;%page 2 is perception, page 3 is wm
    subject.avgContDiff(i,:,1) = cond2ContDiffMeanVec;%page 1 is baseline
    
     %% Contrast Plotting %%
    
    %%% Contrast Graph Overview %%%
    figure(i)
    set(gcf, 'Name', sprintf('Contrast Statistics over 5 Contrasts for %s versus %s Trials',condition1,condition2));
    subplot(2,6,1)
    
    %Cond 1 - WM/P
    plot(cond1orgData(:,4));
    hold on
    plot(cond1orgTE(:,3));
    ylim([0 1]);
    title(sprintf('Estimated Vs. Actual Contrast for %s',condition1));
    hold on
    plot(1:length(cond1meanVec),cond1meanVec);
    legend('Estimated Contrast','Actual Contrast','Avg. Estimated Contrast per Level');
    
    %Cond 2 - BL
    subplot(2,6,7)
    plot(cond2orgData(:,4));
    hold on
    plot(cond2orgTE(:,3));
    ylim([0 1]);
    title(sprintf('Estimated Vs. Actual Contrast for %s',condition2));
    hold on
    plot(1:length(cond2meanVec),cond2meanVec);
    legend('Estimated Contrast','Actual Contrast','Avg. Estimated Contrast per Level');
    
    % Histogram axis limits
    [cond1bins1,~] = histcounts(cond1cont1Data(:,4),10);
    [cond1bins2,~] = histcounts(cond1cont2Data(:,4),10);
    [cond1bins3,~] = histcounts(cond1cont3Data(:,4),10);
    [cond1bins4,~] = histcounts(cond1cont4Data(:,4),10);
    [cond1bins5,~] = histcounts(cond1cont5Data(:,4),10);
    [cond2bins1,~] = histcounts(cond2cont1Data(:,4),10);
    [cond2bins2,~] = histcounts(cond2cont2Data(:,4),10);
    [cond2bins3,~] = histcounts(cond2cont3Data(:,4),10);
    [cond2bins4,~] = histcounts(cond2cont4Data(:,4),10);
    [cond2bins5,~] = histcounts(cond2cont5Data(:,4),10);
    histmax = max([max(cond1bins1) max(cond1bins2) max(cond1bins3) max(cond1bins4) max(cond1bins5) max(cond2bins1) max(cond2bins2) max(cond2bins3) max(cond2bins4) max(cond2bins5)]);
    
    %%% Condition 1 Contrast Histograms %%%
    
    subplot (2,6,2)
    hist(cond1cont1Data(:,4),10);
    hold on
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    line([cond1cont1TE(:,3) cond1cont1TE(:,3)],ylim, 'Linewidth',1,'Color','r');
    hold on
    line([unique(cond1cont1mean) unique(cond1cont1mean)],ylim, 'Linewidth',3,'Color','g');
    hold off
    title(sprintf('Hist of Estimated Contrast for %.2f Contrastin %s Trials ',cond1cont1TE(1,3),condition1));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    
    subplot (2,6,3)
    hist(cond1cont2Data(:,4),10);
    hold on
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    line([cond1cont2TE(:,3) cond1cont2TE(:,3)],ylim, 'Linewidth',1,'Color','r');
    hold on
    line([unique(cond1cont2mean) unique(cond1cont2mean)],ylim, 'Linewidth',3,'Color','g');
    hold off
    title(sprintf('Hist of Estimated Contrast for %.2f Contrast in %s Trials',cond1cont2TE(1,3),condition1));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    
    subplot (2,6,4)
    hist(cond1cont3Data(:,4),10);
    hold on
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    line([cond1cont3TE(:,3) cond1cont3TE(:,3)],ylim, 'Linewidth',1,'Color','r');
    hold on
    line([unique(cond1cont3mean) unique(cond1cont3mean)],ylim, 'Linewidth',3,'Color','g');
    hold off
    title(sprintf('Hist of Estimated Contrast for %.2f Contrast in %s Trials ',cond1cont3TE(1,3),condition1));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    
    subplot (2,6,5)
    hist(cond1cont4Data(:,4),10);
    hold on
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    line([cond1cont4TE(:,3) cond1cont4TE(:,3)],ylim, 'Linewidth',1,'Color','r');
    hold on
    line([unique(cond1cont4mean) unique(cond1cont4mean)],ylim, 'Linewidth',3,'Color','g');
    hold off
    title(sprintf('Hist of Estimated Contrast for %.2f Contrast in %s Trials ',cond1cont4TE(1,3),condition1));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    
    subplot (2,6,6)
    hist(cond1cont5Data(:,4),10);
    hold on
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    line([cond1cont5TE(:,3) cond1cont5TE(:,3)],ylim, 'Linewidth',1,'Color','r');
    hold on
    line([unique(cond1cont5mean) unique(cond1cont5mean)],ylim, 'Linewidth',3,'Color','g');
    hold off
    title(sprintf('Hist of Estimated Contrast for %.2f Contrastin %s Trials ',cond1cont5TE(1,3),condition1));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    
    %%% Condition 2 Contrast Graphs %%%
    
    subplot (2,6,8)
    hist(cond2cont1Data(:,4),10);
    hold on
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    line([cond2cont1TE(:,3) cond2cont1TE(:,3)],ylim, 'Linewidth',1,'Color','r');
    hold on
    line([unique(cond2cont1mean) unique(cond2cont1mean)],ylim, 'Linewidth',3,'Color','g');
    hold off
    title(sprintf('Hist of Estimated Contrast for %.2f Contrastin %s Trials ',cond2cont1TE(1,3),condition2));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    
    subplot (2,6,9)
    hist(cond2cont2Data(:,4),10);
    hold on
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    line([cond2cont2TE(:,3) cond2cont2TE(:,3)],ylim, 'Linewidth',1,'Color','r');
    hold on
    line([unique(cond2cont2mean) unique(cond2cont2mean)],ylim, 'Linewidth',3,'Color','g');
    hold off
    title(sprintf('Hist of Estimated Contrast for %.2f Contrast in %s Trials',cond2cont2TE(1,3),condition2));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    
    subplot (2,6,10)
    hist(cond2cont3Data(:,4),10);
    hold on
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    line([cond2cont3TE(:,3) cond2cont3TE(:,3)],ylim, 'Linewidth',1,'Color','r');
    hold on
    line([unique(cond2cont3mean) unique(cond2cont3mean)],ylim, 'Linewidth',3,'Color','g');
    hold off
    title(sprintf('Hist of Estimated Contrast for %.2f Contrast in %s Trials',cond2cont3TE(1,3),condition2));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    
    subplot (2,6,11)
    hist(cond2cont4Data(:,4),10);
    hold on
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    line([cond2cont4TE(:,3) cond2cont4TE(:,3)],ylim, 'Linewidth',1,'Color','r');
    hold on
    line([unique(cond2cont4mean) unique(cond2cont4mean)],ylim, 'Linewidth',3,'Color','g');
    hold off
    title(sprintf('Hist of Estimated Contrast for %.2f Contrast in %s Trials',cond2cont4TE(1,3),condition2));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    
    subplot (2,6,12)
    hist(cond2cont5Data(:,4),10);
    hold on
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    line([cond2cont5TE(:,3) cond2cont5TE(:,3)],ylim, 'Linewidth',1,'Color','r');
    hold on
    line([unique(cond2cont5mean) unique(cond2cont5mean)],ylim, 'Linewidth',3,'Color','g');
    hold off
    title(sprintf('Hist of Estimated Contrast for %.2f Contrast in %s Trials',cond2cont5TE(1,3),condition2));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    
    % Printing Data %
    fprintf('\nTRIAL %d\n',i)
    
    fprintf('\nContrast Stats for %s Trials:\n',condition1)
    fprintf('\nThe mean of responses for contrast %f was %f. The absolute difference is: %f',cond1cont1TE(1,3),cond1cont1mean(1), abs(cond1cont1TE(1,3)-cond1cont1mean(1)));
    fprintf('\nThe mean of responses for contrast %f was %f. The absolute difference is: %f',cond1cont2TE(1,3),cond1cont2mean(1), abs(cond1cont2TE(1,3)-cond1cont2mean(1)));
    fprintf('\nThe mean of responses for contrast %f was %f. The absolute difference is: %f',cond1cont3TE(1,3),cond1cont3mean(1), abs(cond1cont3TE(1,3)-cond1cont3mean(1)));
    fprintf('\nThe mean of responses for contrast %f was %f. The absolute difference is: %f',cond1cont4TE(1,3),cond1cont4mean(1), abs(cond1cont4TE(1,3)-cond1cont4mean(1)));
    fprintf('\nThe mean of responses for contrast %f was %f. The absolute difference is: %f\n',cond1cont5TE(1,3),cond1cont5mean(1), abs(cond1cont5TE(1,3)-cond1cont5mean(1)));
    
    fprintf('\nContrast Stats for %s Trials:\n',condition2)
    fprintf('\nThe mean of responses for contrast %f was %f. The absolute difference is: %f',cond2cont1TE(1,3),cond2cont1mean(1), abs(cond2cont1TE(1,3)-cond2cont1mean(1)));
    fprintf('\nThe mean of responses for contrast %f was %f. The absolute difference is: %f',cond2cont2TE(1,3),cond2cont2mean(1), abs(cond2cont2TE(1,3)-cond2cont2mean(1)));
    fprintf('\nThe mean of responses for contrast %f was %f. The absolute difference is: %f',cond2cont3TE(1,3),cond2cont3mean(1), abs(cond2cont3TE(1,3)-cond2cont3mean(1)));
    fprintf('\nThe mean of responses for contrast %f was %f. The absolute difference is: %f',cond2cont4TE(1,3),cond2cont4mean(1), abs(cond2cont4TE(1,3)-cond2cont4mean(1)));
    fprintf('\nThe mean of responses for contrast %f was %f. The absolute difference is: %f\n',cond2cont5TE(1,3),cond2cont5mean(1), abs(cond2cont5TE(1,3)-cond2cont5mean(1)));
    
    
    %% Graph means and estimations on one graph %%  
    figure(1+2*runsCompleted)
    subplot(2,2,i)
    set(gcf, 'Name',('Perceived Contrast versus center contrast'));
    xlabel('Center Contrast')
    ylabel('Perceived Contrast')
    xlim([0.1 0.8])
    ylim([0.1 0.8])
    loglog(p.centerContrast,unique(cond1meanVec),'-o') %perception condition
    hold on
    loglog(p.centerContrast,unique(cond2meanVec),'-o') %baseline condition
    loglog([0.1 0.8],[0.1 0.8], ':')
    legend(sprintf('%s Trials',condition1),sprintf('%s Trials',condition2))
    
    subject.avgEstContrast(i,:,1) = unique(cond2meanVec); % save baseline
    subject.avgEstContrast(i,:,2) = unique(cond1meanVec); % save perception
    
    %% LOCATION %%
    
    %Y Limit for accurate comparisons
    yLocMax = max([max(cond1cont1Data(:,2)) max(cond1cont2Data(:,2)) max(cond1cont3Data(:,2)) max(cond1cont4Data(:,2)) max(cond1cont5Data(:,2)) max(cond2cont1Data(:,2)) max(cond2cont2Data(:,2)) max(cond2cont3Data(:,2)) max(cond2cont4Data(:,2)) max(cond2cont5Data(:,2))])+5;
    figure(i+4)
    set(gcf, 'Name',('Location Statistics over 5 Contrasts'));
    
    % General Location Estimation and Trendline
    subplot(2,6,1)
    plot(cond1orgData(:,2))
    hold on
    ylim([0 yLocMax])
    xlabel('Lowest to Highest Contrast Trials')
    ylabel('Estimated vs. Actual Location Difference (°)')
    title(sprintf('Location Difference for all Contrasts for %s Trials',condition1))
    fit = polyfit(1:length(cond1orgData(:,2)),cond1orgData(:,2)',1);
    plot(polyval(fit,1:length(cond1orgData(:,2)))); 
    hold off
    legend('Abs value of difference in location','Average Location Difference Trendline')
    
    subplot(2,6,7)
    plot(cond2orgData(:,2))
    hold on
    ylim([0 yLocMax])
    xlabel('Lowest to Highest Contrast Trials')
    ylabel('Estimated vs. Actual Location Difference (°)')
    title(sprintf('Location Difference for all Contrasts for %s Trials',condition2))
    fit = polyfit(1:length(cond2orgData(:,2)),cond2orgData(:,2)',1);
    plot(polyval(fit,1:length(cond2orgData(:,2)))); 
    hold off
    legend('Abs value of difference in location','Average Location Difference Trendline')
    
    %% Location Difference Graphs %%
    
    fprintf('\nLocation Stats for %s Trials:\n',condition1)
    %Contrast 1 Location Difference
    subplot(2,6,2)
    plot(cond1cont1Data(:,2))
    xlim([1 length(cond1cont1Data)])
    ylim([0 yLocMax])
    xlabel('Trial Number')
    ylabel('Estimated vs. Actual Location Difference (°)')
    title(sprintf('Location Difference for %.2f Contrast',cond1cont1TE(1,3)))
    hold on
    plot(repmat(mean(cond1cont1Data(:,2)),1,length(cond1cont1Data(:,2))))
    fprintf('\nThe average location difference for a Contrast of %.2f is %.2f\n',cond1cont1TE(1,3),mean(cond1cont1Data(:,2)));
    legend('Abs value of difference in location','Average Location Difference')

     %Contrast 2 Location Difference
    subplot(2,6,3)
    plot(cond1cont2Data(:,2))
    xlim([1 length(cond1cont2Data)])
    ylim([0 yLocMax])
    xlabel('Trial Number')
    ylabel('Estimated vs. Actual Location Difference (°)')
    title(sprintf('Location Difference for %.2f Contrast',cond1cont2TE(1,3)))
    hold on
    plot(repmat(mean(cond1cont2Data(:,2)),1,length(cond1cont2Data(:,2))))
    fprintf('The average location difference for a Contrast of %.2f is %.2f\n',cond1cont2TE(1,3),mean(cond1cont2Data(:,2)));
    legend('Abs value of difference in location','Average Location Difference')
    
     %Contrast 3 Location Difference
    subplot(2,6,4)
    plot(cond1cont3Data(:,2))
    xlim([1 length(cond1cont3Data)])
    ylim([0 yLocMax])
    xlabel('Trial Number')
    ylabel('Estimated vs. Actual Location Difference (°)')
    title(sprintf('Location Difference for %.2f Contrast',cond1cont3TE(1,3)))
    hold on
    plot(repmat(mean(cond1cont3Data(:,2)),1,length(cond1cont3Data(:,2))))
    fprintf('The average location difference for a Contrast of %.2f is %.2f\n',cond1cont3TE(1,3),mean(cond1cont3Data(:,2)));
    legend('Abs value of difference in location','Average Location Difference')
    
     %Contrast 4 Location Difference
    subplot(2,6,5)
    plot(cond1cont4Data(:,2))
    xlim([1 length(cond1cont4Data)])
    ylim([0 yLocMax])
    xlabel('Trial Number')
    ylabel('Estimated vs. Actual Location Difference (°)')
    title(sprintf('Location Difference for %.2f Contrast',cond1cont4TE(1,3)))
    hold on
    plot(repmat(mean(cond1cont4Data(:,2)),1,length(cond1cont4Data(:,2))))
    fprintf('The average location difference for a Contrast of %.2f is %.2f\n',cond1cont4TE(1,3),mean(cond1cont4Data(:,2)));
    legend('Abs value of difference in location','Average Location Difference')
    
     %Contrast 5 Location Difference
    subplot(2,6,6)
    plot(cond1cont5Data(:,2))
    xlim([1 length(cond1cont5Data)])
    ylim([0 yLocMax])
    xlabel('Trial Number')
    ylabel('Estimated vs. Actual Location Difference (°)')
    title(sprintf('Location Difference for %.2f Contrast',cond1cont5TE(1,3)))
    hold on
    plot(repmat(mean(cond1cont5Data(:,2)),1,length(cond1cont5Data(:,2))))
    fprintf('The average location difference for a Contrast of %.2f is %.2f\n',cond1cont5TE(1,3),mean(cond1cont5Data(:,2)));
    legend('Abs value of difference in location','Average Location Difference')

    % Condition 2 Location Graphs %
    fprintf('\nLocation Stats for %s Trials:\n',condition2)
    
    %Contrast 1 Location Difference
    subplot(2,6,8)
    plot(cond2cont1Data(:,2))
    xlim([1 length(cond2cont1Data)])
    ylim([0 yLocMax])
    xlabel('Trial Number')
    ylabel('Estimated vs. Actual Location Difference (°)')
    title(sprintf('Location Difference for %.2f Contrast',cond2cont1TE(1,3)))
    hold on
    plot(repmat(mean(cond2cont1Data(:,2)),1,length(cond2cont1Data(:,2))))
    fprintf('\nThe average location difference for a Contrast of %.2f is %.2f\n',cond2cont1TE(1,3),mean(cond2cont1Data(:,2)));
    legend('Abs value of difference in location','Average Location Difference')

     %Contrast 2 Location Difference
    subplot(2,6,9)
    plot(cond2cont2Data(:,2))
    xlim([1 length(cond2cont2Data)])
    ylim([0 yLocMax])
    xlabel('Trial Number')
    ylabel('Estimated vs. Actual Location Difference (°)')
    title(sprintf('Location Difference for %.2f Contrast',cond2cont2TE(1,3)))
    hold on
    plot(repmat(mean(cond2cont2Data(:,2)),1,length(cond2cont2Data(:,2))))
    fprintf('The average location difference for a Contrast of %.2f is %.2f\n',cond2cont2TE(1,3),mean(cond2cont2Data(:,2)));
    legend('Abs value of difference in location','Average Location Difference')
    
     %Contrast 3 Location Difference
    subplot(2,6,10)
    plot(cond2cont3Data(:,2))
    xlim([1 length(cond2cont3Data)])
    ylim([0 yLocMax])
    xlabel('Trial Number')
    ylabel('Estimated vs. Actual Location Difference (°)')
    title(sprintf('Location Difference for %.2f Contrast',cond2cont3TE(1,3)))
    hold on
    plot(repmat(mean(cond2cont3Data(:,2)),1,length(cond2cont3Data(:,2))))
    fprintf('The average location difference for a Contrast of %.2f is %.2f\n',cond2cont3TE(1,3),mean(cond2cont3Data(:,2)));
    legend('Abs value of difference in location','Average Location Difference')
    
     %Contrast 4 Location Difference
    subplot(2,6,11)
    plot(cond2cont4Data(:,2))
    xlim([1 length(cond2cont4Data)])
    ylim([0 yLocMax])
    xlabel('Trial Number')
    ylabel('Estimated vs. Actual Location Difference (°)')
    title(sprintf('Location Difference for %.2f Contrast',cond2cont4TE(1,3)))
    hold on
    plot(repmat(mean(cond2cont4Data(:,2)),1,length(cond2cont4Data(:,2))))
    fprintf('The average location difference for a Contrast of %.2f is %.2f\n',cond2cont4TE(1,3),mean(cond2cont4Data(:,2)));
    legend('Abs value of difference in location','Average Location Difference')
    
     %Contrast 5 Location Difference
    subplot(2,6,12)
    plot(cond2cont5Data(:,2))
    xlim([1 length(cond2cont5Data)])
    ylim([0 yLocMax])
    xlabel('Trial Number')
    ylabel('Estimated vs. Actual Location Difference (°)')
    title(sprintf('Location Difference for %.2f Contrast',cond2cont5TE(1,3)))
    hold on
    plot(repmat(mean(cond2cont5Data(:,2)),1,length(cond2cont5Data(:,2))))
    fprintf('The average location difference for a Contrast of %.2f is %.2f\n',cond2cont5TE(1,3),mean(cond2cont5Data(:,2)));
    legend('Abs value of difference in location','Average Location Difference')
    
    
    %% STATISTICS WITHIN A RUN %%
    fprintf('\nCorrelation Statistics\n')
    
     rCont = corrcoef(data(:,5),data(:,6)); %correlation coefficient between contrast difference and response time
     fprintf('\nThe correlation coefficient between the difference in contrast and response time is %.2f',rCont(1,2))
     rLoc = corrcoef(data(:,2),data(:,3)); %correlation coefficient between location difference and repsonse time
     fprintf('\nThe correlation coefficient between the difference in location and response time is %.2f',rCont(1,2))
     rContLoc = corrcoef(data(:,2),data(:,5));
     fprintf('\nThe correlation coefficient between contrast difference and lcation difference is %.2f\n',rContLoc(1,2)) 
    end
end
    %% TESTING BETWEEN RUNS %%
    
fprintf('\nTESTING BETWEEN RUNS\n\n');
if strcmp(p.experiment,'test')==1
   fprintf('(Test Experiment: only comparing perception trials to each other)\n\n')
end
if runsCompleted == 4
   run1Stats = theData(1);
        if any(run1Stats.p.trialEvents(:,1) == 1)
            run1Stats.condition = 'Perception';
        elseif any(run1Stats.p.trialEvents(:,1) == 2)
            run1Stats.condition = 'Working Memory';
        else
            error('Problem with setting up conditions from trial events')
        end
   run2Stats = theData(2);
        if any(run2Stats.p.trialEvents(:,1) == 1)
        run2Stats.condition = 'Perception';
        elseif any(run2Stats.p.trialEvents(:,1) == 2)
            run2Stats.condition = 'Working Memory';
        else
            error('Problem with setting up conditions from trial events')
        end
   run3Stats = theData(3);
        if any(run3Stats.p.trialEvents(:,1) == 1)
            run3Stats.condition = 'Perception';
        elseif any(run3Stats.p.trialEvents(:,1) == 2)
            run3Stats.condition = 'Working Memory';
        else
            error('Problem with setting up conditions from trial events')
        end
   run4Stats = theData(4);
        if any(run4Stats.p.trialEvents(:,1) == 1)
            run4Stats.condition = 'Perception';
        elseif any(run4Stats.p.trialEvents(:,1) == 2)
            run4Stats.condition = 'Working Memory';
        else
            error('Problem with setting up conditions from trial events')
        end
    possCombs = nchoosek(4,2);
    run1Data = run1Stats.p.cond1Data;%doing first condition data - only compares the working memory or perception conditions to each other
    run2Data = run2Stats.p.cond1Data;
    run3Data = run3Stats.p.cond1Data;
    run4Data = run4Stats.p.cond1Data;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 6 different combinations of comparison for trials 1-4
        % 3 different potential matchups - p & p, wm & wm, wm & p
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Trial 1 vs. Trial 2
        [results.trials12.h,results.trials12.p] = ttest(run1Data(:,5),run2Data(:,5));
        if results.trials12.h == 1
            fprintf('\nThe estimated contrast differences were statistically significant between %s (trial 1) and %s (trial 2) with a calculated probability of %0.3f\n',run1Stats.condition, run2Stats.condition, results.trials12.p);
        else
            fprintf('The estimated contrast between %s (trial 1) and %s (trial 2) is not statistically significant\n', run1Stats.condition, run2Stats.condition)
        end
             if sum(strcmp(run1Stats.condition,'Perception') || strcmp(run2Stats.condition,'Perception')) == 0
                 results.trials12.condition = 0; % =0 is both working memory conditioned trials
             elseif sum(strcmp(run1Stats.condition,'Perception') || strcmp(run2Stats.condition,'Perception')) == 1
                 results.trials12.condition = 1; % =1 is one working memory trial, and one perception trial
             elseif sum(strcmp(run1Stats.condition,'Perception') || strcmp(run2Stats.condition,'Perception')) == 2 
                 results.trials12.condition = 2; % =2 is both perception conditioned trials
             else
                 error('Error with condition pairings.');
             end

         % Trial 2 vs. Trial 3
        [results.trials23.h,results.trials23.p] = ttest(run2Data(:,5),run3Data(:,5));
        if results.trials23.h == 1
            fprintf('\nThe estimated contrast differences were statistically significant between %s (trial 2) and %s (trial 3) with a calculated probability of %0.3f\n',run2Stats.condition, run3Stats.condition, results.trials23.p);
        else
            fprintf('The estimated contrast between %s (trial 2) and %s (trial 3) is not statistically significant\n', run2Stats.condition, run3Stats.condition)
        end
             if sum(strcmp(run2Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 0
                 results.trials23.condition = 0; % =0 is both working memory conditioned trials
             elseif sum(strcmp(run2Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 1
                 results.trials23.condition = 1; % =1 is one working memory trial, and one perception trial
             elseif sum(strcmp(run2Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 2 
                 results.trials23.condition = 2; % =2 is both perception conditioned trials
             else
                 error('Error with condition pairings.');
             end
 
         % Trial 3 vs. Trial 4
        [results.trials34.h,results.trials34.p] = ttest(run3Data(:,5),run4Data(:,5));
        if results.trials34.h == 1
             fprintf('\nThe estimated contrast differences were statistically significant between %s (trial 3) and %s (trial 4) with a calculated probability of %0.3f\n',run3Stats.condition, run4Stats.condition, results.trials34.p);
        else
             fprintf('The estimated contrast between %s (trial 3) and %s (trial 4) is not statistically significant\n', run3Stats.condition, run4Stats.condition)
        end
             if sum(strcmp(run3Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 0
                 results.trials34.condition = 0; % =0 is both working memory conditioned trials
             elseif sum(strcmp(run3Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 1
                 results.trials34.condition = 1; % =1 is one working memory trial, and one perception trial
             elseif sum(strcmp(run3Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 2 
                 results.trials34.condition = 2; % =2 is both perception conditioned trials
             else
                 error('Error with condition pairings.');
             end

         % Trial 4 vs. Trial 1
        [results.trials41.h,results.trials41.p] = ttest(run4Data(:,5),run1Data(:,5));
        if results.trials41.h == 1
            fprintf('\nThe estimated contrast differences were statistically significant between %s (trial 4) and %s (trial 1) with a calculated probability of %0.3f\n',run4Stats.condition, run1Stats.condition, results.trials41.p);
        else
            fprintf('The estimated contrast between %s (trial 4) and %s (trial 1) is not statistically significant\n', run4Stats.condition, run1Stats.condition)
        end
             if sum(strcmp(run4Stats.condition,'Perception') || strcmp(run1Stats.condition,'Perception')) == 0
                 results.trials41.condition = 0; % =0 is both working memory conditioned trials
             elseif sum(strcmp(run4Stats.condition,'Perception') || strcmp(run1Stats.condition,'Perception')) == 1
                 results.trials41.condition = 1; % =1 is one working memory trial, and one perception trial
             elseif sum(strcmp(run4Stats.condition,'Perception') || strcmp(run1Stats.condition,'Perception')) == 2 
                 results.trials41.condition = 2; %=2 is both perception conditioned trials
             else
                 error('Error with condition pairings.');
             end
 
         % Trial 1 vs. Trial 3
        [results.trials13.h,results.trials13.p] = ttest(run1Data(:,5),run3Data(:,5));
        if results.trials13.h == 1
            fprintf('\nThe estimated contrast differences were statistically significant between %s (trial 1) and %s (trial 3) with a calculated probability of %0.3f\n',run1Stats.condition, run3Stats.condition, results.trials23.p);
        else
            fprintf('The estimated contrast between %s (trial 2) and %s (trial 3) is not statistically significant\n', run1Stats.condition, run3Stats.condition)
        end
             if sum(strcmp(run1Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 0
                 results.trials13.condition = 0; % =0 is both working memory conditioned trials
             elseif sum(strcmp(run1Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 1
                 results.trials13.condition= 1; % =1 is one working memory trial, and one perception trial
             elseif sum(strcmp(run1Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 2 
                 results.trials23.condition = 2; % =2 is both perception conditioned trials
             else
                 error('Error with condition pairings.');
             end

        % Trial 2 vs. Trial 4
        [results.trials24.h,results.trials24.p] = ttest(run2Data(:,5),run4Data(:,5));
        if results.trials24.h == 1
            fprintf('\nThe estimated contrast differences were statistically significant between %s (trial 2) and %s (trial 4) with a calculated probability of %0.3f\n',run2Stats.condition, run4Stats.condition, results.trials24.p);
        else
            fprintf('The estimated contrast between %s (trial 2) and %s (trial 4) is not statistically significant\n', run2Stats.condition, run4Stats.condition)
        end
             if sum(strcmp(run2Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 0
                 results.trials24.condition = 0; % =0 is both working memory conditioned trials
             elseif sum(strcmp(run2Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 1
                 results.trials24.condition = 1; % =1 is one working memory trial, and one perception trial
             elseif sum(strcmp(run2Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 2 
                 results.trials24.condition = 2; % =2 is both perception conditioned trials
             else
                 error('Error with condition pairings.');
             end

%% GRAPHS OF SAVED OUT DATA %%
figure
loglog(p.centerContrast,mean(subject.avgEstContrast(:,:,2),1),'-o') %perception condition
hold on
loglog(p.centerContrast,mean(subject.avgEstContrast(:,:,1),1),'-o') %baseline condition
loglog([0.1 0.8],[0.1 0.8], 'k--')
legend(sprintf('%s Trials',condition1),sprintf('%s Trials',condition2))
set(gcf, 'Name',('Perceived Contrast versus center contrast'));
xticks([0.1 0.8]); yticks([0.1 0.8]);
xticklabels({'10','80'});yticklabels({'10','80'});
xlabel('Center Contrast')
ylabel('Perceived Contrast')
xlim([0.1 0.8])
ylim([0.1 0.8])

figure
loglog(p.centerContrast,mean(subject.avgDiffLoc(:,:,2),1),'-o') %perception second page
hold on
loglog(p.centerContrast,mean(subject.avgDiffLoc(:,:,1),1),'-o') % baseline first page
legend(sprintf('%s Trials',condition1),sprintf('%s Trials',condition2))
set(gcf, 'Name',('Location Difference versus Center Contrast'));
xlabel('Center Contrast')
ylabel('Difference in Location')

figure
loglog(p.centerContrast,mean(subject.avgContDiff(:,:,2),1),'-o') %perception second page
hold on
loglog(p.centerContrast,mean(subject.avgContDiff(:,:,1),1),'-o') % baseline first page
legend(sprintf('%s Trials',condition1),sprintf('%s Trials',condition2))
set(gcf, 'Name',('Contrast Difference versus Center Contrast'));
xlabel('Center Contrast')
ylabel('Difference in Contrast')
end