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

        figure(1)
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
    hold off
    title(sprintf('Hist of Estimated Contrast for %.2f Contrast',cond1cont1TE(1,3)));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    
    subplot (2,6,3)
    hist(cond1cont2Data(:,4),10);
    hold on
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    line([cond1cont2TE(:,3) cond1cont2TE(:,3)],ylim, 'Linewidth',1,'Color','r');
    hold off
    title(sprintf('Hist of Estimated Contrast for %.2f Contrast',cond1cont2TE(1,3)));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    
    subplot (2,6,4)
    hist(cond1cont3Data(:,4),10);
    hold on
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    line([cond1cont3TE(:,3) cond1cont3TE(:,3)],ylim, 'Linewidth',1,'Color','r');
    hold off
    title(sprintf('Hist of Estimated Contrast for %.2f Contrast',cond1cont3TE(1,3)));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    subplot (2,6,5)
    hist(cond1cont4Data(:,4),10);
    hold on
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    line([cond1cont4TE(:,3) cond1cont1TE(:,3)],ylim, 'Linewidth',1,'Color','r');
    hold off
    title(sprintf('Hist of Estimated Contrast for %.2f Contrast',cond1cont4TE(1,3)));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    subplot (2,6,6)
    hist(cond1cont5Data(:,4),10);
    hold on
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    line([cond1cont5TE(:,3) cond1cont5TE(:,3)],ylim, 'Linewidth',1,'Color','r');
    hold off
    title(sprintf('Hist of Estimated Contrast for %.2f Contrast',cond1cont5TE(1,3)));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    %%% Condition 2 Contrast Graphs %%%
    
    end
end