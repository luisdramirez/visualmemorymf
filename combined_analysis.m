%% Analysis designed to evaluate and show visual representations of data from one single run %%
% Preliminary data loading and setup %
clear all;
close all;
expDir=pwd;
dataDir = 'data_master';
p.experiment = 'test';
p.subject = 'JS';
cd(dataDir)

if exist(['data_visualmemorymf_' p.experiment '_' p.subject '.mat'],'file') ~= 0
    load(['data_visualmemorymf_' p.experiment '_' p.subject '.mat']); %works for HC, test, and regular trials
    theDatacurr = theData(1); %index to trial number
    p = theDatacurr.p; data = theDatacurr.data; t = theDatacurr.t;
    data = cell2mat(struct2cell(data));
    data = data';
    [trials,params] = size(p.trialEvents);
else
    error('data file does not exist')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %Test_HC and 1 Contrast Runs%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for i = 1:length(theData) - this displays all at same time on same graphs
if p.numContrasts == 1 && p.numConditions == 2
  %%  Organize Trial Events & Data %%
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
    [bins1,edges1] = histcounts(cond1Data(:,4),20);
    [bins2,edges2] = histcounts(cond2Data(:,4),20);
    histmax = max([max(bins1) max(bins2)]);
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
    fprintf('The average response for %s trials was a contrast of %.4f, with a standard deviation of %.4f.\n The absolute difference between the average and actual contrast is %f\n',cond1Name,cond1Mean, cond1Std,abs(cond1Mean-p.centerContrast));
    fprintf('The average response for %s trials was a contrast of %.4f, with a standard deviation of %.4f.\n The absolute difference between the average and actual contrast is %f\n',cond2Name,cond2Mean, cond2Std, abs(cond2Mean-p.centerContrast));

    %Percent Error of Contrast Estimations
    cond1PE = zeros(size(1:length(cond1Data(:,4))));
    cond2PE = zeros(size(1:length(cond2Data(:,4))));
    for i = 1:length(cond1Data(:,4))
        cond1PE(i) = (abs((cond1Data(i,4)-cond1Results(i,3))/cond1Results(i,3)))*100;
        cond2PE(i) = (abs((cond2Data(i,4)-cond2Results(i,3))/cond2Results(i,3)))*100;
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
        fprintf('The reported contrasts for %s and %s were not statistically significant\n',cond1Name,cond2Name);
    elseif h == 1
        fprintf('The reported contrasts for %s and %s were statistically significant from each other, with a p value of %.3f\n',cond1Name,cond2Name,p);
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
        fprintf('%s had a higher estimated contrast than %s\n',cond1Name,cond2Name);
        fprintf('with a percent difference of %.3f%\n\n',(cond1Mean-cond2Mean)*100);
    elseif cond1Mean < cond2Mean
        fprintf('%s had a higher estimated contrast than %s\n',cond2Name,cond1Name);
        fprintf('with a percent difference of %.3f%\n\n',abs(cond1Mean-cond2Mean)*100);
    end
%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Test/5 Contrast Runs%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif p.numContrasts == 5 && p.numConditions == 2
  %for i = 1:length(theData)
%       theDatacurr = theData(i); %index to trial number
    p = theDatacurr.p; data = theDatacurr.data; t = theDatacurr.t;
    data = cell2mat(struct2cell(data));
    %% Setup of seperate trial Events and data matrices, & misc. statistical data  %%
        %1 - PERCEPTION
        %2 - WORKING MEMORY
        %3 - BASELINE
    % NOTE: For each subject, index into their data file will give only 1 condition
    
    sortedTE = sortrows(p.trialEvents,[1,3]); %sorts based off of conditions, then based off of contrasts, (based off of trial num?)
    
    col6 = p.trialEvents(:,6);
    data = [data' col6];
    [dataTrials,dataParams] = size(data);
    sortedData = sortrows(data,7);
    
    if any(p.trialEvents(:,1) == 1)
        condition = 'Perception';
    elseif any(p.trialEvents(:,1) == 2)
        condition = 'Working Memory';
    elseif any(p.trialEvents(:,1) ==3)
        condition = 'Baseline';
    end
    
    %Seperate trial events based on contrast level
    p.trialEvents = [p.trialEvents (1:length(p.trialEvents))']; %adds run number to seventh column in p.trialEvents for sorting purposes
    TE1 = p.trialEvents(p.trialEvents(:,6)==1,:);
    TE2 = p.trialEvents(p.trialEvents(:,6)==2,:);
    TE3 = p.trialEvents(p.trialEvents(:,6)==3,:);
    TE4 = p.trialEvents(p.trialEvents(:,6)==4,:);
    TE5 = p.trialEvents(p.trialEvents(:,6)==5,:);
    orgTE = [TE1; TE2; TE3; TE4; TE5];
    
    data1 = zeros(length(TE1),dataParams);
    data2 = zeros(length(TE2),dataParams);
    data3 = zeros(length(TE3),dataParams);
    data4 = zeros(length(TE4),dataParams);
    data5 = zeros(length(TE5),dataParams);
    for i = 1:length(TE1)
        data1(i,:) = data(TE1(i,7),:);
        data2(i,:) = data(TE2(i,7),:);
        data3(i,:) = data(TE3(i,7),:);
        data4(i,:) = data(TE4(i,7),:);
        data5(i,:) = data(TE5(i,7),:);
    end
    %Means
    mean1 = ones([1 length(data1)])*mean(data1(:,4));
    mean2 = ones([1 length(data2)])*mean(data2(:,4));
    mean3 = ones([1 length(data3)])*mean(data3(:,4));
    mean4 = ones([1 length(data4)])*mean(data4(:,4));
    mean5 = ones([1 length(data5)])*mean(data5(:,4));
    meanVec = [mean1 mean2 mean3 mean4 mean5];
        
    %% Plot Estimated Contrasts vs. Actual Contrasts %%
    figure(1)
    set(gcf, 'Name', sprintf('Contrast Statistics over 5 Contrasts for %s',condition));
    subplot(2,6,1)
    plot(sortedData(:,4));
    hold on
    plot(orgTE(:,3));
    ylim([0 1]);
    title('Estimated Vs. Actual Contrast')
    hold on
    plot(1:length(meanVec),meanVec);
    legend('Estimated Contrast','Actual Contrast','Avg. Estimated Contrast per Level');
    
    %Display Trendline
    subplot(2,6,7)
    fit = polyfit([1:length(sortedData(:,4))],sortedData(:,4)',5); %5 polynomial for number of contrasts - visual trendline
    plot(polyval(fit,1:length(sortedData(:,4)))); %Line of best fit for the Estimated Contrast
    ylim([0 1])
    hold on
    plot(orgTE(:,3));
    hold on
    title('Trendline of results over Contrast Step')
    legend('Trendline to Fifth Degree','Actual Contrast')
%     contrastvec = unique(orgTE(:,3))';
%     middlevec = [(p.numTrials/p.numContrasts)/2 (p.numTrials/p.numContrasts)/2+(p.numTrials/p.numContrasts) (p.numTrials/p.numContrasts)/2+((p.numTrials/p.numContrasts)*2) (p.numTrials/p.numContrasts)/2+((p.numTrials/p.numContrasts)*3) (p.numTrials/p.numContrasts)/2+((p.numTrials/p.numContrasts)*4)];
%     plot(middlevec, contrastvec,'-o');
    
    
    % Print statistical data
    fprintf('\nThe mean of responses for contrast %f was %f. The absolute difference is: %f\n\n',TE1(1,3),mean1(1), abs(TE1(1,3)-mean1(1)));
    fprintf('The mean of responses for contrast %f was %f. The absolute difference is: %f\n\n',TE2(1,3),mean2(1),abs(TE2(1,3)-mean2(1)));
    fprintf('The mean of responses for contrast %f was %f. The absolute difference is: %f\n\n',TE3(1,3),mean3(1),abs(TE3(1,3)-mean3(1)));
    fprintf('The mean of responses for contrast %f was %f. The absolute difference is: %f\n\n',TE4(1,3),mean4(1),abs(TE4(1,3)-mean4(1)));
    fprintf('The mean of responses for contrast %f was %f. The absolute difference is: %f\n\n',TE5(1,3),mean5(1),abs(TE5(1,3)-mean5(1)));
    
    %% Histograms of Responses for Each Contrast Step %%
    
    % Histogram axis limits
    [bins1,edges1] = histcounts(data1(:,4),10);
    [bins2,edges2] = histcounts(data2(:,4),10);
    [bins3,edges3] = histcounts(data3(:,4),10);
    [bins4,edges4] = histcounts(data4(:,4),10);
    [bins5,edges5] = histcounts(data5(:,4),10);
    histmax = max([max(bins1) max(bins2) max(bins3) max(bins4) max(bins5)]);
    
    % First Contrast Hist
    subplot(2,6,2)
    hist(data1(:,4),10);
    hold on
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    line([TE1(:,3) TE1(:,3)],ylim, 'Linewidth',1,'Color','r');
    title(sprintf('Hist of Estimated Contrast for %.2f Contrast',TE1(1,3)));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    
    % Second Contrast Hist
    subplot(2,6,3)
    hist(data2(:,4),10);
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    hold on
    line([TE2(:,3) TE2(:,3)],ylim, 'Linewidth',1,'Color','r');
    title(sprintf('Hist of Estimated Contrast for %.2f contrast',TE2(1,3)));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    
    % Third Contrast Hist
    subplot(2,6,4)
    hist(data3(:,4),10);
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    hold on
    line([TE3(:,3) TE3(:,3)],ylim, 'Linewidth',1,'Color','r');
    title(sprintf('Hist of Estimated Contrast for %.2f contrast',TE3(1,3)));
    legend('Estimated Contrasts', 'Actual Contrast Level'); 
    
    % Fourth Contrast Hist
    subplot(2,6,5)
    hist(data4(:,4),10);
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    hold on
    line([TE4(:,3) TE4(:,3)],ylim, 'Linewidth',1,'Color','r');
    title(sprintf('Hist of Estimated Contrast for %.2f contrast',TE4(1,3)));
    legend('Estimated Contrasts', 'Actual Contrast Level'); 
    
    % Fifth Contrast Hist
    subplot(2,6,6)
    hist(data5(:,4),10);
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    hold on
    line([TE5(:,3) TE5(:,3)],ylim, 'Linewidth',1,'Color','r');
    title(sprintf('Hist of Estimated Contrast for %.2f contrast',TE5(1,3)));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    
    %% Percent Error of Contrast Responses %%
    
    %Percent Error Calculations (includes location PE)
    contrast1PE = repmat(zeros(size(1:length(data1(:,4)))),2,1);
    contrast2PE = repmat(zeros(size(1:length(data2(:,4)))),2,1);
    contrast3PE = repmat(zeros(size(1:length(data3(:,4)))),2,1);
    contrast4PE = repmat(zeros(size(1:length(data4(:,4)))),2,1);
    contrast5PE = repmat(zeros(size(1:length(data5(:,4)))),2,1);
    for i = 1:(length(sortedTE)/p.numContrasts)
        contrast1PE(1,i) = (abs((data1(i,4)-TE1(i,3))/TE1(i,3)))*100;
        contrast1PE(2,i) = (abs((data1(i,2)/TE1(i,2))))*100;
        contrast2PE(1,i) = (abs((data2(i,4)-TE2(i,3))/TE2(i,3)))*100;
        contrast2PE(2,i) = (abs((data2(i,2)/TE2(i,2))))*100;
        contrast3PE(1,i) = (abs((data3(i,4)-TE3(i,3))/TE3(i,3)))*100;
        contrast3PE(2,i) = (abs((data3(i,2)/TE3(i,2))))*100;
        contrast4PE(1,i) = (abs((data4(i,4)-TE4(i,3))/TE4(i,3)))*100;
        contrast4PE(2,i) = (abs((data4(i,2)/TE4(i,2))))*100;
        contrast5PE(1,i) = (abs((data5(i,4)-TE5(i,3))/TE5(i,3)))*100;
        contrast5PE(2,i) = (abs((data5(i,2)/TE5(i,2))))*100;
    end
   yLimPE = max([max(contrast1PE) max(contrast2PE) max(contrast3PE) max(contrast4PE) max(contrast5PE) ]) + 5;
    
    %% Plotting Percent Errors for 5 Contrasts %%
    
    % Contrast 1 Percent Error
    subplot(2,6,8)
    plot(contrast1PE(1,:));
    xlim([1 length(contrast1PE)])
    ylim([0 yLimPE])
    title(sprintf('%.2f Contrast Percent Error',TE1(1,3)))
    hold on
    plot(repmat(mean(contrast1PE(1,:)),1,length(contrast1PE)))
    legend('Percent Error of Each Estimation','Average Percent Error')
    xlabel('Trial Number')
    ylabel('Percent Error')
    
    % Contrast 2 Percent Error
    subplot(2,6,9)
    plot(contrast2PE(1,:));
    xlim([1 length(contrast2PE)])
    ylim([0 yLimPE])
    title(sprintf('%.2f Contrast Percent Error',TE2(1,3)))
    hold on
    plot(repmat(mean(contrast2PE(1,:)),1,length(contrast2PE)))
    legend('Percent Error of Each Estimation','Average Percent Error')
    xlabel('Trial Number')
    ylabel('Percent Error')
    
    % Contrast 3 Percent Error
    subplot(2,6,10)
    plot(contrast3PE(1,:));
    xlim([1 length(contrast3PE)])
    ylim([0 yLimPE])
    title(sprintf('%.2f Contrast Percent Error',TE3(1,3)))
    hold on
    plot(repmat(mean(contrast3PE(1,:)),1,length(contrast3PE)))
    legend('Percent Error of Each Estimation','Average Percent Error')
    xlabel('Trial Number')
    ylabel('Percent Error')
    
    % Contrast 4 Percent Error
    subplot(2,6,11)
    plot(contrast4PE(1,:));
    xlim([1 length(contrast4PE)])
    ylim([0 yLimPE])
    title(sprintf('%.2f Contrast Percent Error',TE4(1,3)))
    hold on
    plot(repmat(mean(contrast4PE(1,:)),1,length(contrast4PE)))
    legend('Percent Error of Each Estimation','Average Percent Error')
    xlabel('Trial Number')
    ylabel('Percent Error')
    
    % Contrast 5 Percent Error
    subplot(2,6,12)
    plot(contrast5PE(1,:));
    xlim([1 length(contrast5PE)])
    ylim([0 yLimPE])
    title(sprintf('%.2f Contrast Percent Error',TE5(1,3)))
    hold on
    plot(repmat(mean(contrast5PE(1,:)),1,length(contrast5PE)))
    legend('Percent Error of Each Estimation','Average Percent Error')
    xlabel('Trial Number')
    ylabel('Percent Error')
    
    
    %% Location Analysis %%

    % Location Difference in relation to differnt contrasts - any discrepancies
    % between average difference by contrast level?
    
    %Statistical Information
%     locDiff = sortedData(:,2);
%     [h,p,ci,stats] = ttest(sortedData(:,2),mean(sortedData(:,2))); % h

    %Y Limit for accurate comparisons
    yLocMax = max([max(data1(:,2)) max(data2(:,2)) max(data3(:,2)) max(data4(:,2)) max(data5(:,2))])+5;
    figure(2)
    set(gcf, 'Name', sprintf('Location Statistics over 5 Contrasts for %s',condition));

    % General Location Differences and Trendline
    subplot(2,6,1)
    plot(sortedData(:,2))
    hold on
    ylim([0 yLocMax])
    xlabel('Lowest to Highest Contrast Trials')
    ylabel('Estimated vs. Actual Location Difference (°)')
    title('Location Difference for all Contrasts')
    fit = polyfit(1:length(sortedData(:,2)),sortedData(:,2)',1);
    plot(polyval(fit,1:length(sortedData(:,2)))); 
    legend('Abs value of difference in location','Average Location Difference Trendline')

    %Contrast 1 Location Difference
    subplot(2,6,2)
    plot(data1(:,2))
    xlim([1 length(data1)])
    ylim([0 yLocMax])
    xlabel('Trial Number')
    ylabel('Estimated vs. Actual Location Difference (°)')
    title(sprintf('Location Difference for %.2f Contrast',TE1(1,3)))
    hold on
    plot(repmat(mean(data1(:,2)),1,length(data1(:,2))))
    fprintf('\nThe average location difference for a Contrast of %.2f is %.2f\n',TE1(1,3),mean(data1(:,2)));
    legend('Abs value of difference in location','Average Location Difference')

    %Contrast 2 Location Difference
    subplot(2,6,3)
    plot(data2(:,2))
    xlim([1 length(data2)])
    ylim([0 yLocMax])
    xlabel('Trial Number')
    ylabel('Estimated vs. Actual Location Difference (°)')
    title(sprintf('Location Difference for %.2f Contrast',TE2(1,3)))
    hold on
    plot(repmat(mean(data2(:,2)),1,length(data2(:,2))))
    fprintf('\nThe average location difference for a Contrast of %.2f is %.2f\n',TE2(1,3),mean(data2(:,2)));
    legend('Abs value of difference in location','Average Location Difference')

    %Contrast 3 Location Difference
    subplot(2,6,4)
    plot(data3(:,2))
    xlim([1 length(data3)])
    ylim([0 yLocMax])
    xlabel('Trial Number')
    ylabel('Estimated vs. Actual Location Difference (°)')
    title(sprintf('Location Difference for %.2f Contrast',TE3(1,3)))
    hold on
    plot(repmat(mean(data3(:,2)),1,length(data3(:,2))))
    fprintf('\nThe average location difference for a Contrast of %.2f is %.2f\n',TE3(1,3),mean(data3(:,2)));
    legend('Abs value of difference in location','Average Location Difference')

    %Contrast 4 Location Difference
    subplot(2,6,5)
    plot(data4(:,2))
    xlim([1 length(data4)])
    ylim([0 yLocMax])
    xlabel('Trial Number')
    ylabel('Estimated vs. Actual Location Difference (°)')
    title(sprintf('Location Difference for %.2f Contrast',TE4(1,3)))
    hold on
    plot(repmat(mean(data4(:,2)),1,length(data4(:,2))))
    fprintf('\nThe average location difference for a Contrast of %.2f is %.2f\n',TE4(1,3),mean(data4(:,2)));
    legend('Abs value of difference in location','Average Location Difference')

    %Contrast 5 Location Difference
    subplot(2,6,6)
    plot(data5(:,2))
    xlim([1 length(data5)])
    ylim([0 yLocMax])
    xlabel('Trial Number')
    ylabel('Estimated vs. Actual Location Difference (°)')
    title(sprintf('Location Difference for %.2f Contrast',TE5(1,3)))
    hold on
    plot(repmat(mean(data5(:,2)),1,length(data5(:,2))))
    fprintf('\nThe average location difference for a Contrast of %.2f is %.2f\n',TE5(1,3),mean(data5(:,2)));
    legend('Abs value of difference in location','Average Location Difference')


%% Percent Error in Location Estimations %%

    yLimPE = max([max(contrast1PE(2,:)) max(contrast2PE(2,:)) max(contrast3PE(2,:)) max(contrast4PE(2,:)) max(contrast5PE(2,:))])+5;

    % Location Percent Error, Contrast 1 
    subplot(2,6,8)
    plot(contrast1PE(2,:));
    xlim([1 length(contrast1PE)])
    ylim([0 yLimPE])
    title(sprintf('Location Difference Percent Error for %.2f',TE1(1,3)))
    hold on
    plot(repmat(mean(contrast1PE(2,:)),1,length(contrast1PE)))
    legend('Percent Error of Each Estimation','Average Percent Error')
    xlabel('Trial Number')
    ylabel('Percent Error')
    legend('Percent Error','Average Percent Error')

    % Location Percent Error, Contrast 2 
    subplot(2,6,9)
    plot(contrast2PE(2,:));
    xlim([1 length(contrast2PE)])
    ylim([0 yLimPE])
    title(sprintf('Location Difference Percent Error for %.2f',TE2(1,3)))
    hold on
    plot(repmat(mean(contrast2PE(2,:)),1,length(contrast2PE)))
    legend('Percent Error of Each Estimation','Average Percent Error')
    xlabel('Trial Number')
    ylabel('Percent Error')
    legend('Percent Error','Average Percent Error')
    
    % Location Percent Error, Contrast 3 
    subplot(2,6,10)
    plot(contrast3PE(2,:));
    xlim([1 length(contrast3PE)])
    ylim([0 yLimPE])
    title(sprintf('Location Difference Percent Error for %.2f',TE3(1,3)))
    hold on
    plot(repmat(mean(contrast3PE(2,:)),1,length(contrast3PE)))
    legend('Percent Error of Each Estimation','Average Percent Error')
    xlabel('Trial Number')
    ylabel('Percent Error')
    legend('Percent Error','Average Percent Error')
    
    % Location Percent Error, Contrast 4 
    subplot(2,6,11)
    plot(contrast4PE(2,:));
    xlim([1 length(contrast4PE)])
    ylim([0 yLimPE])
    title(sprintf('Location Difference Percent Error for %.2f',TE4(1,3)))
    hold on
    plot(repmat(mean(contrast4PE(2,:)),1,length(contrast4PE)))
    legend('Percent Error of Each Estimation','Average Percent Error')
    xlabel('Trial Number')
    ylabel('Percent Error')
    legend('Percent Error','Average Percent Error')
    
    % Location Percent Error, Contrast 5 
    subplot(2,6,12)
    plot(contrast5PE(2,:));
    xlim([1 length(contrast5PE)])
    ylim([0 yLimPE])
    title(sprintf('Location Difference Percent Error for %.2f',TE5(1,3)))
    hold on
    plot(repmat(mean(contrast5PE(2,:)),1,length(contrast5PE)))
    legend('Percent Error of Each Estimation','Average Percent Error')
    xlabel('Trial Number')
    ylabel('Percent Error')
    legend('Percent Error','Average Percent Error')
    
    %% Timing Parameters %%
    contrastTime = data(:,6) - data(:,3);
    locationTime = data(:,3);
    timingData = [contrastTime locationTime];
    lengthTD = length(timingData);
%     for i = 1:lengthTD
%         if any(timingData(i,:) > 10)
%             timingData(i,:) = [];
%         end
%     end
%     for i = 1:length(timingData)
%         m = timingData(i,1);
%         n = timingData(i,2);
%         if m >10 || n>10
%             timingData(i,:) = [];
%         end
%     end
%   trying to figure out how to take outliers out (over 10ish seconds?)   
    figure(3)
    set(gcf, 'Name', sprintf('Timing Statistics for %s',condition));
    bar(1:2,mean(timingData));
    hold on
    scatter(ones(length(timingData),1),timingData(:,1))
    scatter(2*ones(length(timingData),1),timingData(:,2))
    set(gca,'XtickLabel',{'Contrast RT' , 'LocationRT'})
    
    %% Organize statistical information to compare against other trials %%
    stats.data = data;
    stats.data1 = data1;
    stats.data2 = data2;
    stats.data3 = data3;
    stats.data4 = data4;
    stats.data5 = data5;
    stats.sortedData = sortedData;
    stats.trialEvents = p.trialEvents;
    stats.orgTE = orgTE;
    stats.mean1 = mean1;
    stats.mean2 = mean2;
    stats.mean3 = mean3;
    stats.mean4 = mean4;
    stats.mean5 = mean5;
    stats.condition = condition;
    theData(p.runNumber).stats = stats;
    save(['data_visualmemorymf_' p.experiment '_' p.subject '.mat'], 'theData')
%   end
    else
        disp('Number of contrasts or conditions does not correspond to experiment design.')
end
cd(expDir)