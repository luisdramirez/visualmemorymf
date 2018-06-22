%% Begin %%
% Preliminary data loading and setup %
clear all;
close all;

load('data_visualmemorymf_test_JS.mat') %works for HC, test, and regular trials
theData = theData(2); %index to trial number
p = theData.p; data = theData.data; t = theData.t;
data = cell2mat(struct2cell(data));
data = data';
[trials,params] = size(p.trialEvents);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %Test_HC and 1 contrast runs%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if p.numContrasts == 1 && p.numConditions == 2
    %% Organize Trial Events & Data %%
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

    figure
    set(gcf, 'Name', sprintf('%s Versus %s Visual Analysis',cond1Name,cond2Name));
    subplot(2,3,1)
    plot(cond1Data(:,4)); %Estimated Contrast
    hold on
    plot(1:length(cond1Results(:,3)),cond1Results(:,3));
%     myfittype = fittype('a+b*log(x)','dependent',{'y'}, 'independent',{'x'},'coefficients' ,{'a,b'});
%     myfit = fit(1:length(pData(:,4),pData(:,4)',myfittype));
%     plot(myfit); %%%% LOG FIT %%%%
    fit = polyfit(1:length(cond1Data(:,4)),cond1Data(:,4)',1);
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
    fprintf('\nThe average response for %s trials was a contrast of %.4f, with a standard deviation of %.4f.\n\n The absolute difference between the average and actual contrast is %f\n\n',cond1Name,cond1Mean, cond1Std,abs(cond1Mean-p.centerContrast));
    fprintf('The average response for %s trials was a contrast of %.4f, with a standard deviation of %.4f.\n\n The absolute difference between the average and actual contrast is %f\n\n',cond2Name,cond2Mean, cond2Std, abs(cond2Mean-p.centerContrast));

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
%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Test/5 Contrast Runs%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif p.numContrasts == 5 && p.numConditions == 2
   %% Setup of seperate trial Events and data matrices, & misc. statistical data  %%
        %1 - PERCEPTION
        %2 - WORKING MEMORY
        %3 - BASELINE
    % NOTE: For each subject, index into their data file will give only 1 condition
    sortedTE = sortrows(p.trialEvents,3); %sorts rows based on contrasts
    col6 = p.trialEvents(:,6);
    data = [data col6];
    [dataTrials,dataParams] = size(data);
    sortedData = sortrows(data,6);
    
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
    orgData = [data1; data2; data3; data4; data5];
    plot(orgData(:,4));
    hold on
    mean1 = ones([1 length(data1)])*mean(data1(:,4));
    mean2 = ones([1 length(data2)])*mean(data2(:,4));
    mean3 = ones([1 length(data3)])*mean(data3(:,4));
    mean4 = ones([1 length(data4)])*mean(data4(:,4));
    mean5 = ones([1 length(data5)])*mean(data5(:,4));
    meanVec = [mean1 mean2 mean3 mean4 mean5];
        
    %% Plot Estimated Contrasts vs. Actual Contrasts %%
    figure
    set(gcf, 'Name', sprintf('Contrast Variance over 5 Contrast steps for %s',condition));
    subplot(2,4,1)
    plot(orgData(:,4));
    hold on
    plot(orgTE(:,3));
    ylim([0 1]);
    title('Estimated Vs. Actual Contrast')
    hold on
    plot(1:(length(mean1)+length(mean2)+length(mean3)+length(mean4)+length(mean5)),meanVec);
    legend('Estimated Contrast','Actual Contrast','Avg. Estimated Contrast per Level');
    
    %Display Trendline
    subplot(2,4,2)
    fit = polyfit([1:length(orgData(:,4))],orgData(:,4)',5); %5 polynomial for number of contrasts - trendline
    plot(polyval(fit,1:length(orgData(:,4)))); %Line of best fit for the Estimated Contrast
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
    subplot(2,4,3)
    hist(data1(:,4),10);
    hold on
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    line([TE1(:,3) TE1(:,3)],ylim, 'Linewidth',1,'Color','r');
    title(sprintf('Hist of Estimated Contrast for %.2f Contrast',TE1(1,3)));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    
    % Second Contrast Hist
    subplot(2,4,4)
    hist(data2(:,4),10);
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    hold on
    line([TE2(:,3) TE2(:,3)],ylim, 'Linewidth',1,'Color','r');
    title(sprintf('Hist of Estimated Contrast for %.2f contrast',TE2(1,3)));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    
    % Third Contrast Hist
    subplot(2,4,5)
    hist(data3(:,4),10);
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    hold on
    line([TE3(:,3) TE3(:,3)],ylim, 'Linewidth',1,'Color','r');
    title(sprintf('Hist of Estimated Contrast for %.2f contrast',TE3(1,3)));
    legend('Estimated Contrasts', 'Actual Contrast Level'); 
    
    % Fourth Contrast Hist
    subplot(2,4,6)
    hist(data4(:,4),10);
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    hold on
    line([TE4(:,3) TE4(:,3)],ylim, 'Linewidth',1,'Color','r');
    title(sprintf('Hist of Estimated Contrast for %.2f contrast',TE4(1,3)));
    legend('Estimated Contrasts', 'Actual Contrast Level'); 
    
    % Fifth Contrast Hist
    subplot(2,4,7)
    hist(data5(:,4),10);
    xlim([0 1]);
    ylim([0 histmax+0.5]);
    hold on
    line([TE5(:,3) TE5(:,3)],ylim, 'Linewidth',1,'Color','r');
    title(sprintf('Hist of Estimated Contrast for %.2f contrast',TE5(1,3)));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    else
        disp('Number of contrasts or conditions does not correspond to experiment design.')
end