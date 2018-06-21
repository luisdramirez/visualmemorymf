clear all;
close all;

load('data_visualmemorymf_test_JS.mat') % test_HC or test or regular (?)
theData = theData(2); %index to trial number
p = theData.p; data = theData.data; t = theData.t;
data = cell2mat(struct2cell(data));
data = data';
[trials,params] = size(p.trialEvents);

%HC testing, baseline vs perception
if p.numContrasts == 1 && p.numConditions == 2
    %% Organize Data
% perception trials
pResults = p.trialEvents(p.trialEvents(:,1)==1,:);
pData = data(p.trialEvents(:,1)==1,:);
% baseline trials
blResults = p.trialEvents(p.trialEvents(:,1)==3,:);
blData = data(p.trialEvents(:,1)==3,:);
    %%  Perception Estimated Contrast in comparison to Actual Contrast  %%
    
    % Compare Perception and Baseline on same y-limit
    minContrast = min([min(pData(:,4)) min(blData(:,4))]);
    maxContrast = max([max(pData(:,4)) max(blData(:,4))]);
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
    set(gcf, 'Name', 'Visual Working Memory Versus Perception');
    subplot(2,3,1)
    plot(pData(:,4)); %Estimated Contrast
    hold on
    plot(1:length(pResults(:,3)),pResults(:,3));
    myfittype = fittype('a+b*log(x)','dependent',{'y'}, 'independent',{'x'},'coefficients' ,{'a,b'});
    myfit = fit(1:length(pData(:,4),pData(:,4)',myfittype));
    plot(myfit);
    
%     fit = polyfit(1:length(pData(:,4)),pData(:,4)',1);
%     plot(polyval(fit,1:length(pData(:,4)))); %Line of best fit for the Estimated Contrast
    xlim([1 length(pData(:,4))]);
    ylim([yMin yMax]);
    title('Perception Contrast Compared to Actual Contrast')
    legend('Subject Contrast Estimation','Actual Contrast', 'Average Contrast Estimation')
    xlabel('Trial Number')
    ylabel('Contrast')

    %% Baseline Contrast in comparison to Actual Contrast  %%
    subplot(2,3,4)
    plot(blData(:,4)); %Estimated Contrast
    hold on
    plot(1:length(blResults(:,3)),blResults(:,3));
    fit = polyfit(1:length(blData(:,4)),blData(:,4)',1);
    plot(polyval(fit,1:length(blData(:,4)))); %Line of best fit for the Estimated Contrast
    xlim([1 length(blData(:,4))]);
    ylim([yMin yMax]);
    title('Baseline Contrast Estimations Compared to actual Contrast');
    legend('Subject Contrast Estimation','Actual Contrast', 'Average Contrast Estimation')
    xlabel('Trial Number')
    ylabel('Contrast')

    %Histogram of perception contrast estimations
    subplot(2,3,2)
    [pBins,pEdges] = histcounts(pData(:,4),20);
    [cOBins,cOEdges] = histcounts(blData(:,4),20);
    histmax = max([max(pBins) max(cOBins)]);
    pHist = histogram(pData(:,4),20);
    xlim([0 1]);
    ylim([0 histmax]);
    title('Histogram of Perception Contrast Estimations');
    line([p.centerContrast p.centerContrast],ylim, 'Linewidth',1,'Color','r');
    lgd = legend('Bins','Actual Contrast');
    lgd.Location = 'northwest';
    xlabel('Contrast')
    ylabel('Frquency of response')

    %Histogram of baseline contrast estimations
    subplot(2,3,5)
    blHist = histogram(blData(:,4),20);
    xlim([0 1]);
    ylim([0 histmax+1]);
    title('Histogram of Baseline Contrast Estimations');
    hold on
    line([p.centerContrast p.centerContrast],ylim, 'Linewidth',1,'Color','r');
    lgd = legend('Bins','Actual Contrast');
    lgd.Location = 'northwest';
    xlabel('Contrast')
    ylabel('Frequency of response')


    %Means & Standard Deviations of estimated contrasts
    pStd = std(pData(:,4));
    pMean = mean(pData(:,4));
    blStd =std(blData(:,4));
    blMean = mean(blData(:,4));
    fprintf('The average response for Perception trials was a contrast of %.4f, with a standard deviation of %.4f.\n The absolute difference between the average and actual contrast is %f\n',pMean, pStd,abs(pMean-p.centerContrast));
    fprintf('The average response for Baseline trials was a contrast of %.4f, with a standard deviation of %.4f.\n The absolute difference between the average and actual contrast is %f\n',blMean, blStd, abs(blMean-p.centerContrast));

    %Percent Error of Contrast Estimations
    pPE = zeros(size(1:length(pData(:,4))));
    blPE = zeros(size(1:length(blData(:,4))));


    for i = 1:length(pData(:,4))
        pPE(i) = (abs((pData(i,4)-pResults(i,3))/pResults(i,3)))*100;
        blPE(i) = (abs((blData(i,4)-blResults(i,3))/blResults(i,3)))*100;
    end

    maxPE = max([max(pPE) max(blPE)]);
    if maxPE + 5 >= 100
        yLimPE = 100;
    else
        yLimPE = (maxPE + 5);
    end

    subplot(2,3,3)
    plot(pPE);
    xlim([1 length(pData)])
    ylim([0 yLimPE])
    title('perception Contrast Percent Error')
    hold on
    plot(repmat(mean(pPE),1,length(pData)))
    legend('Percent Error of Each Estimation','Average Percent Error')
    xlabel('Trial Number')
    ylabel('Percent Error')

    subplot(2,3,6)
    plot(blPE);
    xlim([1 length(blData)])
    ylim([0 yLimPE])
    title('Baseline Percent Error')
    hold on
    plot(repmat(mean(blPE),1,length(blData)))
    legend('Percent Error of Each Estimation','Average Percent Error')
    xlabel('Trial Number')
    ylabel('Percent Error')
  %% 5 contrasts
elseif p.numContrasts == 5 && p.numConditions == 2
    %Pull & Organize data
    %5 contrasts, 30 repititions
    %column 1 of p.trialEvents: Conditions
        %1 - PERCEPTION
        %2 - WORKING MEMORY
        %3 - BASELINE
    %Organize graph into different contrasts for comparisons
    sortedTE = sortrows(p.trialEvents,3); %sorts rows based on contrasts
    col6 = p.trialEvents(:,6);
    data = [data col6];
    [dataTrials,dataParams] = size(data);
    sortedData = sortrows(data,6);

    p.trialEvents = [p.trialEvents (1:length(p.trialEvents))']; %adds run number to seventh column in p.trialEvents
    firstTE = p.trialEvents(p.trialEvents(:,6)==1,:);
    secondTE = p.trialEvents(p.trialEvents(:,6)==2,:);
    thirdTE = p.trialEvents(p.trialEvents(:,6)==3,:);
    fourthTE = p.trialEvents(p.trialEvents(:,6)==4,:);
    fifthTE = p.trialEvents(p.trialEvents(:,6)==5,:);
    orgTE = [firstTE; secondTE; thirdTE; fourthTE; fifthTE];
        
    % Plot estimated versus actual contrast

    figure
    subplot(2,3,1)
    plot(orgTE(:,3));
    ylim([0 1]);
    title('Estimated Vs. Actual Contrast')
    hold on

    firstdata = zeros(length(firstTE),dataParams);
    seconddata = zeros(length(secondTE),dataParams);
    thirddata = zeros(length(thirdTE),dataParams);
    fourthdata = zeros(length(fourthTE),dataParams);
    fifthdata = zeros(length(fifthTE),dataParams);
    for i = 1:length(firstTE)
        firstdata(i,:) = data(firstTE(i,7),:);
        seconddata(i,:) = data(secondTE(i,7),:);
        thirddata(i,:) = data(thirdTE(i,7),:);
        fourthdata(i,:) = data(fourthTE(i,7),:);
        fifthdata(i,:) = data(fifthTE(i,7),:);
    end
    orgData = [firstdata; seconddata; thirddata; fourthdata; fifthdata];
    plot(orgData(:,4));
    hold on
    firstMean = (ones([1 length(firstdata)])*mean(firstdata(:,4)));
    secondMean = ones([1 length(seconddata)])*mean(seconddata(:,4));
    thirdMean = ones([1 length(thirddata)])*mean(thirddata(:,4));
    fourthMean = ones([1 length(fourthdata)])*mean(fourthdata(:,4));
    fifthMean = ones([1 length(fifthdata)])*mean(fifthdata(:,4));
    meanVector = [firstMean secondMean thirdMean fourthMean fifthMean];
    plot(1:(length(firstMean)+length(secondMean)+length(thirdMean)+length(fourthMean)+length(fifthMean)),meanVector);
    legend('Estimated Contrast','Actual Contrast','Avg. Estimated Contrast per Level');

    fprintf('\nThe mean of responses for contrast %f was %f. The absolute difference is: %f\n\n',firstTE(1,3),firstMean(1), abs(firstTE(1,3)-firstMean(1)));
    fprintf('The mean of responses for contrast %f was %f. The absolute difference is: %f\n\n',secondTE(1,3),secondMean(1),abs(secondTE(1,3)-secondMean(1)));
    fprintf('The mean of responses for contrast %f was %f. The absolute difference is: %f\n\n',thirdTE(1,3),thirdMean(1),abs(thirdTE(1,3)-thirdMean(1)));
    fprintf('The mean of responses for contrast %f was %f. The absolute difference is: %f\n\n',fourthTE(1,3),fourthMean(1),abs(fourthTE(1,3)-fourthMean(1)));
    fprintf('The mean of responses for contrast %f was %f. The absolute difference is: %f\n\n',fifthTE(1,3),fifthMean(1),abs(fifthTE(1,3)-fifthMean(1)));
    %% Histograms of responses for each contrast level
    subplot(2,3,2)
    hist(firstdata(:,4),20);
    hold on
    xlim([0 1]);
    ylim([0 6]);
    line([firstTE(:,3) firstTE(:,3)],ylim, 'Linewidth',1,'Color','r');
    title(sprintf('Hist of Estimated Contrast for %.2f contrast',firstTE(1,3)));
    legend('Estimated Contrasts', 'Actual Contrast Level');

    subplot(2,3,3)
    hist(seconddata(:,4),20);
    xlim([0 1]);
    ylim([0 6]);
    hold on
    line([secondTE(:,3) secondTE(:,3)],ylim, 'Linewidth',1,'Color','r');
    title(sprintf('Hist of Estimated Contrast for %.2f contrast',secondTE(1,3)));
    legend('Estimated Contrasts', 'Actual Contrast Level');

    subplot(2,3,4)
    hist(thirddata(:,4),20);
    xlim([0 1]);
    ylim([0 6]);
    hold on
    line([thirdTE(:,3) thirdTE(:,3)],ylim, 'Linewidth',1,'Color','r');
    title(sprintf('Hist of Estimated Contrast for %.2f contrast',thirdTE(1,3)));
    legend('Estimated Contrasts', 'Actual Contrast Level');

    subplot(2,3,5)
    hist(fourthdata(:,4),20);
    xlim([0 1]);
    ylim([0 6]);
    hold on
    line([fourthTE(:,3) fourthTE(:,3)],ylim, 'Linewidth',1,'Color','r');
    title(sprintf('Hist of Estimated Contrast for %.2f contrast',fourthTE(1,3)));
    legend('Estimated Contrasts', 'Actual Contrast Level');

    subplot(2,3,6)
    hist(fifthdata(:,4),20);
    xlim([0 1]);
    ylim([0 6]);
    hold on
    line([fifthTE(:,3) fifthTE(:,3)],ylim, 'Linewidth',1,'Color','r');
    title(sprintf('Hist of Estimated Contrast for %.2f contrast',fifthTE(1,3)));
    legend('Estimated Contrasts', 'Actual Contrast Level');
    else
        disp('Number of contrasts or conditions does not correspond to experiment design.')
end