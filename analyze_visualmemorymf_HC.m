clear all;
close all;

p.experiment = 'test_HC';
p.subject = 'JS';

dataDir = 'data';
cd(dataDir);
if exist(['data_visualmemorymf_' p.experiment '_' p.subject '.mat'],'file') ~= 0
    load(['data_visualmemorymf_' p.experiment '_' p.subject '.mat']);
    p = theData(2).p; data=theData(2).data;
else
    error('Subject file does not exist.')
end


%is there a formula for suppression?

%Histogram of Perception/baseline Results
%Only works for the hardcoded version
[trials,params] = size(p.trialEvents);
pResults = zeros(trials,params);
pData = zeros(trials,params);
blResults = zeros(trials,params); 
blData = zeros(trials,params);
data = cell2mat(struct2cell(data));
data = data';
%%
for i = 1:length(p.trialEvents)
    if p.trialEvents(i,1) == 1
        pResults(i,:) = p.trialEvents(i,:); %Creates new matrix of trials that only had the wm condition
        pData(i,:) = data(i,:);
    elseif p.trialEvents(i,1) == 3
        blResults(i,:) = p.trialEvents(i,:); %Matrix that only tests the center only results
        blData(i,:) = data(i,:);
    else 
        disp('error') % change this process when start testing the perception condition
    end
end
%%
pResults(~any(pResults,2),:) = []; 
pData(~any(pData,2),:) = [];
blResults(~any(blResults,2),:) = [];
blData(~any(blData,2),:) = [];
%%
%Y Limits to compare results of working memory and center only on same
%bounds
%When using multiple contrasts,it is better to set the bounds to [0 1]
minContrast = min([min(pData(:,4)) min(blData(:,4))]);
maxContrast = max([max(pData(:,4)) max(blData(:,4))]);
%lower bound
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
%Perception Estimated Contrast in comparison to Actual Contrast
subplot(2,3,1)
plot(pData(:,4)); %EstimatedContrast is fourth column
hold on
plot(1:length(pResults(:,3)),pResults(:,3));
fit = polyfit(1:length(pData(:,4)),pData(:,4)',1);
plot(polyval(fit,1:length(pData(:,4)))); %Line of best fit for the Estimated Contrast
xlim([1 length(pData(:,4))]);
ylim([yMin yMax]);
title('Perception Contrast Compared to Actual Contrast')
legend('Subject Contrast Estimation','Actual Contrast', 'Average Contrast Estimation')
xlabel('Trial Number')
ylabel('Contrast')

%Center Only Contrast in comparison to actual contrast
subplot(2,3,4)
plot(blData(:,4)); %EstimatedContrast is fourth column
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

%Histogram of working memory contrast estimations
subplot(2,3,2)
[wmBins,wmEdges] = histcounts(pData(:,4),20);
[cOBins,cOEdges] = histcounts(blData(:,4),40);
histmax = max([max(wmBins) max(cOBins)]);
wmHist = histogram(pData(:,4),20);
xlim([0 1]);
ylim([0 histmax]);
title('Histogram of Perception Contrast Estimations');
line([p.centerContrast p.centerContrast],ylim, 'Linewidth',1,'Color','r');
lgd = legend('Bins','Actual Contrast');
lgd.Location = 'northwest';
xlabel('Contrast')
ylabel('Frquency of response')

%Histogram of center only contrast estimations
subplot(2,3,5)
coHist = histogram(blData(:,4),40);
xlim([0 1]);dataDir = 'data';
ylim([0 histmax+1]);
title('Histogram of Baseline Contrast Estimations');
hold on
line([p.centerContrast p.centerContrast],ylim, 'Linewidth',1,'Color','r');
lgd = legend('Bins','Actual Contrast');
lgd.Location = 'northwest';
xlabel('Contrast')
ylabel('Frequency of response')


%Means & Standard Deviations of estimated contrasts
wmStd = std(pData(:,4));
wmMean = mean(pData(:,4));
cOStd =std(blData(:,4));
cOMean = mean(blData(:,4));
fprintf('The average response for perception trials was a contrast of %.4f, with a standard deviation of %.4f.\n',wmMean, wmStd);
fprintf('The average response for baseline trials was a contrast of %.4f, with a standard deviation of %.4f.\n',cOMean, cOStd);

%Percent Error of Contrast Estimations
wmPE = zeros(size(1:length(pData(:,4))));
cOPE = zeros(size(1:length(blData(:,4))));


for i = 1:length(pData(:,4))
    wmPE(i) = (abs((pData(i,4)-pResults(i,3))/pResults(i,3)))*100;
    cOPE(i) = (abs((blData(i,4)-blResults(i,3))/blResults(i,3)))*100;
end

maxPE = max([max(wmPE) max(cOPE)]);
if maxPE + 5 >= 100
    yLimPE = 100;
else
    yLimPE = (maxPE + 5);
end

subplot(2,3,3)
plot(wmPE);
xlim([1 length(pData)])
ylim([0 yLimPE])
title('perception Contrast Percent Error')
hold on
plot(repmat(mean(wmPE),1,length(pData)))
legend('Percent Error of Each Estimation','Average Percent Error')
xlabel('Trial Number')
ylabel('Percent Error')

subplot(2,3,6)
plot(cOPE);
xlim([1 length(blData)])
ylim([0 yLimPE])
title('Baseline Percent Error')
hold on
plot(repmat(mean(cOPE),1,length(blData)))
legend('Percent Error of Each Estimation','Average Percent Error')
xlabel('Trial Number')
ylabel('Percent Error')
