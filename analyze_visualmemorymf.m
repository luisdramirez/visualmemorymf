%% Pull & Organize data

load('data_visualmemorymf_test_JS.mat')
theData = theData(2);
p = theData.p; data = theData.data; t = theData.t;
%5 contrasts, 30 repititions
%column 1 of p.trialEvents: Conditions
    %1 - PERCEPTION
    %2 - WORKING MEMORY
    %3 - BASELINE
% Organize graph into different contrasts for comparisons

data = (cell2mat(struct2cell(data)))';  
sortedTE = sortrows(p.trialEvents,3); %sorts rows based on contrasts
col6 = p.trialEvents(:,6);
data = [data col6];
[dataTrials,dataParams] = size(data);
sortedData = sortrows(data,6);


p.trialEvents = [p.trialEvents (1:length(p.trialEvents))'];
[trials, params] = size(p.trialEvents);
    firstTE = zeros(trials,params);
    secondTE = zeros(trials,params);
    thirdTE = zeros(trials,params);
    fourthTE = zeros(trials,params);
    fifthTE = zeros(trials,params);
for i = 1:length(p.trialEvents)
    if p.trialEvents(i,6) == 1
        firstTE(i,:) = p.trialEvents(i,:);
    elseif p.trialEvents(i,6) == 2
        secondTE(i,:) = p.trialEvents(i,:);
    elseif p.trialEvents(i,6) == 3
        thirdTE(i,:) = p.trialEvents(i,:);
    elseif p.trialEvents(i,6) == 4
        fourthTE(i,:) = p.trialEvents(i,:);
    elseif p.trialEvents(i,6) == 5
        fifthTE(i,:) = p.trialEvents(i,:);
    end
end

firstTE(~any(firstTE,2),:) = [];
secondTE(~any(secondTE,2),:) = [];
thirdTE(~any(thirdTE,2),:) = [];
fourthTE(~any(fourthTE,2),:) = [];
fifthTE(~any(fifthTE,2),:) = [];
orgTE = [firstTE; secondTE; thirdTE; fourthTE; fifthTE];

%% Plot estimated versus actual contrast

figure
subplot(2,3,1)
plot(orgTE(:,3));
ylim([0 1]);
title('Estimated Vs. Actual Contrast')
hold on

% prevent HC - make this a loop that makes the first-fifth data vecs based
% on numContrasts
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

