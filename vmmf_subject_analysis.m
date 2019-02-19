%%% vmmf_subject_analysis.m
%% Prepare data and working directories
clear all; close all; clc;

subjPlots = 0;
groupPlots = 1;
superPlots = 0;

expDir=pwd;
dataDir='data_master';
subjects = 1:9;
% Load in data
cd(dataDir)
load('visualmemory_subjectsRan')
subjectProfile=struct('SubjectName',[] ,'Order', [], 'Condition', [], 'TheData',[],'OrganizedData',[],'LocationError', [],'ContrastEstimate',[]);

subjectProfile.SubjectName = nan(1,length(subjects)); % initialize subject names
subjectProfile.SubjectName = subjects;
subjectProfile.Order = nan(1,length(subjects)); % initialize order
subjectProfile.TheData = cell(1,length(subjects)); % intialize theData holders
subjectProfile.OrganizedData = cell(1,length(subjects)); % initialized organize data
subjectProfile.LocationError = nan(1,length(subjects)); % initialize location error (Average location error for each subject)
subjectProfile.ContrastEstimate = nan(length(subjects),5,3); % initialize contrast estimate (Average contrast estimate for each condition)

% Load in theData
totalNumTrials = 0;
for currSubj=1:numel(subjectProfile.SubjectName)
    if exist(['data_vmmf_' num2str(subjectProfile.SubjectName(currSubj)) '.mat'],'file') ~= 0
        load(['data_vmmf_' num2str(subjectProfile.SubjectName(currSubj))]);
    end
    subjectProfile.Order(currSubj) = strcmp(visualmemory_subjectsRan(2,currSubj),'a'); % report order
    subjectProfile.Condition(:,currSubj) = theData(1).p.trialSchedule; %subject condition schedule
    subjectProfile.TheData{currSubj} = theData; %subject data
    totalNumTrials = totalNumTrials + (numel(subjectProfile.TheData{currSubj})*subjectProfile.TheData{currSubj}(1).p.numTrials);
end
cd(expDir)

%% Data Organization
centerContrast = unique(subjectProfile.TheData{1}(1).p.trialEvents(:,3));
allLocationError = [];
Conditions = [];
Orders = [];
Contrasts = [];
probeOffset = [];
contrastError = [];
contrastEstimates = [];
index_probeSelective = [];

% Organize Contrast and Location  Estimates
for currSubj = 1:numel(subjectProfile.SubjectName)
    ContrastData = struct('Perception',{[] [] [] [] []},'WorkingMemory',{[] [] [] [] []},'Baseline',{[] [] [] [] []}); % 5 slots for each contrast
    ContrastData_PS = struct('Perception',{[] [] [] [] []},'WorkingMemory',{[] [] [] [] []},'Baseline',{[] [] [] [] []}); % 5 slots for each contrast
    dataFields = fieldnames(ContrastData); %Store contrastData field names
    
    % Initialize location error and condition tracking matrices
    locationErrorMat = nan(length(subjectProfile.TheData{currSubj}(1).p.trialEvents),numel(subjectProfile.TheData{currSubj}));
    conditionMat = nan(length(subjectProfile.TheData{currSubj}(1).p.trialEvents),numel(subjectProfile.TheData{currSubj}));
    
    for currRun = 1:numel(subjectProfile.TheData{currSubj})
        % Grab the condition of each trial in current run
        conditionMat(:,currRun) = subjectProfile.TheData{currSubj}(currRun).p.trialEvents(:,1);
        
        % Keep track of current order
        if currRun <= 4
            currOrder = subjectProfile.Order(currSubj);
        else
            currOrder = ~subjectProfile.Order(currSubj);
        end
        
        % Go through each field in allData structure (Perception, Working
        % Memory, Baseline) and get the mean estimate for a contrast level
        for currField = 1:numel(dataFields)
            for currContrast = 1:numel(centerContrast)
                relevantTrials = subjectProfile.TheData{currSubj}(currRun).p.trialEvents(:,1) == currField & subjectProfile.TheData{currSubj}(currRun).p.trialEvents(:,3) == centerContrast(currContrast);
                ContrastData(currContrast).(dataFields{currField})(currRun,1) = mean(subjectProfile.TheData{currSubj}(currRun).data.EstimatedContrast(relevantTrials)); %Contrast Estimate
                ContrastData(currContrast).(dataFields{currField})(currRun,2) = currOrder; % order of current run   
                
                %Take probe selective contrast data - probes not equal to
                %0.1 or 0.75
                relevantTrials_probeSelective = subjectProfile.TheData{currSubj}(currRun).p.trialEvents(:,1) == currField ...
                    & subjectProfile.TheData{currSubj}(currRun).p.trialEvents(:,3) == centerContrast(currContrast) & subjectProfile.TheData{currSubj}(currRun).p.trialEvents(:,5) ~= 0.75 ...
                    & subjectProfile.TheData{currSubj}(currRun).p.trialEvents(:,5) ~= 0.1; %Takes same data, minus trials that include 0.1 or 0.75 contrast
                ContrastData_PS(currContrast).(dataFields{currField})(currRun,1) = mean(subjectProfile.TheData{currSubj}(currRun).data.EstimatedContrast(relevantTrials_probeSelective)); %Contrast Estimate
                ContrastData_PS(currContrast).(dataFields{currField})(currRun,2) = currOrder; % order of current run   
            end
        end
        
        % Reformat location data and determine probe offset from target
        locationErrorMat(:,currRun) = subjectProfile.TheData{currSubj}(currRun).p.trialEvents(:,2) - subjectProfile.TheData{currSubj}(currRun).data.EstimatedLocation(:);
        
        % Format location space to be -180 to 180
        locationErrorMat(locationErrorMat > 180) = locationErrorMat(locationErrorMat > 180) - 360;
        locationErrorMat(locationErrorMat < -180) = locationErrorMat(locationErrorMat < -180) + 360;
        probeOffset = [probeOffset; subjectProfile.TheData{currSubj}(currRun).p.trialEvents(:,2)-subjectProfile.TheData{currSubj}(currRun).p.trialEvents(:,4)];
        probeOffset(probeOffset > 180) = probeOffset(probeOffset > 180) - 360;
        probeOffset(probeOffset < -180) = probeOffset(probeOffset < -180) + 360;
        
        % Store contrast error
        contrastError = [contrastError; subjectProfile.TheData{currSubj}(currRun).data.DifferenceContrast(:)];
        Orders = [Orders; repmat(currOrder,numel(relevantTrials),1)];
        Contrasts = [Contrasts; subjectProfile.TheData{currSubj}(currRun).p.trialEvents(:,3)];
        contrastEstimates = [contrastEstimates; subjectProfile.TheData{currSubj}(currRun).data.EstimatedContrast(:)];
    end
    
    % Save out organized data structure
    subjectProfile.OrganizedData{currSubj} = ContrastData;
    % Save out mean location error
    subjectProfile.LocationError(currSubj) = mean(mean(locationErrorMat));
    %Save out mean contrast estimates
    for currField = 1:numel(dataFields)
        for currContrast = 1:numel(centerContrast)
            subjectProfile.ContrastEstimate(currSubj,currContrast,currField) = nanmean(ContrastData(currContrast).(dataFields{currField})(:,1),1);
            subjectProfile.ContrastEstimateOrder0(currSubj,currContrast,currField) = nanmean(ContrastData(currContrast).(dataFields{currField})...
                (ContrastData(currContrast).(dataFields{currField})(:,2) == 0),1); %Takes the average contrast of order 0 trials
            subjectProfile.ContrastEstimateOrder1(currSubj,currContrast,currField) = nanmean(ContrastData(currContrast).(dataFields{currField})...
                (ContrastData(currContrast).(dataFields{currField})(:,2) == 1),1); %Takes the average contrast of order 1 trials
            ContrastEstimate_PS(currSubj,currContrast,currField) = nanmean(ContrastData_PS(currContrast).(dataFields{currField})(:,1),1);
        end
    end

 
    
    % Generate Location Data Matrix
    allLocationError= [allLocationError; locationErrorMat(:)];
    Conditions = [Conditions; conditionMat(:)];
    
    % Individual Subject Plots
    if subjPlots
        % Location report distributions
        figure
        nbins = 100;
        histogram(locationErrorMat(:,1:4),nbins)
        if numel(subjectProfile.TheData{currSubj}) > 4
            hold on
            histogram(locationErrorMat(:,5:8),nbins)
        end
        title(['Location Error S' num2str(currSubj)])
        legend(num2str(subjectProfile.Order(currSubj)), num2str(~subjectProfile.Order(currSubj)))
        hold off
        
        % Center contrast vs perceived contrast
        plotContrasts = 100*round(centerContrast,2);
        figure
        cvp = loglog(plotContrasts, squeeze(subjectProfile.ContrastEstimate(currSubj,:,:)));
        hold on
        %         errorbar(plotContrasts,squeeze(subjectProfile.ContrastEstimate(currSubj,:,:)),)
        loglog(plotContrasts,plotContrasts,'--k')
        set(cvp, {'color'}, {[0 0 1]; [1 0 0]; [0 0 0]});
        set(cvp,'LineWidth',2)
        xticks(plotContrasts); yticks(plotContrasts);
        xticklabels({plotContrasts}); yticklabels({plotContrasts});
        xlabel('Center Contrast'); ylabel('Perceived Contrast');
        set(gca,'TickDir','out','XTick',plotContrasts,'YTick',plotContrasts,...
            'XTickLabel',plotContrasts,'YTickLabel',plotContrasts,...
            'TickDir','out');
        legend('Perception','Working Memory','Baseline')
        title(['Center vs. Perceived Contrasts S' num2str(currSubj)])
        hold off
    end
end

LocationMat = [allLocationError, contrastError, probeOffset];
offsets = -180:180;
offsets = offsets(:);

% mean contrast and location error for each probeOffset


% Calculate the mean and error for contrast estimates of each order



% Supression Index %
% Perception: first page.
suppressionIndex.Collapsed(1:numel(subjectProfile.SubjectName),:,1) = ((subjectProfile.ContrastEstimate(:,:,1)) - (subjectProfile.ContrastEstimate(:,:,3))) ...
    ./((subjectProfile.ContrastEstimate(:,:,1)) + (subjectProfile.ContrastEstimate(:,:,3)));
suppressionIndex.Order0(1:numel(subjectProfile.SubjectName),:,1) = ((subjectProfile.ContrastEstimateOrder0(:,:,1)) - (subjectProfile.ContrastEstimateOrder0(:,:,3))) ...
    ./((subjectProfile.ContrastEstimateOrder0(:,:,1)) + (subjectProfile.ContrastEstimateOrder0(:,:,3)));
suppressionIndex.Order1(1:numel(subjectProfile.SubjectName),:,1) = ((subjectProfile.ContrastEstimateOrder1(:,:,1)) - (subjectProfile.ContrastEstimateOrder1(:,:,3))) ...
    ./((subjectProfile.ContrastEstimateOrder1(:,:,1)) + (subjectProfile.ContrastEstimateOrder1(:,:,3)));
%Working Memory: second page.
suppressionIndex.Collapsed(1:numel(subjectProfile.SubjectName),:,2) = ((subjectProfile.ContrastEstimate(:,:,2)) - (subjectProfile.ContrastEstimate(:,:,3))) ...
    ./((subjectProfile.ContrastEstimate(:,:,2)) + (subjectProfile.ContrastEstimate(:,:,3)));
suppressionIndex.Order0(1:numel(subjectProfile.SubjectName),:,2) = ((subjectProfile.ContrastEstimateOrder0(:,:,2)) - (subjectProfile.ContrastEstimateOrder0(:,:,3))) ...
    ./((subjectProfile.ContrastEstimateOrder0(:,:,2)) + (subjectProfile.ContrastEstimateOrder0(:,:,3)));
suppressionIndex.Order1(1:numel(subjectProfile.SubjectName),:,2) = ((subjectProfile.ContrastEstimateOrder1(:,:,2)) - (subjectProfile.ContrastEstimateOrder1(:,:,3))) ...
    ./((subjectProfile.ContrastEstimateOrder1(:,:,2)) + (subjectProfile.ContrastEstimateOrder1(:,:,3)));




%% Plots
%%% Fit a gaussian to each plot
% Y=Amplitude*exp(-0.5*((X-Mean)/SD)^2)
%fminsearch(function,initial parameters, options, xdata, ydata)

% Super Subject Plots
if superPlots
    nbins = 100;
    %     % Total Location Error for baseline condition
    %     figure
    %     histfit(allLocationError(Conditions == 3),nbins,'normal')
    %     hold on
    %
    %     title(" 'Super Subject' Location Error")
    %     xlabel('Location Error'); ylabel('Frequency')
    %
    %     %Location Error split up by order, for baseline condition
    %     figure
    %     histfit(allLocationError(Conditions == 3 & Orders == 1),nbins)
    %     hold on
    %     histfit(allLocationError(Conditions == 3 & Orders == 0),nbins)
    %     title("'Super Subject' Location Error by Order")
    %     legend({'Loc-Con','Con-Loc'})
    %     xlabel('Location Error'); ylabel('Frequency')
    
    % Contrast estimate for baseline condition
    for currContrast = 1:numel(centerContrast)
        currInd = Conditions == 3 & Contrasts == centerContrast(currContrast);
        currData = contrastEstimates(currInd);
        currMean = mean(currData);
        currSD = std(currData);
        currGauss = exp(-0.5*((currData-currMean)/currSD).^2);
        
        figure
        histogram(currData,nbins);
        title([" 'Super Subject' Contrast Error (" num2str(centerContrast(currContrast)) "% Contrast)"])
        xlabel('Contrast Estimate'); ylabel('Frequency')
        figure; plot(currGauss,'k')
        
        
        %Contrast Error split up by order, for baseline condition
        %     figure
        %     h2 = histfit(log(currData),nbins);
        % %     h2(1).FaceColor = [0, 0.4470, 0.7410];
        %     h2(2).Color = [0, 0, 0];
        %     hold on
        %     h3 = histfit(log(currData),nbins);
        % %     h3(1).FaceColor = [0.6350, 0.0780, 0.1840];
        % %     h3(2).Color = [1 0 0];
        %     title(["'Super Subject' Contrast Error Order (" num2str(centerContrast(currContrast)) "% Contrast)"])
        %     legend({'Loc-Con','Loc-Con Fit','Con-Loc','Con-Loc Fit'})
        %     xlabel('Contrast Estimate'); ylabel('Frequency')
    end
    
    %     % Contrast Error as a function of probe offset, for baseline
    %     condition
    %     figure
    %     bar(offsets, );
    %     title('Probe Offset vs. Contrast Error')
    %     xlabel('Probe Offset'); ylabel('Contrast Error');
    %
    %     % Location Error as a function of probe offset, for baseline
    %     condition
    %     figure
    %     bar(offsets,);
    %     title('Probe Offset vs. Location Error')
    %     xlabel('Probe Offset'); ylabel('Location Error');
end




% %Group Plots
 if groupPlots
%     % Total Contrast Estimates
%     plotContrasts = 100*round(centerContrast,2);
%     figure
%     cvp_group = loglog(plotContrasts, 100.*squeeze(mean(subjectProfile.ContrastEstimate)));
%     hold on
%     errorbar(repmat(plotContrasts,1,3),100.*squeeze(mean(subjectProfile.ContrastEstimate)),100*squeeze(ContrastEstimate_ste),'k','LineStyle','None')
%     loglog(plotContrasts,plotContrasts,'--k')
%     set(cvp_group, {'color'}, {[0 0 1]; [1 0 0]; [0 0 0]});
%     set(cvp_group,'LineWidth',2)
%     xticks(plotContrasts); yticks(plotContrasts);
%     xticklabels({plotContrasts}); yticklabels({plotContrasts});
%     xlabel('Center Contrast'); ylabel('Perceived Contrast');
%     set(gca,'TickDir','out','XTick',plotContrasts,'YTick',plotContrasts,...
%         'XTickLabel',plotContrasts,'YTickLabel',plotContrasts,...
%         'TickDir','out');
%     legend('Perception','Working Memory','Baseline','Location','NorthWest')
%     title(['Center vs. Perceived Contrast'])
%     hold off
%
%     % Order 1 Contrast Estimates
%     figure
%     cvp_order1 = loglog(plotContrasts, 100.*());
%     hold on
%     errorbar(repmat(plotContrasts,1,3),100.*(),100.*(),'k','LineStyle','None')
%     loglog(plotContrasts,plotContrasts,'--k')
%     set(cvp_order1, {'color'}, {[0 0 1]; [1 0 0]; [0 0 0]});
%     set(cvp_order1,'LineWidth',2)
%     xticks(plotContrasts); yticks(plotContrasts);
%     xticklabels({plotContrasts}); yticklabels({plotContrasts});
%     xlabel('Center Contrast'); ylabel('Perceived Contrast');
%     set(gca,'TickDir','out','XTick',plotContrasts,'YTick',plotContrasts,...
%         'XTickLabel',plotContrasts,'YTickLabel',plotContrasts,...
%         'TickDir','out');
%     legend('Perception','Working Memory','Baseline','Location','NorthWest')
%     title(['Center vs. Perceived Contrast Order 1'])
%     hold off
%
%     % Order 2 Contrast Estimates
%     figure
%     cvp_order2 = loglog(plotContrasts, 100.*());
%     hold on
%     errorbar(repmat(plotContrasts,1,3),100.*(),100*(),'k','LineStyle','None')
%     loglog(plotContrasts,plotContrasts,'--k')
%     set(cvp_order2, {'color'}, {[0 0 1]; [1 0 0]; [0 0 0]});
%     set(cvp_order2,'LineWidth',2)
%     xticks(plotContrasts); yticks(plotContrasts);
%     xticklabels({plotContrasts}); yticklabels({plotContrasts});
%     xlabel('Center Contrast'); ylabel('Perceived Contrast');
%     set(gca,'TickDir','out','XTick',plotContrasts,'YTick',plotContrasts,...
%         'XTickLabel',plotContrasts,'YTickLabel',plotContrasts,...
%         'TickDir','out');
%     legend('Perception','Working Memory','Baseline','Location','NorthWest')
%     title(['Center vs. Perceived Contrast Order 2'])
%     hold off


% Total Contrast Estimates for Selective probes - exlcudes 0.75 and 0.1
% probe trials
    plotContrasts = 100*round(centerContrast,2);
    figure;
    cvp_group = loglog(plotContrasts, 100.*squeeze(mean(ContrastEstimate_PS)));
    hold on
    errorbar(repmat(plotContrasts,1,3),100.*squeeze(mean(ContrastEstimate_PS)),100*squeeze(std(ContrastEstimate_PS)./sqrt(size(ContrastEstimate_PS,1))),'k','LineStyle','None')
    loglog(plotContrasts,plotContrasts,'--k')
    set(cvp_group, {'color'}, {[0 0 1]; [1 0 0]; [0 0 0]});
    set(cvp_group,'LineWidth',2)
    xticks(plotContrasts); yticks(plotContrasts);
    xticklabels({plotContrasts}); yticklabels({plotContrasts});
    xlabel('Center Contrast'); ylabel('Perceived Contrast');
    set(gca,'TickDir','out','XTick',plotContrasts,'YTick',plotContrasts,...
        'XTickLabel',plotContrasts,'YTickLabel',plotContrasts,...
        'TickDir','out');
    legend('Perception','Working Memory','Baseline','Location','NorthWest')
    title(['Center vs. Perceived Contrast, without 0.75 or 0.1 Probe Trials'])
    hold off

%Suppression Index%
%Collapsed over both orders
hold off
figure('Color', [1 1 1]);
set(gcf, 'Name', sprintf('Suppression Index: Perception vs. Working Memory'));
hold all;
%Perception
errorbar(centerContrast', mean(suppressionIndex.Collapsed(:,:,1),1), std(suppressionIndex.Collapsed(:,:,1),1)...
    /sqrt(size(visualmemory_subjectsRan,2)), 'r','LineWidth',3);
%Working Memory
errorbar(centerContrast', mean(suppressionIndex.Collapsed(:,:,2),1), std(suppressionIndex.Collapsed(:,:,2),1)...
    /sqrt(size(visualmemory_subjectsRan,2)), 'b','LineWidth',3);
plot(repmat(centerContrast', [numel(subjectProfile.SubjectName) 1]), (suppressionIndex.Collapsed(:,:,1)),'r.','MarkerSize',10);
plot(repmat(centerContrast', [numel(subjectProfile.SubjectName) 1]), (suppressionIndex.Collapsed(:,:,2)),'b.','MarkerSize',10);
hold on
plot([0.09 0.8],[0 0],':k')
ylabel('Suppression Index (surround-nosurround)/(surround+nosurround)'); 
xlabel('Contrast (%)');
xlim([0.09 0.8]);
ylim([-0.3 0.3]);
axis square;
legend({'Perception Condition','Working Memory Condition'})

% Split between orders%
figure('Color', [1 1 1]);

subplot(1,2,1)
set(gcf, 'Name', sprintf('Suppression Index: Perception vs. Working Memory'));
hold all;
%Perception
errorbar(centerContrast', nanmean(suppressionIndex.Order0(:,:,1),1), nanstd(suppressionIndex.Order0(:,:,1),1)...
    /sqrt(numel(subjectProfile.SubjectName)), 'r','LineWidth',3)
%Working Memory
errorbar(centerContrast', nanmean(suppressionIndex.Order0(:,:,2),1), nanstd(suppressionIndex.Order0(:,:,2),1)...
    /sqrt(size(numel(subjectProfile.SubjectName)), 'b','LineWidth',3)); hold all;
plot(repmat(centerContrast', [numel(subjectProfile.SubjectName) 1]), (suppressionIndex.Order0(:,:,1)),'r.','MarkerSize',10);
plot(repmat(centerContrast', [numel(subjectProfile.SubjectName) 1]), (suppressionIndex.Order0(:,:,2)),'b.','MarkerSize',10);
hold all
plot([0.09 0.8],[0 0],':k')
ylabel('Suppression Index (surround-nosurround)/(surround+nosurround)'); 
xlabel('Contrast (%)');
xlim([0.09 0.8]);
ylim([-0.3 0.3]);
title('Order 0: Contrast, then Location');
axis square;
legend({'Perception Condition','Working Memory Condition'})

subplot(1,2,2)
hold all;
%Perception
errorbar(centerContrast', mean(suppressionIndex.Order1(:,:,1),1), std(suppressionIndex.Order1(:,:,1),1)...
    /sqrt(numel(subjectProfile.SubjectName)), 'r','LineWidth',3);
%Working Memory
errorbar(centerContrast', mean(suppressionIndex.Order1(:,:,2),1), std(suppressionIndex.Order1(:,:,2),1)...
    /sqrt(numel(subjectProfile.SubjectName)), 'b','LineWidth',3);
plot(repmat(centerContrast', [numel(subjectProfile.SubjectName) 1]), (suppressionIndex.Order1(:,:,1)),'r.','MarkerSize',10);
plot(repmat(centerContrast', [numel(subjectProfile.SubjectName) 1]), (suppressionIndex.Order1(:,:,2)),'b.','MarkerSize',10);
hold all
plot([0.09 0.8],[0 0],':k')
ylabel('Suppression Index (surround-nosurround)/(surround+nosurround)'); 
xlabel('Contrast (%)');
xlim([0.09 0.8]);
ylim([-0.3 0.3]);
title('Order 1: Location, then Contrast');
axis square;
legend({'Perception Condition','Working Memory Condition'})




 end
