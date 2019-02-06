%%% vmmf_subject_analysis.m
%% Prepare data and working directories
clear all; close all; clc;

subjPlots = 0;
groupPlots = 0;
superPlots = 0;

expDir=pwd;
dataDir='data_master';

% Load in data
cd(dataDir)
load('visualmemory_subjectsRan')
subjectProfile=struct('SubjectName',[] ,'Order', [], 'Condition', [], 'TheData',[],'OrganizedData',[],'LocationError', [],'ContrastEstimate',[]);

subjectProfile.SubjectName = cell(1,10); % initialize subject names
subjectProfile.Order = nan(1,10); % initialize order
subjectProfile.TheData = cell(1,10); % intialize theData holders
subjectProfile.OrganizedData = cell(1,10); % initialized organize data
subjectProfile.LocationError = nan(1,10); % initialize location error (Average location error for each subject)
subjectProfile.ContrastEstimate = nan(10,5,3); % initialize contrast estimate (concactanated contrast estimates for all subjects, contrasts, and conditions)

subjectProfile.SubjectName = cellfun(@str2double,{visualmemory_subjectsRan{1,:}});
for currSubj=1:numel(subjectProfile.SubjectName)
    if exist(['data_vmmf_00' num2str(subjectProfile.SubjectName(currSubj)) '.mat'],'file') ~= 0
        load(['data_vmmf_00' num2str(subjectProfile.SubjectName(currSubj))]);
    end
    subjectProfile.Order(currSubj) = strcmp(visualmemory_subjectsRan(2,currSubj),'a'); % report order
    subjectProfile.Condition(:,currSubj) = theData(1).p.trialSchedule; %subject condition schedule
    subjectProfile.TheData{currSubj} = theData; %subject data
end
cd(expDir)
%% Data Organization
centerContrast = unique(subjectProfile.TheData{1}(1).p.trialEvents(:,3));
allLocationData = [];
allConditions = [];
allOrder = [];
order1LocationData = [];
order2LocationData = [];
order1ContrastEstimate = struct('Data',[],'Mean',[],'STD',[]);
order2ContrastEstimate = struct('Data',[],'Mean',[],'STD',[]);
order1ContrastEstimate.Data = nan(10,5,3);
order2ContrastEstimate.Data = nan(10,5,3);
probeOffset = struct('Data', [],'Mean',[],'STD',[]);
contrastError =  struct('Data', [],'Mean',[],'STD',[]);

% allData(currContrast).(dataFields{currField})(currRun,3)  = subjectProfile.TheData{currSubj}(currRun).p.trialEvents(relevantTrials,2) - subjectProfile.TheData{currSubj}(currRun).data.EstimatedLocation(relevantTrials); %Location Error
% allData(currContrast).(dataFields{currField})(currRun,4) =  subjectProfile.TheData{currSubj}(currRun).p.trialEvents(relevantTrials,2) - subjectProfile.TheData{currSubj}(currRun).p.trialEvents(relevantTrials,4); %Probe Offset (target - probe)


% Organize Contrast and Location  Estimates
for currSubj = 1:numel(subjectProfile.SubjectName)
    ContrastData = struct('Perception',{[] [] [] [] []},'WorkingMemory',{[] [] [] [] []},'Baseline',{[] [] [] [] []}); % 5 slots for each contrast
    dataFields = fieldnames(ContrastData);
    locationData = nan(length(subjectProfile.TheData{currSubj}(1).p.trialEvents),numel(subjectProfile.TheData{currSubj}));
    conditionMat = nan(length(subjectProfile.TheData{currSubj}(1).p.trialEvents),numel(subjectProfile.TheData{currSubj}));
    for currRun = 1:numel(subjectProfile.TheData{currSubj})
        %grab the condition of each trial in current run
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
            end
        end
        
        % Reformat location data and determine probe offset from target
        locationData(:,currRun) = subjectProfile.TheData{currSubj}(currRun).p.trialEvents(:,2) - subjectProfile.TheData{currSubj}(currRun).data.EstimatedLocation(:);
        locationData(locationData > 180) = locationData(locationData > 180) - 360;
        locationData(locationData < -180) = locationData(locationData < -180) + 360;
        probeOffset.Data = [probeOffset.Data; subjectProfile.TheData{currSubj}(currRun).p.trialEvents(:,2)-subjectProfile.TheData{currSubj}(currRun).p.trialEvents(:,4)];
        contrastError.Data = [contrastError.Data; subjectProfile.TheData{currSubj}(currRun).data.DifferenceContrast(:)];
    end
    
    % Save out organized data structure 
    subjectProfile.OrganizedData{currSubj} = ContrastData;
    % Save out mean location error 
    subjectProfile.LocationError(currSubj) = mean(mean(locationData));
    %Save out mean contrast estimates
    for currField = 1:numel(dataFields)
        for currContrast = 1:numel(centerContrast)
            subjectProfile.ContrastEstimate(currSubj,currContrast,currField) = nanmean(ContrastData(currContrast).(dataFields{currField})(:,1),1);
        end
    end
    
    % Split location data by order
    allLocationData= [allLocationData; locationData(:)];
    allConditions = [allConditions; conditionMat(:)];
    allOrder = [allOrder; repmat(currOrder,numel(locationData),1)];
    
    if subjectProfile.Order(currSubj)
        order1LocationData = [order1LocationData; locationData(:)];
    else
        order2LocationData = [order2LocationData; locationData(:)];
    end
       
    if subjPlots
        % Location report distributions
        figure
        nbins = 100;
        histogram(locationData(:,1:4),nbins)
        if numel(subjectProfile.TheData{currSubj}) > 4
            hold on
            histogram(locationData(:,5:8),nbins)
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

%Format probeOffset space to be -180 to 180
probeOffset.Data(probeOffset.Data > 180) = probeOffset.Data(probeOffset.Data > 180) - 360;
probeOffset.Data(probeOffset.Data < -180) = probeOffset.Data(probeOffset.Data < -180) + 360;
offsets = -180:180;
offsets = offsets(:);

% Initialize contrast and location error structures
contrastError.Mean = nan(length(unique(probeOffset.Data)),1);
locationError.Mean = nan(length(unique(probeOffset.Data)),1);
contrastError.STD = nan(length(unique(probeOffset.Data)),1);
locationError.STD = nan(length(unique(probeOffset.Data)),1);

% Calculate contrast and location error avgs and stds
for currOffset = 1:length(offsets)
    contrastError.Mean(currOffset) = mean(contrastError.Data(probeOffset.Data == offsets(currOffset)));
    locationError.Mean(currOffset) = mean(allLocationData(probeOffset.Data == offsets(currOffset)));
    contrastError.STD(currOffset) = std(contrastError.Data(probeOffset.Data == offsets(currOffset)));
    locationError.STD(currOffset) = std(allLocationData(probeOffset.Data == offsets(currOffset)));
end

% Initialize order structure variables
order1ContrastEstimate.Mean = nanmean(order1ContrastEstimate.Data);
order1ContrastEstimate.STE = nanstd(order2ContrastEstimate.Data)/numel(subjectProfile.SubjectName);
order2ContrastEstimate.Mean = nanmean(order2ContrastEstimate.Data);
order2ContrastEstimate.STE = nanstd(order2ContrastEstimate.Data)/numel(subjectProfile.SubjectName);

% Separate Contrast Estimates by Order
ContrastEstimate_ste  = std(subjectProfile.ContrastEstimate)/length(subjectProfile.SubjectName);
for currSubj = 1:numel(subjectProfile.SubjectName)
    for currField = 1:numel(dataFields)
        for currContrast = 1:numel(centerContrast)
            order1ContrastEstimate.Data(currSubj,currContrast,currField) = nanmean(subjectProfile.OrganizedData{currSubj}(currContrast).(dataFields{currField})(subjectProfile.OrganizedData{currSubj}(currContrast).(dataFields{currField})(:,2) == 1,1));
            order2ContrastEstimate.Data(currSubj,currContrast,currField) = nanmean(subjectProfile.OrganizedData{currSubj}(currContrast).(dataFields{currField})(subjectProfile.OrganizedData{currSubj}(currContrast).(dataFields{currField})(:,2) == 0,1));
        end
    end
end
%% Plots

% Super Subject Plots
if superPlots
    nbins = 100;
    % Total Location Error
    figure
    histogram(allLocationData,nbins)
    title(" 'Super Subject' Location Error")
    xlabel('Location Error'); ylabel('Frequency')
    
    %Location Error split up by order
    figure
    histogram(order1LocationData,nbins)
    hold on
    histogram(order2LocationData,nbins)
    title("'Super Subject' Location Error by Order")
    legend({'Loc-Con','Con-Loc'})
    xlabel('Location Error'); ylabel('Frequency')
    
    % Contrast Error as a function of probe offset
    figure
    bar(offsets, contrastError.Mean);
    hold on
    plot(offsets,repmat(mean(contrastError.Mean),length(offsets),1),'k','LineWidth',2)
    title('Probe Offset vs. Contrast Error')
    xlabel('Probe Offset'); ylabel('Contrast Error');
    
    % Location Error as a function of probe offset
    figure
    bar(offsets, locationError.Mean);
    hold on
    plot(offsets,repmat(mean(locationError.Mean),length(offsets),1),'k','LineWidth',2)
    title('Probe Offset vs. Location Error')
    xlabel('Probe Offset'); ylabel('Location Error');
end

%Group Plots
if groupPlots
    % Total Contrast Estimates
    plotContrasts = 100*round(centerContrast,2);
    figure
    cvp_group = loglog(plotContrasts, 100.*squeeze(mean(subjectProfile.ContrastEstimate)));
    hold on
    errorbar(repmat(plotContrasts,1,3),100.*squeeze(mean(subjectProfile.ContrastEstimate)),100*squeeze(ContrastEstimate_ste),'k','LineStyle','None')
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
    title(['Center vs. Perceived Contrast'])
    hold off
    
    % Order 1 Contrast Estimates
    figure
    cvp_order1 = loglog(plotContrasts, 100.*squeeze(order1ContrastEstimate.mean));
    hold on
    errorbar(repmat(plotContrasts,1,3),100.*squeeze(order1ContrastEstimate.mean),100*squeeze(order1ContrastEstimate.ste),'k','LineStyle','None')
    loglog(plotContrasts,plotContrasts,'--k')
    set(cvp_order1, {'color'}, {[0 0 1]; [1 0 0]; [0 0 0]});
    set(cvp_order1,'LineWidth',2)
    xticks(plotContrasts); yticks(plotContrasts);
    xticklabels({plotContrasts}); yticklabels({plotContrasts});
    xlabel('Center Contrast'); ylabel('Perceived Contrast');
    set(gca,'TickDir','out','XTick',plotContrasts,'YTick',plotContrasts,...
        'XTickLabel',plotContrasts,'YTickLabel',plotContrasts,...
        'TickDir','out');
    legend('Perception','Working Memory','Baseline','Location','NorthWest')
    title(['Center vs. Perceived Contrast Order 1'])
    hold off
    
    % Order 2 Contrast Estimates
    figure
    cvp_order2 = loglog(plotContrasts, 100.*squeeze(order2ContrastEstimate.mean));
    hold on
    errorbar(repmat(plotContrasts,1,3),100.*squeeze(order2ContrastEstimate.mean),100*squeeze(order2ContrastEstimate.ste),'k','LineStyle','None')
    loglog(plotContrasts,plotContrasts,'--k')
    set(cvp_order2, {'color'}, {[0 0 1]; [1 0 0]; [0 0 0]});
    set(cvp_order2,'LineWidth',2)
    xticks(plotContrasts); yticks(plotContrasts);
    xticklabels({plotContrasts}); yticklabels({plotContrasts});
    xlabel('Center Contrast'); ylabel('Perceived Contrast');
    set(gca,'TickDir','out','XTick',plotContrasts,'YTick',plotContrasts,...
        'XTickLabel',plotContrasts,'YTickLabel',plotContrasts,...
        'TickDir','out');
    legend('Perception','Working Memory','Baseline','Location','NorthWest')
    title(['Center vs. Perceived Contrast Order 2'])
    hold off
end
