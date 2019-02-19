%%% vmmf_subject_analysis.m
%% Prepare data and working directories
clear all; close all; clc;
scrsz = get(groot,'ScreenSize');

subjPlots = 0;
groupPlots = 1;
superPlots = 1;

expDir=pwd;
dataDir='data_master';
subjects = 1:10;
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

% Organize Contrast and Location  Estimates
for currSubj = 1:numel(subjectProfile.SubjectName)
    ContrastData = struct('Perception',{[] [] [] [] []},'WorkingMemory',{[] [] [] [] []},'Baseline',{[] [] [] [] []}); % 5 slots for each contrast
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

%% Probe Effect Analysis
probeMat = [allLocationError, contrastError, probeOffset];
offsets = -180:180;
offsets = offsets(:);
offsetRange  = 45/2;

for nContrast = 1:numel(centerContrast)
    currContrast = centerContrast(nContrast);
    currIndx = Contrasts == currContrast & Conditions == 3;
    currProbeOffset = probeOffset(currIndx);
    currLocationError = allLocationError(currIndx);
    currContrastError = contrastError(currIndx);
    
    
end

%% Distribution Analysis
% Fit Gaussian distribution to all participant data
% Y=Amplitude*exp(-0.5*((X-Mean)/SD)^2)
for nCondition = 1:numel(dataFields)
    for nContrast = 1:numel(centerContrast)
        currInd = Conditions == nCondition & Contrasts == centerContrast(nContrast);
        currData = contrastEstimates(currInd);
        edges_con = 0:0.05:1;
        [N_bl_con(nContrast,:),x_con] = histcounts(currData,edges_con);
        N_bl_con(nContrast,:) = N_bl_con(nContrast,:)./max(N_bl_con(nContrast,:));
        xvalues_con = (x_con(1:end-1)+x_con(2:end))/2;
    end
    
    % fit all subjects data simultaneous
    % mean, amplitude, sigma
    startValues = [centerContrast' repmat(0.1, [1 numel(centerContrast)])];
    options = optimset('MaxFunEvals', 5000.*numel(startValues), 'MaxIter', 5000.*numel(startValues));
    
    [est_params_tmp, r2_bl] = fminsearch('mygauss_allContrasts', startValues, options, N_bl_con, xvalues_con);
    tmp = reshape(est_params_tmp, [5 2])';
    est_params_bl = [tmp(1,:); ones(1,5); tmp(2,:)];
    
    % Get estimates for vwm when both mean and width are allowed to freely vary
    
    
    % When the width is fixed by taking the estimates from the perception
    % condtion - does this explain the data as good as allow it to vary?
    
    % Generate gaussian fits from estimated params
    
    figure('Position',[50 50 scrsz(3) scrsz(3)/2], 'color', [1 1 1], 'name', [dataFields{nCondition} ' Contrast Estimates'])
    for nContrast = 1:numel(centerContrast)
        x_fit = 0.01:0.01:1;
        y_est_bl = est_params_bl(2,nContrast)*exp(-(x_fit-est_params_bl(1,nContrast)).^2/(2*(est_params_bl(3,nContrast)^2)));
        xvalues_con = (x_con(1:end-1)+x_con(2:end))/2;
        subplot(1,5,nContrast)
        bar(xvalues_con, N_bl_con(nContrast,:));
        hold on
        plot(x_fit, y_est_bl, 'r', 'LineWidth', 2)
        box off; title([num2str(round(100*centerContrast(nContrast))) '%']); xlabel('Contrast (%)'); ylabel('Normalized count of perceived contrast');
        axis square
        %plot parameters alone 
    end
end
%% Supression Index
% Perception: first page.
suppressionIndex.Collapsed(1:size(visualmemory_subjectsRan,2),:,1) = ((subjectProfile.ContrastEstimate(:,:,1)) - (subjectProfile.ContrastEstimate(:,:,3))) ...
    ./((subjectProfile.ContrastEstimate(:,:,1)) + (subjectProfile.ContrastEstimate(:,:,3)));
suppressionIndex.Order0(1:size(visualmemory_subjectsRan,2),:,1) = ((subjectProfile.ContrastEstimateOrder0(:,:,1)) - (subjectProfile.ContrastEstimateOrder0(:,:,3))) ...
    ./((subjectProfile.ContrastEstimateOrder0(:,:,1)) + (subjectProfile.ContrastEstimateOrder0(:,:,3)));
suppressionIndex.Order1(1:size(visualmemory_subjectsRan,2),:,1) = ((subjectProfile.ContrastEstimateOrder1(:,:,1)) - (subjectProfile.ContrastEstimateOrder1(:,:,3))) ...
    ./((subjectProfile.ContrastEstimateOrder1(:,:,1)) + (subjectProfile.ContrastEstimateOrder1(:,:,3)));
%Working Memory: second page.
suppressionIndex.Collapsed(1:size(visualmemory_subjectsRan,2),:,2) = ((subjectProfile.ContrastEstimate(:,:,2)) - (subjectProfile.ContrastEstimate(:,:,3))) ...
    ./((subjectProfile.ContrastEstimate(:,:,2)) + (subjectProfile.ContrastEstimate(:,:,3)));
suppressionIndex.Order0(1:size(visualmemory_subjectsRan,2),:,2) = ((subjectProfile.ContrastEstimateOrder0(:,:,2)) - (subjectProfile.ContrastEstimateOrder0(:,:,3))) ...
    ./((subjectProfile.ContrastEstimateOrder0(:,:,2)) + (subjectProfile.ContrastEstimateOrder0(:,:,3)));
suppressionIndex.Order1(1:size(visualmemory_subjectsRan,2),:,2) = ((subjectProfile.ContrastEstimateOrder1(:,:,2)) - (subjectProfile.ContrastEstimateOrder1(:,:,3))) ...
    ./((subjectProfile.ContrastEstimateOrder1(:,:,2)) + (subjectProfile.ContrastEstimateOrder1(:,:,3)));


%% Plots

% Super Subject Plots
if superPlots
    % Total Location Error for baseline condition
    for nCondition = 1:numel(dataFields)
        figure('Position',[50 50 scrsz(3) scrsz(3)/2], 'color', [1 1 1], 'name', [dataFields{nCondition} ' Location Error'])
        for nContrast = 1:numel(centerContrast)
            currInd = Conditions == nCondition & Contrasts == centerContrast(nContrast);
            edges_loc = -180:180;
            [N_bl_loc(nContrast,:),x_loc] = histcounts(allLocationError(currInd),edges_loc);
            N_bl_loc(nContrast,:) = N_bl_loc(nContrast,:)./max(N_bl_loc(nContrast,:));
            xvalues_loc = (x_loc(1:end-1)+x_loc(2:end))/2;
            
            subplot(1,5,nContrast)
            bar(xvalues_loc, N_bl_loc(nContrast,:));
            hold on
            title([num2str(round(100*centerContrast(nContrast))) '%'])
            xlabel('Location Error'); ylabel('Frequency')
        end
    end
    
    %     %Location Error split up by order, for baseline condition
    %     figure
    %     histfit(allLocationError(Conditions == 3 & Orders == 1),nbins)
    %     hold on
    %     histfit(allLocationError(Conditions == 3 & Orders == 0),nbins)
    %     title("'Super Subject' Location Error by Order")
    %     legend({'Loc-Con','Con-Loc'})
    %     xlabel('Location Error'); ylabel('Frequency')
end


%%Group Plots
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
    plot(repmat(centerContrast', [size(visualmemory_subjectsRan,2) 1]), (suppressionIndex.Collapsed(:,:,1)),'r.','MarkerSize',10);
    plot(repmat(centerContrast', [size(visualmemory_subjectsRan,2) 1]), (suppressionIndex.Collapsed(:,:,2)),'b.','MarkerSize',10);
    hold on
    plot([0.09 0.8],[0 0],':k')
    ylabel('Suppression Index (surround-nosurround)/(surround+nosurround)');
    xlabel('Contrast (%)');
    xlim([0.09 0.8]);
    ylim([-0.15 0.15]);
    axis square;
    legend({'Perception Condition','Working Memory Condition'})
    
    % Split between orders%
    figure('Color', [1 1 1]);
    
    subplot(1,2,1)
    set(gcf, 'Name', sprintf('Suppression Index: Perception vs. Working Memory'));
    hold all;
    %Perception
    errorbar(centerContrast', nanmean(suppressionIndex.Order0(:,:,1),1), nanstd(suppressionIndex.Order0(:,:,1),1)...
        /sqrt(size(visualmemory_subjectsRan,2)), 'r','LineWidth',3)
    %Working Memory
    errorbar(centerContrast', nanmean(suppressionIndex.Order0(:,:,2),1), nanstd(suppressionIndex.Order0(:,:,2),1)...
        /sqrt(size(visualmemory_subjectsRan,2)), 'b','LineWidth',3); hold all;
    plot(repmat(centerContrast', [size(visualmemory_subjectsRan,2) 1]), (suppressionIndex.Order0(:,:,1)),'r.','MarkerSize',10);
    plot(repmat(centerContrast', [size(visualmemory_subjectsRan,2) 1]), (suppressionIndex.Order0(:,:,2)),'b.','MarkerSize',10);
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
        /sqrt(size(visualmemory_subjectsRan,2)), 'r','LineWidth',3);
    %Working Memory
    errorbar(centerContrast', mean(suppressionIndex.Order1(:,:,2),1), std(suppressionIndex.Order1(:,:,2),1)...
        /sqrt(size(visualmemory_subjectsRan,2)), 'b','LineWidth',3);
    plot(repmat(centerContrast', [size(visualmemory_subjectsRan,2) 1]), (suppressionIndex.Order1(:,:,1)),'r.','MarkerSize',10);
    plot(repmat(centerContrast', [size(visualmemory_subjectsRan,2) 1]), (suppressionIndex.Order1(:,:,2)),'b.','MarkerSize',10);
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
