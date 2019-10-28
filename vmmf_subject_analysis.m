%%% vmmf_subject_analysis.m
%% Prepare data and working directories
clear variables; close all; clc;
scrsz = get(groot,'ScreenSize');

subjPlots = 1;
groupPlots = 1;
superPlots = 1;
normalizationModelFitPlots = 0;

perceptionIndex = 1;
workingMemoryIndex = 2;
baselineIndex = 3;

expDir=pwd;
dataDir='data_master';
subjects = [1:9];
% Load in data
cd(dataDir)
load('visualmemory_subjectsRan')
subjectProfile=struct('SubjectName',[] ,'Order', [], 'Condition', [], 'TheData',[],'OrganizedData',[],'LocationError', [],'ContrastEstimate',[]);

subjectProfile.SubjectName = subjects;
subjectProfile.Order = nan(1,length(subjects)); % initialize report order
subjectProfile.TheData = cell(1,length(subjects)); % intialize theData holders
subjectProfile.OrganizedData = cell(1,length(subjects)); % initialized organize data
subjectProfile.LocationError = nan(1,length(subjects)); % initialize location error (Average location error for each subject)
subjectProfile.ContrastEstimate = nan(length(subjects),5,3); % initialize contrast estimate (Average contrast estimate for each condition)

% Load in theData
totalNumTrials = 0;
for ns=1:numel(subjectProfile.SubjectName)
    if exist(['data_vmmf_' num2str(subjectProfile.SubjectName(ns)) '.mat'],'file') ~= 0
        load(['data_vmmf_' num2str(subjectProfile.SubjectName(ns))]);
    end
    subjectProfile.Order(ns) = strcmp(visualmemory_subjectsRan(2,ns),'a'); % report order (1=location-contrast; 0=contrast-location)
    subjectProfile.Condition(:,ns) = theData(1).p.trialSchedule; %subject condition schedule
    subjectProfile.TheData{ns} = theData; %subject data
    totalNumTrials = totalNumTrials + (numel(subjectProfile.TheData{ns})*subjectProfile.TheData{ns}(1).p.numTrials);
end
cd(expDir)

%% Data Organization
centerContrast = unique(subjectProfile.TheData{1}(1).p.trialEvents(:,3));
plotContrasts = 100.*round(centerContrast,2);
allLocationError = [];
Conditions = [];
Orders = [];
Contrasts = [];
probeOffset = [];
contrastError = [];
contrastEstimates = [];
pooledContrastEstimatesMat = [];
index_probeSelective = [];
locationErrorMatSuper = [];
conditionMatSuper = [];
locationAveragesByCond = zeros(3,length(subjects));

subplotDIMs = subjects(~(rem(length(subjects), subjects)));

if mod(length(subplotDIMs),2) ~= 0
    subplotX = median(subplotDIMs);
    subplotY=subplotX;
else
    subplotX = subplotDIMs(length(subplotDIMs)/2);
    subplotY = subplotDIMs(length(subplotDIMs)/2+1);
end

% Organize Contrast and Location  Estimates
for ns = 1:numel(subjectProfile.SubjectName)
    structSize = cell(1,5);
    ContrastData = struct('Perception',structSize,'WorkingMemory',structSize,'P_Baseline',structSize,'W_Baseline',structSize); % 5 slots for each contrast
    ContrastData_PS = struct('Perception',structSize,'WorkingMemory',structSize,'P_Baseline',structSize,'W_Baseline',structSize); % 5 slots for each contrast; PS = probe selective.
    dataFields = fieldnames(ContrastData); %Store contrastData field names
    
    % Initialize location error and condition tracking matrices
    locationError = nan(length(subjectProfile.TheData{ns}(1).p.trialEvents),numel(subjectProfile.TheData{ns}));
    conditionMat = nan(length(subjectProfile.TheData{ns}(1).p.trialEvents),numel(subjectProfile.TheData{ns}));
    
    for currRun = 1:numel(subjectProfile.TheData{ns})
        % Grab the condition of each trial in current run
        conditionMat(:,currRun) = subjectProfile.TheData{ns}(currRun).p.trialEvents(:,1);
        
        % Keep track of current order
        if currRun <= 4
            runIndx = currRun;
            currOrder = subjectProfile.Order(ns);
        else
            runIndx = currRun-4;
            currOrder = ~subjectProfile.Order(ns);
        end
        
        % Go through each field in allData structure (Perception, Working
        % Memory, Baseline) and get the mean estimate for a contrast level
        for nf = 1:numel(dataFields)
            for nc = 1:numel(centerContrast)
                if nf < 3 % this is a simple.temporary fix to bypass the fact that there is no "condition 4" to distinguish P_Baseline from W_Baseline
                    currCond = nf;
                    relevantTrials = subjectProfile.TheData{ns}(currRun).p.trialEvents(:,1) == currCond & subjectProfile.TheData{ns}(currRun).p.trialEvents(:,3) == centerContrast(nc);
                else
                    currCond = 3; % once the "W_Basline" field is reached, still use '3' as the current condition rather then currField (which would be 4 when this is accessed)
                    relevantTrials = subjectProfile.TheData{ns}(currRun).p.trialEvents(:,1) == currCond & subjectProfile.TheData{ns}(currRun).p.trialEvents(:,3) == centerContrast(nc);
                    if ~(sum(relevantTrials) > 0 && sum(conditionMat(:,currRun) == nf-2))
                        relevantTrials(:) = 0;
                    end
                end
                
                if sum(relevantTrials) > 0
                    ContrastData(nc).(dataFields{nf}) = ...
                        [ContrastData(nc).(dataFields{nf}); ...
                        [subjectProfile.TheData{ns}(currRun).data.EstimatedContrast(relevantTrials)' currOrder*ones(sum(relevantTrials),1)]];
                end
                
                %Take probe selective contrast data - probes not equal to
                %0.1 or 0.75
                relevantTrials_probeSelective = subjectProfile.TheData{ns}(currRun).p.trialEvents(:,1) == currCond ...
                    & subjectProfile.TheData{ns}(currRun).p.trialEvents(:,3) == centerContrast(nc) & subjectProfile.TheData{ns}(currRun).p.trialEvents(:,5) ~= 0.75 ...
                    & subjectProfile.TheData{ns}(currRun).p.trialEvents(:,5) ~= 0.1; %Takes same data, minus trials that include 0.1 or 0.75 contrast
                if sum(relevantTrials_probeSelective) > 0
                    ContrastData_PS(nc).(dataFields{nf}) = ...
                        [ContrastData_PS(nc).(dataFields{nf}); ...
                        [subjectProfile.TheData{ns}(currRun).data.EstimatedContrast(relevantTrials_probeSelective)' currOrder*ones(sum(relevantTrials_probeSelective),1)]];
                end
            end
        end
        
        % Reformat location data and determine probe offset from target
        locationError(:,currRun) = subjectProfile.TheData{ns}(currRun).p.trialEvents(:,2) - subjectProfile.TheData{ns}(currRun).data.EstimatedLocation(:);
        
        % Format location space to be -180 to 180
        locationError(locationError > 180) = locationError(locationError > 180) - 360;
        locationError(locationError < -180) = locationError(locationError < -180) + 360;
        probeOffset = [probeOffset; subjectProfile.TheData{ns}(currRun).p.trialEvents(:,2)-subjectProfile.TheData{ns}(currRun).p.trialEvents(:,4)];
        probeOffset(probeOffset > 180) = probeOffset(probeOffset > 180) - 360;
        probeOffset(probeOffset < -180) = probeOffset(probeOffset < -180) + 360;
        
        % Store contrast error
        contrastError = [contrastError; subjectProfile.TheData{ns}(currRun).data.DifferenceContrast(:)];
        Orders = [Orders; repmat(currOrder,numel(relevantTrials),1)];
        Contrasts = [Contrasts; subjectProfile.TheData{ns}(currRun).p.trialEvents(:,3)];
        contrastEstimates = [contrastEstimates; subjectProfile.TheData{ns}(currRun).data.EstimatedContrast(:)];
        
        % Store separate location error
        if currOrder % is location-contrast
            locationError1(:,runIndx,ns) = locationError(:,currRun);
        else %contrast-location
            locationError0(:,runIndx,ns) = locationError(:,currRun);
        end
        
    end
    
    % Save out average location error Per Condition Per Subject
    for i = 1:3
        locationAveragesByCond(i,ns) = mean(mean(abs(locationError(conditionMat==i))));
    end
    
    % Save out organized data structure
    subjectProfile.OrganizedData{ns} = ContrastData;
    
    % Save out mean location error
    subjectProfile.LocationError(ns) = mean(mean(abs(locationError)));
    
    %Save out mean contrast estimates
    for nf = 1:numel(dataFields)
        for nc = 1:numel(centerContrast)
            subjectProfile.ContrastEstimate(ns,nc,nf) = nanmean(ContrastData(nc).(dataFields{nf})(:,1),1);
            subjectProfile.ContrastEstimateOrder0(ns,nc,nf) = nanmean(ContrastData(nc).(dataFields{nf})...
                (ContrastData(nc).(dataFields{nf})(:,2) == 0),1); %Takes the average contrast of order 0 trials
            subjectProfile.ContrastEstimateOrder1(ns,nc,nf) = nanmean(ContrastData(nc).(dataFields{nf})...
                (ContrastData(nc).(dataFields{nf})(:,2) == 1),1); %Takes the average contrast of order 1 trials
            ContrastEstimate_PS(ns,nc,nf) = nanmean(ContrastData_PS(nc).(dataFields{nf})(:,1),1);
            contrastEstimatesMat(:,nc,nf,ns) = ContrastData(nc).(dataFields{nf})(:,1);
        end
    end
    pooledContrastEstimatesMat = [pooledContrastEstimatesMat; contrastEstimatesMat(:,:,:,ns)];
    
    % Generate Location Data Matrix
    allLocationError= [allLocationError; locationError(:)];
    Conditions = [Conditions; conditionMat(:)];
    
    % Test whether baseline location estimates come from different distributions
    [h_loc(ns,:), p_loc(ns,:)] = ttest(mean(locationError0,3),mean(locationError1,3));
    
    % Test whether baseline contrast estimates come from different distributions
    for nc = 1:numel(centerContrast)
        P_Baseline(:,nc) = ContrastData(nc).P_Baseline(:,1);
        W_Baseline(:,nc) = ContrastData(nc).W_Baseline(:,1);
    end
    [h_bl(ns,:),p_bl(ns,:)] = ttest(P_Baseline,W_Baseline);
    
    %set(cvp, {'Color'}, {[0 0 1]; [1 0 0]; [0 0 0]});
    % Individual Subject Plots
    nbins = 100;
    if subjPlots
        % Location report distributions (comparing orders here)
        figure(1)
        subplot(subplotX,subplotY,ns)
        if subjectProfile.Order(ns)
            histogram(locationError(:,1:4),nbins,'Normalization','count','BinWidth',5)
            if numel(subjectProfile.TheData{ns}) > 4
                hold on
                histogram(locationError(:,5:8),nbins,'Normalization','count','BinWidth',5) % plot second order
            end
            box off; set(gca, 'TickDir','out','ColorOrder',[1 0 0; 0 0 1]); xlim([-180 180]);
            title(['S' num2str(ns)]);
        elseif ~subjectProfile.Order(ns)
            histogram(locationError(:,5:8),nbins,'Normalization','count','BinWidth',5)
            if numel(subjectProfile.TheData{ns}) > 4
                hold on
                histogram(locationError(:,1:4),nbins,'Normalization','count','BinWidth',5) % plot second order
            end
            box off; set(gca, 'TickDir','out','ColorOrder',[1 0 0; 0 0 1]); xlim([-180 180]);
            title(['S' num2str(ns)]);
        end
        
        % Center contrast vs perceived contrast
        figure(2)
        subplot(subplotX,subplotY,ns)
        errorbar(plotContrasts,100.*subjectProfile.ContrastEstimate(ns,:,1),100.*std(contrastEstimatesMat(:,:,1,ns))./sqrt(numel(subjectProfile.TheData{ns})),'.-','CapSize',0);
        hold on
        errorbar(plotContrasts,100.*subjectProfile.ContrastEstimate(ns,:,2),100.*std(contrastEstimatesMat(:,:,2,ns))./sqrt(numel(subjectProfile.TheData{ns})),'.-','CapSize',0);
        errorbar(plotContrasts,100.*subjectProfile.ContrastEstimate(ns,:,3),100.*std(contrastEstimatesMat(:,:,3,ns))./sqrt(numel(subjectProfile.TheData{ns})),'.--','CapSize',0);
        errorbar(plotContrasts,100.*subjectProfile.ContrastEstimate(ns,:,4),100.*std(contrastEstimatesMat(:,:,4,ns))./sqrt(numel(subjectProfile.TheData{ns})),'.--','CapSize',0);
        line([plotContrasts(1) plotContrasts(end)],[plotContrasts(1) plotContrasts(end)],'Color','k')
        xticks(plotContrasts); yticks(plotContrasts);
        xticklabels({plotContrasts}); yticklabels({plotContrasts});
        set(gca,'XScale','log','YScale','log','ColorOrder',[1 0 0; 0 0 1; 0.5 0 0; 0 0 0.5],'TickDir','out','XTick',plotContrasts,'YTick',plotContrasts,...
            'XTickLabel',plotContrasts,'YTickLabel',plotContrasts,...
            'TickDir','out'); box off; xlim([min(plotContrasts) max(plotContrasts)]);
        title(['S' num2str(ns)])
    end
    
end

if subjPlots
    figure(1)
    legend('Loc-Con', 'Con-Loc')
    xlabel('Location Error'); ylabel('Count');
    %     set(gcf, 'Position',scrsz);
    suptitle('Location Error')
    
    figure(2)
    legend('Sim.','Seq.','Sim. BL','Seq. BL','Location','NorthWest')
    suptitle('Presented v. Perceived Contrast')
    xlabel('Center Contrast'); ylabel('Perceived Contrast')
    %     set(gcf, 'Position',scrsz);
end

if groupPlots
    %     figure(3) % Group Location Error
    %     nbins = 100;
    %     histogram(mean(locationErrorMat0,3),nbins,'Normalization','count','BinWidth',5)
    %     if numel(subjectProfile.TheData{currSubj}) > 4
    %         hold on
    %         histogram(mean(locationErrorMat1,3),nbins,'Normalization','count','BinWidth',5) % plot second order
    %     end
    %     box off; set(gca, 'TickDir','out','ColorOrder',[1 0 0; 0 0 1]); xlim([-180 180]);
    %     title('Group Location Error');
    %     legend('Con-Loc','Loc-Con')
    
    figure % Group Contrast Estimates
    errorbar(plotContrasts,100.*squeeze(mean(subjectProfile.ContrastEstimate(:,:,1))),100.*std(subjectProfile.ContrastEstimate(:,:,1))./sqrt(length(subjects)),'.-','CapSize',0);
    hold on
    errorbar(plotContrasts,100.*squeeze(mean(subjectProfile.ContrastEstimate(:,:,2))),100.*std(subjectProfile.ContrastEstimate(:,:,2))./sqrt(length(subjects)),'.-','CapSize',0);
    errorbar(plotContrasts,100.*squeeze(mean(subjectProfile.ContrastEstimate(:,:,3))),100.*std(subjectProfile.ContrastEstimate(:,:,3))./sqrt(length(subjects)),'.--','CapSize',0);
    errorbar(plotContrasts,100.*squeeze(mean(subjectProfile.ContrastEstimate(:,:,4))),100.*std(subjectProfile.ContrastEstimate(:,:,4))./sqrt(length(subjects)),'.--','CapSize',0);
    line([plotContrasts(1) plotContrasts(end)],[plotContrasts(1) plotContrasts(end)],'Color','k')
    xticks(plotContrasts); yticks(plotContrasts);
    xticklabels({plotContrasts}); yticklabels({plotContrasts});
    set(gca,'XScale','log','YScale','log','ColorOrder',[1 0 0; 0 0 1; 0.5 0 0; 0 0 0.5],'TickDir','out','XTick',plotContrasts,'YTick',plotContrasts,...
        'XTickLabel',plotContrasts,'YTickLabel',plotContrasts,...
        'TickDir','out'); box off; xlim([min(plotContrasts) max(plotContrasts)]); ylim([min(plotContrasts) max(plotContrasts)]);
    legend('Sim.','Seq.','Sim. BL','Seq. BL','Location','NorthWest')
    title('Group: Presented v. Perceived Contrast')
    xlabel('Center Contrast'); ylabel('Perceived Contrast')
end

if superPlots
    figure % Group Contrast Estimates
    errorbar(plotContrasts,100.*mean(pooledContrastEstimatesMat(:,:,1)),100.*std(pooledContrastEstimatesMat(:,:,1))./sqrt(length(subjects)),'.-','CapSize',0);
    hold on
    errorbar(plotContrasts,100.*mean(pooledContrastEstimatesMat(:,:,2)),100.*std(pooledContrastEstimatesMat(:,:,2))./sqrt(length(subjects)),'.-','CapSize',0);
    errorbar(plotContrasts,100.*mean(pooledContrastEstimatesMat(:,:,3)),100.*std(pooledContrastEstimatesMat(:,:,3))./sqrt(length(subjects)),'.--','CapSize',0);
    errorbar(plotContrasts,100.*mean(pooledContrastEstimatesMat(:,:,4)),100.*std(pooledContrastEstimatesMat(:,:,4))./sqrt(length(subjects)),'.--','CapSize',0);
    line([plotContrasts(1) plotContrasts(end)],[plotContrasts(1) plotContrasts(end)],'Color','k')
    xticks(plotContrasts); yticks(plotContrasts);
    xticklabels({plotContrasts}); yticklabels({plotContrasts});
    set(gca,'XScale','log','YScale','log','ColorOrder',[1 0 0; 0 0 1; 0.5 0 0; 0 0 0.5],'TickDir','out','XTick',plotContrasts,'YTick',plotContrasts,...
        'XTickLabel',plotContrasts,'YTickLabel',plotContrasts,...
        'TickDir','out'); box off; xlim([min(plotContrasts) max(plotContrasts)]); ylim([min(plotContrasts) max(plotContrasts)]);
    legend('Sim.','Seq.','Sim. BL','Seq. BL','Location','NorthWest')
    title('Pooled: Presented v. Perceived Contrast')
    xlabel('Center Contrast'); ylabel('Perceived Contrast')
end

%% Contrast Estimate Distribution Fits
% Individual subject contrast estimate distrubution fits

plot_indx = 0:5:length(subjects)*length(centerContrast);
figure %Contrast Estimates + Fits
for ns = 1:numel(subjects)
    for nc = 1:length(centerContrast)
        subplot(length(subjects),length(centerContrast),plot_indx(ns)+nc)
        
        hold on
        for nf = 1:numel(dataFields)
            
            tmp_hist = histogram(subjectProfile.OrganizedData{ns}(nc).(dataFields{nf})(:,1),'Normalization','pdf');
            tmp_hist.BinEdges = logspace(log10(.01),log10(1),25);
            set(gca, 'XScale','log','XTick',0:0.1:1, 'XTickLabelRotation', -45)
            xlim([0 1])
            ydata = tmp_hist.Values;
            xdata = logspace(log10(.01),log10(1),49);
            startValues = [log10(centerContrast(nc)) 0.1 1];
            options = optimset('MaxFunEvals', 5000.*numel(startValues), 'MaxIter', 5000.*numel(startValues));
            
            [est_params_tmp(nf).params(ns,nc,:), est_params_tmp(nf).sse(ns,nc)] = fminsearch('mygauss', startValues, options, ydata, log10(xdata(2:2:end)));
            est_y = est_params_tmp(nf).params(ns,nc,3)*exp(-(log10(xdata)-est_params_tmp(nf).params(ns,nc,1)).^2/(2*(est_params_tmp(nf).params(ns,nc,2)^2)));
            plot(xdata,est_y)
            conEstDists.x(:,nc,nf,ns) = xdata;
            conEstDists.y(:,nc,nf,ns) = est_y;
        end
        line([centerContrast(nc) centerContrast(nc)],[0 max(ydata)],'Color','r')
        set(gca,'TickDir','out','ColorOrder',[1 0 0; 1 0 0; 0 0 1; 0 0 1; 0.5 0 0; 0.5 0 0; 0 0 0.5; 0 0 0.5]); box off; %title([num2str(plotContrasts(nc))])
    end
end
legend('Sim.', 'Sim. Fit', 'Seq','Seq. Fit', 'Sim. BL', 'Sim. BL Fit','Seq. BL','Seq. BL Fit')
xlabel('Presented contrast'); ylabel('PDF');
suptitle('Contrast Estimates + Fits (subject)')

figure % Contrast Estimate Fits
for ns = 1:numel(subjects)
    for nc = 1:length(centerContrast)
        subplot(length(subjects),length(centerContrast),plot_indx(ns)+nc)
        hold on
        for nf = 1:numel(dataFields)
            if nf <= 2
                plot(100.*conEstDists.x(:,nc,nf,ns),conEstDists.y(:,nc,nf,ns),'-')
            else
                plot(100.*conEstDists.x(:,nc,nf,ns),conEstDists.y(:,nc,nf,ns),'--')
            end
        end
        set(gca,'TickDir','out','XScale','log','XTick',0:10:100,'ColorOrder',[ 1 0 0; 0 0 1; 0.5 0 0; 0 0 0.5], 'XTickLabelRotation', -45); box off; %title([num2str(plotContrasts(nc))])
        line([plotContrasts(nc) plotContrasts(nc)],[0 max(max(squeeze(max(conEstDists.y(:,nc,:,ns)))))],'Color','k','LineStyle','-')
    end
end
legend('Sim. Fit','Seq. Fit', 'Sim. BL Fit','Seq. BL Fit')
xlabel('Presented contrast'); ylabel('PDF');
suptitle('Contrast Estimate Fits (subject)')

% Group Contrast Estimate Fits
conEstDists.groupx = mean(conEstDists.x,4);
conEstDists.groupy = mean(conEstDists.y,4);

figure 
for nc = 1:length(centerContrast)
    subplot(1,length(centerContrast),nc)
    hold on
    for nf = 1:numel(dataFields)
        if nf <= 2
            plot(100.*conEstDists.groupx(:,nc,nf),conEstDists.groupy(:,nc,nf),'-','LineWidth',1.5)
        else
            plot(100.*conEstDists.groupx(:,nc,nf),conEstDists.groupy(:,nc,nf),'--','LineWidth',1.5)
        end
    end
    set(gca,'TickDir','out','XScale','log','XTick',0:10:100,'ColorOrder',[1 0 0; 0 0 1; 0.5 0 0; 0 0 0.5], 'XTickLabelRotation', -45); box off; %title([num2str(plotContrasts(nc))])
    line([plotContrasts(nc) plotContrasts(nc)],[0 max(max(squeeze(max(conEstDists.groupy(:,nc,:)))))],'Color','k','LineStyle','-','LineWidth',1.5)
    title(num2str(plotContrasts(nc)))
end
legend('Sim. Fit','Seq. Fit', 'Sim. BL Fit','Seq. BL Fit')
xlabel('Presented contrast'); ylabel('PDF');
suptitle('Contrast Estimate Fits (group)')

% Contrast Estimates Summary Figure
subjcolor = num2cell(jet(numel(subjects)),2);
condColor = [{[1 0 0]}; {[0 0 1]}];

figure('color', [1 1 1])
subplot(1,3,1)
h = errorbar([1:3:15; 1.5:3:16]', [mean(10.^est_params_tmp(1).params(:,:,1)); mean(10.^est_params_tmp(2).params(:,:,1))]',  ...
    ([std(10.^est_params_tmp(1).params(:,:,1)); std(10.^est_params_tmp(2).params(:,:,1))]./sqrt(numel(subjects)))', 'ok', 'Capsize', 0, 'LineStyle', 'None','LineWidth',2,'MarkerSize',8);
set(h, {'color'}, condColor, {'MarkerFaceColor'}, condColor)
hold on

for ns = 1:length(subjects)
    plot([1:3:15; 1.5:3:16],[10.^est_params_tmp(1).params(ns,:,1); 10.^est_params_tmp(2).params(ns,:,1)],'LineStyle','-','Color',subjcolor{ns},'Marker','.','MarkerSize',8,'MarkerEdgeColor',subjcolor{ns})
end
%     plot(repmat([1:3:15],numel(subjects),1),10.^est_params_tmp(1).params(:,:,1),'LineStyle','none','Marker','.','MarkerSize',8,'MarkerEdgeColor','k')
%     plot(repmat([1.5:3:16],numel(subjects),1),10.^est_params_tmp(2).params(:,:,1),'LineStyle','none','Marker','.','MarkerSize',8,'MarkerEdgeColor','k')

xlim([0 16]); title('Mean estimates'); box off
set(gca, 'XTick', [1.25:3:15], 'XTickLabel', plotContrasts, 'TickDir', 'out')
xlabel('Presented contrast'); ylabel('Perceived contrast')

subplot(1,3,2)
h = errorbar([1:3:15; 1.5:3:16]', [mean(est_params_tmp(1).params(:,:,2)); mean(est_params_tmp(2).params(:,:,2))]',  ...
    ([std(est_params_tmp(1).params(:,:,2)); std(est_params_tmp(2).params(:,:,2))]./sqrt(numel(subjects)))', 'ok', 'Capsize', 0, 'LineStyle', 'None','LineWidth',2,'MarkerSize',8);
set(h, {'color'}, condColor, {'MarkerFaceColor'}, condColor)
hold on
for ns = 1:length(subjects)
    plot([1:3:15; 1.5:3:16],[est_params_tmp(1).params(ns,:,2); est_params_tmp(2).params(ns,:,2)],'LineStyle','-','Color',subjcolor{ns},'Marker','.','MarkerSize',8,'MarkerEdgeColor',subjcolor{ns})
end
%     plot(repmat([1:3:15],numel(subjects),1),est_params_tmp(1).params(:,:,2),'LineStyle','none','Marker','.','MarkerSize',8,'MarkerEdgeColor','k')
%     plot(repmat([1.5:3:16],numel(subjects),1),est_params_tmp(2).params(:,:,2),'LineStyle','none','Marker','.','MarkerSize',8,'MarkerEdgeColor','k')
xlim([0 16]); title('Width estimates'); box off
set(gca, 'XTick', [1.25:3:15], 'XTickLabel', plotContrasts, 'TickDir', 'out')
ylabel('Width')

subplot(1,3,3)
h = errorbar([1:3:15; 1.5:3:16]', [mean(est_params_tmp(1).sse(:,:)); mean(est_params_tmp(2).sse(:,:))]',  ...
    ([std(est_params_tmp(1).sse(:,:)); std(est_params_tmp(2).sse(:,:))]./sqrt(numel(subjects)))', 'ok', 'Capsize', 0, 'LineStyle', 'None','LineWidth',2,'MarkerSize',8);
set(h, {'color'}, condColor, {'MarkerFaceColor'}, condColor)
hold on
for ns = 1:length(subjects)
    plot([1:3:15; 1.5:3:16],[est_params_tmp(1).sse(ns,:); est_params_tmp(2).sse(ns,:)],'LineStyle','-','Color',subjcolor{ns},'Marker','.','MarkerSize',8,'MarkerEdgeColor',subjcolor{ns})
end
%     plot(repmat([1:3:15],numel(subjects),1),est_params_tmp(1).sse(:,:),'LineStyle','none','Marker','.','MarkerSize',8,'MarkerEdgeColor','k')
%     plot(repmat([1.5:3:16],numel(subjects),1),est_params_tmp(2).sse(:,:),'LineStyle','none','Marker','.','MarkerSize',8,'MarkerEdgeColor','k')
xlim([0 16]); title('SSE'); box off
set(gca, 'XTick', [1.25:3:15], 'XTickLabel', plotContrasts, 'TickDir', 'out')
ylabel('SSE')

% Pooled Contrast Estimates Distribution Fitting
for nc = 1:length(centerContrast)
    subplot(1,length(centerContrast),nc)
    
    hold on
    for nf = 1:numel(dataFields)
        
        tmp_hist = histogram(pooledContrastEstimatesMat(:,nc,nf),'Normalization','pdf');
        tmp_hist.BinEdges = logspace(log10(.01),log10(1),25);
        set(gca, 'XScale','log','XTick',0:0.1:1, 'XTickLabelRotation', -45)
        xlim([0 1])
        ydata = tmp_hist.Values;
        xdata = logspace(log10(.01),log10(1),49);
        startValues = [log10(centerContrast(nc)) 0.1 1];
        options = optimset('MaxFunEvals', 5000.*numel(startValues), 'MaxIter', 5000.*numel(startValues));
        
        [pooled_est_params_tmp(nf).params(nc,:), pooled_est_params_tmp(nf).sse(nc)] = fminsearch('mygauss', startValues, options, ydata, log10(xdata(2:2:end)));
        pooled_est_y = pooled_est_params_tmp(nf).params(nc,3)*exp(-(log10(xdata)-pooled_est_params_tmp(nf).params(nc,1)).^2/(2*(pooled_est_params_tmp(nf).params(nc,2)^2)));
        plot(xdata,pooled_est_y)
        pooled_conEstDists.x(:,nc,nf) = xdata;
        pooled_conEstDists.y(:,nc,nf) =pooled_est_y;
    end
    line([centerContrast(nc) centerContrast(nc)],[0 max(ydata)],'Color','r')
    set(gca,'TickDir','out','ColorOrder',[1 0 0; 1 0 0; 0 0 1; 0 0 1; 0.5 0 0; 0.5 0 0; 0 0 0.5; 0 0 0.5]); box off; %title([num2str(plotContrasts(nc))])
end
legend('Sim.', 'Sim. Fit', 'Seq','Seq. Fit', 'Sim. BL', 'Sim. BL Fit','Seq. BL','Seq. BL Fit')
xlabel('Presented contrast'); ylabel('PDF');
suptitle('Contrast Estimates + Fits (pooled)')

figure % Just the fits
for nc = 1:length(centerContrast)
    subplot(1,length(centerContrast),nc)
    hold on
    for nf = 1:numel(dataFields)
        if nf <= 2
            plot(100.*pooled_conEstDists.x(:,nc,nf),pooled_conEstDists.y(:,nc,nf),'-','LineWidth',1.5)
        else
            plot(100.*pooled_conEstDists.x(:,nc,nf),pooled_conEstDists.y(:,nc,nf),'--','LineWidth',1.5)
        end
    end
    set(gca,'TickDir','out','XScale','log','XTick',0:10:100,'ColorOrder',[1 0 0; 0 0 1; 0.5 0 0; 0 0 0.5], 'XTickLabelRotation', -45); box off; %title([num2str(plotContrasts(nc))])
    line([plotContrasts(nc) plotContrasts(nc)],[0 max(max(squeeze(max(pooled_conEstDists.y(:,nc,:)))))],'Color','k','LineStyle','-','LineWidth',1.5)
    title(num2str(plotContrasts(nc)))
end
legend('Sim. Fit','Seq. Fit', 'Sim. BL Fit','Seq. BL Fit')
xlabel('Presented contrast'); ylabel('PDF');
suptitle('Contrast Estimate Fits (pooled)')
%% Probe Effect Analysis
for n=1:1
    % % Take estimates within and outside a given range and compare the two
    % % distribution
    % % compare sigma and mean
    % probeMat = [allLocationError, contrastError, probeOffset];
    % offsets = -180:180;
    % offsets = offsets(:);
    % offsetRange  = 20;
    % oppOffsetRange = 180-offsetRange;
    % global numContrasts
    % numContrasts =  numel(centerContrast);
    % if superPlots
    %     for nOrder = 0:1
    %         for nCondition = 1:numel(dataFields)
    %             %             offsetLocationEffect = figure('Position',[50 50 scrsz(3) scrsz(3)/2],'name',[dataFields{nCondition} 'Probe Offset Location Effect Order ' num2str(nOrder)]);
    %             %             offsetContrastEffect = figure('Position',[50 50 scrsz(3) scrsz(3)/2],'name',[dataFields{nCondition} 'Probe Offset Contrast Effect Order ' num2str(nOrder)]);
    %             offsetMeanEstimates = figure('Position',[50 50 scrsz(3) scrsz(3)/2],'name',[dataFields{nCondition} 'Probe Offset Mean Estimates Order ' num2str(nOrder)]);
    %             offsetWidthEstimates = figure('Position',[50 50 scrsz(3) scrsz(3)/2],'name',[dataFields{nCondition} 'Probe Offset Width Estimates Order ' num2str(nOrder)]);
    %             for nContrast = 1:numel(centerContrast)
    %                 currContrast = centerContrast(nContrast);
    %                 currIndx = Contrasts == currContrast & Conditions == nCondition & Orders == nOrder;
    %                 tmpProbeOffset = probeOffset(currIndx);
    %                 tmpLocationError = allLocationError(currIndx);
    %                 tmpContrastEstimates= contrastEstimates(currIndx);
    %
    %                 currProbeOffset = unique(tmpProbeOffset);
    %                 insideRange = currProbeOffset; insideRange(insideRange > offsetRange | insideRange < -offsetRange) = [];
    %
    %                 tmpoutsideRange = currProbeOffset;
    %                 outsideRange = tmpoutsideRange(tmpoutsideRange < -oppOffsetRange | tmpoutsideRange > oppOffsetRange);
    %
    %                 inside = abs(tmpProbeOffset) <= offsetRange;
    %                 outside = abs(tmpProbeOffset) >= oppOffsetRange;
    %
    %                 trialsInside(nContrast) = sum(inside);
    %                 trialsOutside(nContrast) = sum(outside);
    %
    %                 currLocationError_inside = tmpLocationError(inside);
    %                 currLocationError_outside = tmpLocationError(outside);
    %
    %                 currContrastEstimates_inside= tmpContrastEstimates(inside);
    %                 currContrastEstimates_outside = tmpContrastEstimates(outside);
    %
    %                 edges_loc = -180:180;
    %                 [N_bl_loc_probe_in(nContrast,:),x_loc] = histcounts(currLocationError_inside,edges_loc);
    %                 N_bl_loc_probe_in(nContrast,:) = N_bl_loc_probe_in(nContrast,:)./max(N_bl_loc_probe_in(nContrast,:));
    %                 [N_bl_loc_probe_out(nContrast,:),x_loc] = histcounts(currLocationError_outside,edges_loc);
    %                 N_bl_loc_probe_out(nContrast,:) = N_bl_loc_probe_out(nContrast,:)./max(N_bl_loc_probe_out(nContrast,:));
    %                 xvalues_loc = (x_loc(1:end-1)+x_loc(2:end))/2;
    %
    %                 edges_con = 0:0.05:1;
    %
    %                 [N_bl_con_probe_in(nContrast,:),x_con_probe] = histcounts(currContrastEstimates_inside,edges_con);
    %                 [N_bl_con_probe_out(nContrast,:),x_con_probe] = histcounts(currContrastEstimates_outside,edges_con);
    %
    %                 tmp_max = max( [max(N_bl_con_probe_in(nContrast,:)), max(N_bl_con_probe_out(nContrast,:))]);
    %
    %                 N_bl_con_probe_in(nContrast,:) = N_bl_con_probe_in(nContrast,:)./tmp_max;
    %                 N_bl_con_probe_out(nContrast,:) = N_bl_con_probe_out(nContrast,:)./tmp_max;
    %                 xvalues_con = (x_con_probe(1:end-1)+x_con_probe(2:end))/2;
    %             end
    %
    %
    %             % fit all subjects data simultaneous
    %             % mean, amplitude, sigma
    %             startValues = [centerContrast' repmat(0.1, [1 numel(centerContrast)])];
    %             options = optimset('MaxFunEvals', 5000.*numel(startValues), 'MaxIter', 5000.*numel(startValues));
    %
    %             [est_params_tmp_in, r2_bl_con_probe_in] = fminsearch('mygauss_allContrasts', startValues, options, N_bl_con_probe_in, xvalues_con);
    %             tmp_in = reshape(est_params_tmp_in, [5 2])';
    %             est_params_bl_con_probe_in = [tmp_in(1,:); ones(1,5); tmp_in(2,:)];
    %
    %             [est_params_tmp_out, r2_blcon_probe_out] = fminsearch('mygauss_allContrasts', startValues, options, N_bl_con_probe_out, xvalues_con);
    %             tmp_out = reshape(est_params_tmp_out, [5 2])';
    %             est_params_bl_con_probe_out = [tmp_out(1,:); ones(1,5); tmp_out(2,:)];
    %
    %             x_fit = 0.01:0.01:1;
    %
    %             for nContrast = 1:numel(centerContrast)
    %                 % Location Error as a function of Probe Offset
    %                 %             figure(offsetLocationEffect)
    %                 %             subplot(1,5, nContrast)
    %                 %             bar(xvalues_loc,N_bl_loc_probe_in(nContrast,:),'FaceAlpha',.5);
    %                 %             hold on
    %                 %             bar(xvalues_loc,N_bl_loc_probe_out(nContrast,:),'FaceAlpha',.5);
    %                 %             title([num2str(round(100*centerContrast(nContrast))) '%'])
    %                 %             xlabel('Probe Offset'); ylabel('Location Error')
    %
    %                 % Contrast Error as a function of Probe Offset
    %                 y_est_bl_con_probe_in = est_params_bl_con_probe_in(2,nContrast)*exp(-(x_fit-est_params_bl_con_probe_in(1,nContrast)).^2/(2*(est_params_bl_con_probe_in(3,nContrast)^2)));
    %                 y_est_bl_con_probe_out = est_params_bl_con_probe_out(2,nContrast)*exp(-(x_fit-est_params_bl_con_probe_out(1,nContrast)).^2/(2*(est_params_bl_con_probe_out(3,nContrast)^2)));
    %
    %                 %                 figure(offsetContrastEffect)
    %                 %                 subplot(1,5,nContrast)
    %                 %                 bar(xvalues_con,N_bl_con_probe_in(nContrast,:),'FaceAlpha',.5);
    %                 %                 hold on
    %                 %                 bar(xvalues_con,N_bl_con_probe_out(nContrast,:),'FaceAlpha',.5);
    %                 %                 hold on
    %                 %                 insideFit = plot(x_fit, y_est_bl_con_probe_in, 'b', 'LineWidth', 2);
    %                 %                 outsideFit = plot(x_fit, y_est_bl_con_probe_out, 'r', 'LineWidth', 2);
    %                 %                 title([{[num2str(round(100*centerContrast(nContrast))) '%']}; {['#trials in=' num2str(trialsInside(nContrast)) ',#trials out=' num2str(trialsOutside(nContrast))]}])
    %                 %                 xlabel('Contrast Estimate'); ylabel('Frequency')
    %                 %                 insideFit.Color(4)=0.5;outsideFit.Color(4)=0.5;
    %
    %                 % plot mu estimates
    %                 figure(offsetMeanEstimates)
    %                 subplot(1,5,nContrast)
    %                 bar(1,est_params_tmp_in(nContrast),'FaceAlpha',.5)
    %                 hold all
    %                 bar(2, est_params_tmp_out(nContrast),'FaceAlpha',.5)
    %                 set(gca, 'XTick', [1 2])
    %                 set(gca, 'XTickLabel', {'Inside' 'Outside'})
    %                 title([{'Mean Estimates'}; {[num2str(round(100*centerContrast(nContrast))) '%']}; {['#trials in=' num2str(trialsInside(nContrast)) ',#trials out=' num2str(trialsOutside(nContrast))]}])
    %                 ylim([0 1])
    %                 ylabel('Contrast Estimate')
    %
    %                 % plot sigma estiamtes
    %                 figure(offsetWidthEstimates)
    %                 subplot(1,5,nContrast)
    %                 bar(1,est_params_tmp_in(nContrast+numel(centerContrast)),'FaceAlpha',.5)
    %                 hold on
    %                 bar(2,est_params_tmp_out(nContrast+numel(centerContrast)),'FaceAlpha',.5)
    %                 set(gca, 'XTick', [1 2])
    %                 set(gca, 'XTickLabel', {'Inside' 'Outside'})
    %                 title([{'Width Estimates'}; {[num2str(round(100*centerContrast(nContrast))) '%']}; {['#trials in=' num2str(trialsInside(nContrast)) ',#trials out=' num2str(trialsOutside(nContrast))]}])
    %                 ylim([0 1])
    %                 ylabel('Contrast Estimate')
    %
    %             end
    %         end
    %     end
    % end
end
%% Supression Index
% Perception: first page.
suppressionIndex.Collapsed(1:numel(subjectProfile.SubjectName),:,1) = ((subjectProfile.ContrastEstimate(:,:,1)) - (subjectProfile.ContrastEstimate(:,:,3))) ...
    ./((subjectProfile.ContrastEstimate(:,:,1)) + (subjectProfile.ContrastEstimate(:,:,3)));
suppressionIndex.Order0(1:numel(subjectProfile.SubjectName),:,1) = ((subjectProfile.ContrastEstimateOrder0(:,:,1)) - (subjectProfile.ContrastEstimateOrder0(:,:,3))) ...
    ./((subjectProfile.ContrastEstimateOrder0(:,:,1)) + (subjectProfile.ContrastEstimateOrder0(:,:,3)));
suppressionIndex.Order1(1:numel(subjectProfile.SubjectName),:,1) = ((subjectProfile.ContrastEstimateOrder1(:,:,1)) - (subjectProfile.ContrastEstimateOrder1(:,:,3))) ...
    ./((subjectProfile.ContrastEstimateOrder1(:,:,1)) + (subjectProfile.ContrastEstimateOrder1(:,:,3)));
suppressionIndex.Collapsed_PS(1:numel(subjectProfile.SubjectName),:,1) = ((ContrastEstimate_PS(:,:,1)) - (ContrastEstimate_PS(:,:,3))) ...
    ./((ContrastEstimate_PS(:,:,1)) + (ContrastEstimate_PS(:,:,3)));
%Working Memory: second page.
suppressionIndex.Collapsed(1:numel(subjectProfile.SubjectName),:,2) = ((subjectProfile.ContrastEstimate(:,:,2)) - (subjectProfile.ContrastEstimate(:,:,4))) ...
    ./((subjectProfile.ContrastEstimate(:,:,2)) + (subjectProfile.ContrastEstimate(:,:,4)));
suppressionIndex.Order0(1:numel(subjectProfile.SubjectName),:,2) = ((subjectProfile.ContrastEstimateOrder0(:,:,2)) - (subjectProfile.ContrastEstimateOrder0(:,:,4))) ...
    ./((subjectProfile.ContrastEstimateOrder0(:,:,2)) + (subjectProfile.ContrastEstimateOrder0(:,:,4)));
suppressionIndex.Order1(1:numel(subjectProfile.SubjectName),:,2) = ((subjectProfile.ContrastEstimateOrder1(:,:,2)) - (subjectProfile.ContrastEstimateOrder1(:,:,4))) ...
    ./((subjectProfile.ContrastEstimateOrder1(:,:,2)) + (subjectProfile.ContrastEstimateOrder1(:,:,4)));
suppressionIndex.Collapsed_PS(1:numel(subjectProfile.SubjectName),:,2) = ((ContrastEstimate_PS(:,:,2)) - (ContrastEstimate_PS(:,:,4))) ...
    ./((ContrastEstimate_PS(:,:,2)) + (ContrastEstimate_PS(:,:,4)));

% T-test
% Suppression Index
% Two sample t-test, if h = 1 then the two populations come from unequal
% means. 5% significance level used (default). Tests the means for each
% contrast estimation (collapsed over subjects).

% Perciption and working memory
[hSI,pSI,ciSI,statsSI] = ttest2(mean(suppressionIndex.Collapsed(:,:,1),1),...
    mean(suppressionIndex.Collapsed(:,:,2),1));

%Perception and [0 0 0 0 0] suppression
[hSI_P,pSI_P,ciSI_P,statsSI_P] = ttest2(mean(suppressionIndex.Collapsed(:,:,1),1),[ 0 0 0 0 0]);

%Working Memory and [0 0 0 0 0] suppression
[hSI_WM,pSI_WM,ciSI_WM,statsSI_WM] = ttest2([0 0 0 0 0],mean(suppressionIndex.Collapsed(:,:,2),1));

%Suppression Index Plot %
%Collapsed over both orders
hold off
figure('Color', [1 1 1]);
set(gcf, 'Name', sprintf('Suppression Index: Simultaenous vs. Sequential'));
hold all;
%Perception
errorbar(plotContrasts', mean(suppressionIndex.Collapsed(:,:,1),1), std(suppressionIndex.Collapsed(:,:,1),1)...
    /sqrt(length(subjects)), 'r','LineWidth',3);
%Working Memory
errorbar(plotContrasts', mean(suppressionIndex.Collapsed(:,:,2),1), std(suppressionIndex.Collapsed(:,:,2),1)...
    /sqrt(length(subjects)), 'b','LineWidth',3);
plot(repmat(plotContrasts', [length(subjects) 1]), (suppressionIndex.Collapsed(:,:,1)),'r.','MarkerSize',10);
plot(repmat(plotContrasts', [length(subjects) 1]), (suppressionIndex.Collapsed(:,:,2)),'b.','MarkerSize',10);
hold on
plot(100*[0.09 0.8],[0 0],':k')
ylabel('Suppression Index (surround-nosurround)/(surround+nosurround)');
xlabel('Contrast (%)');
xlim(100*[0.09 0.8]);
ylim([-0.5 0.5]);
xticks(plotContrasts)
xticklabels(plotContrasts)
axis square;
legend({'Sim. Condition','Seq. Condition'})
set(gca,'XScale','log','TickDir','out')
title('Suppression Index')

%     % Split between orders%
%     figure('Color', [1 1 1]);
%
%     subplot(1,2,1)
%     set(gcf, 'Name', sprintf('Suppression Index: Perception vs. Working Memory'));
%     hold all;
%     %Perception
%     errorbar(centerContrast', nanmean(suppressionIndex.Order0(:,:,1),1), nanstd(suppressionIndex.Order0(:,:,1),1)...
%         /sqrt(length(subjects)), 'r','LineWidth',3)
%     %Working Memory
%     errorbar(centerContrast', nanmean(suppressionIndex.Order0(:,:,2),1), nanstd(suppressionIndex.Order0(:,:,2),1)...
%         /sqrt(length(subjects)), 'b','LineWidth',3); hold all;
%     plot(repmat(centerContrast', [length(subjects) 1]), (suppressionIndex.Order0(:,:,1)),'r.','MarkerSize',10);
%     plot(repmat(centerContrast', [length(subjects) 1]), (suppressionIndex.Order0(:,:,2)),'b.','MarkerSize',10);
%     hold all
%     plot([0.09 0.8],[0 0],':k')
%     ylabel('Suppression Index (surround-nosurround)/(surround+nosurround)');
%     xlabel('Contrast (%)');
%     xlim([0.09 0.8]);
%     ylim([-0.3 0.3]);
%     title('Order 0: Contrast, then Location');
%     axis square;
%     legend({'Perception Condition','Working Memory Condition'})
%
%     subplot(1,2,2)
%     hold all;
%     %Perception
%     errorbar(centerContrast', mean(suppressionIndex.Order1(:,:,1),1), std(suppressionIndex.Order1(:,:,1),1)...
%         /sqrt(size(visualmemory_subjectsRan,2)), 'r','LineWidth',3);
%     %Working Memory
%     errorbar(centerContrast', mean(suppressionIndex.Order1(:,:,2),1), std(suppressionIndex.Order1(:,:,2),1)...
%         /sqrt(length(subjects)), 'b','LineWidth',3);
%     plot(repmat(centerContrast', [length(subjects) 1]), (suppressionIndex.Order1(:,:,1)),'r.','MarkerSize',10);
%     plot(repmat(centerContrast', [length(subjects) 1]), (suppressionIndex.Order1(:,:,2)),'b.','MarkerSize',10);
%     hold all
%     plot([0.09 0.8],[0 0],':k')
%     ylabel('Suppression Index (surround-nosurround)/(surround+nosurround)');
%     xlabel('Contrast (%)');
%     xlim([0.09 0.8]);
%     ylim([-0.5 0.5]);
%     title('Order 1: Location, then Contrast');
%     axis square;
%     legend({'Perception Condition','Working Memory Condition'})

%% Normalization Model Fitting
experiments = {'Sim. Condition','Seq. Memory'};
%Model Fit Setup
surroundContrast = theData(1).p.surroundContrast; %Surround grating has 100% contrast level, =1
C50 = .6;
Wi = 0;
Wi_var = 0;
b = 0;
n = 2;
C_Test = centerContrast;
for e = 1:numel(experiments)
    for subjCount = 1:numel(subjectProfile.SubjectName)
        % Baselines are NOT currently separated. These are the baseline
        % estimates for working memory and perception combined.
        baselineMat(subjCount,:,1) = subjectProfile.ContrastEstimate(subjCount,:,3);
        baselineMat(subjCount,:,2) = subjectProfile.ContrastEstimate(subjCount,:,4);
        
        if e == 1 %Perception Experiment %
            variableMat(subjCount,:) = subjectProfile.ContrastEstimate(subjCount,:,1);
        elseif e == 2  % Working Memory Experiment %
            variableMat(subjCount,:) = subjectProfile.ContrastEstimate(subjCount,:,2);
        end
        
        % Fit individual data with Normalization model
        options = optimset('MaxFunEvals', 10^6, 'MaxIter', 10^6);
        % Use formula from Xing&Heeger 2001 to use baseline data to constrain alpha and n, and look for the inhibitory weight
        % to explain the suppressive influence of the surround
        % Wi: suppression weight --> surround induced
        % normalization
        % C50: inflection point
        % n: nonlinear transducer, determining steepness
        
        indv_r2 = zeros(2,numel(subjectProfile.SubjectName),2);
        startValuesModelFitting = [C50 n Wi_var];
        Data = {centerContrast, baselineMat(subjCount,:,e), variableMat(subjCount,:), surroundContrast};
        [est_params(e,subjCount,:), r2(e,subjCount)] = fminsearch('fitNormalizationModel_contrastMatch', startValuesModelFitting, options, Data);
        
        
        % Fit Y data based off of estimated parameters
        Y_base = (centerContrast.^est_params(e,subjCount,2)) ./ ...
            ((est_params(e,subjCount,1).^est_params(e,subjCount,2)) + (centerContrast.^est_params(e,subjCount,2)));
        
        Y_var = (centerContrast.^est_params(e,subjCount,2)) ./ ...
            ((est_params(e,subjCount,1).^est_params(e,subjCount,2)) + (centerContrast.^est_params(e,subjCount,2)) + (est_params(e,subjCount,3))*(surroundContrast.^est_params(e,subjCount,2)));
        
        % R^2 values for baseline and variable condition
        r2_BL = 1 - (sum((baselineMat(subjCount,:,e) - Y_base').^2) / sum((baselineMat(subjCount,:,e) - mean(baselineMat(subjCount,:,e))).^2));
        r2_V = 1 - (sum((variableMat(subjCount,:) - Y_var').^2) / sum((variableMat(subjCount,:) - mean(variableMat(subjCount,:))).^2));
        ind_r2(e,subjCount,:) = [r2_BL r2_V]';
        
        if normalizationModelFitPlots
            %Difference in perceived contrast compared to center contrast
            
            C_fit = 0.1:0.01:0.8;
            Y_var = (C_fit.^est_params(e,subjCount,2)) ./ ...
                ((est_params(e,subjCount,1).^est_params(e,subjCount,2)) + (C_fit.^est_params(e,subjCount,2)) + (est_params(e,subjCount,3)*(surroundContrast.^est_params(e,subjCount,2))));
            
            %Display Fits
            subplot(2,round(numel(subjectProfile.SubjectName)/2), subjCount)
            hold all
            colormap lines
            if e == 1
                loglog(centerContrast,variableMat(subjCount,:),'Linewidth',1.25)
                loglog(centerContrast,baselineMat(subjCount,:,e),'Linewidth',1.25)
            elseif e == 2
                loglog(centerContrast,variableMat(subjCount,:),'Linewidth',1.25)
                loglog(centerContrast,baselineMat(subjCount,:,e),'Linewidth',1.25)
            end
            loglog(C_Test,Y_base,'Linewidth',2)
            loglog(C_fit, Y_var,'Linewidth',2)
            hold all;
            xlim([0.1 0.8]), ylim([0.1 0.8]); box off
            ylabel('Perceived Contrast')
            xlabel('Center Contrast');
            plot([0.09 0.8], [0.09, 0.8], 'k:');
            if e == 1
                legend({'Perception Est.','Baseline Est.','Baseline Fit','Perception Fit'},'Location','northwest')
            elseif e == 2
                legend({'Working Memory Est.','Baseline Est.','Baseline Fit','Working Mem Fit'},'Location','northwest')
            end
            title(['Subject ' num2str(subjCount)]);
            axis square;
            hold all;
            set(gca,'YScale','log','XScale','log')
        end
    end
end

figure('color', [1 1 1])
subplot(1,4,1)
errorbar(1:2, [mean(est_params(1,:,1)) mean(est_params(2,:,1))], [std(est_params(1,:,1)) std(est_params(2,:,1))]/sqrt(size(est_params,2)), 'ok', ...
    'MarkerSize', 10, 'MarkerFaceColor', [0 0 0], 'CapSize', 0, 'LineWidth', 3)
hold all,
h = scatter(ones(size(est_params,2),1), est_params(1,:,1), 100, jet(size(est_params,2)), 'filled');
scatter(2*ones(size(est_params,2),1), est_params(2,:,1), 100, jet(size(est_params,2)), 'filled')
title('C50 est'), xlim([0 3]), ylim([0 1]); set(gca, 'XTick', 1:2, 'XTickLabel', {'Sim.' 'Seq.'}, 'XTickLabelRotation', -45,'TickDir','out')
box  off;

subplot(1,4,2)
errorbar(1:2, [mean(est_params(1,:,2)) mean(est_params(2,:,2))], [std(est_params(1,:,2)) std(est_params(2,:,2))]/sqrt(size(est_params,2)), 'ok', ...
    'MarkerSize', 10, 'MarkerFaceColor', [0 0 0], 'CapSize', 0, 'LineWidth', 3)
hold all,
h = scatter(ones(size(est_params,2),1), est_params(1,:,2), 100, jet(size(est_params,2)), 'filled');
scatter(2*ones(size(est_params,2),1), est_params(2,:,2), 100, jet(size(est_params,2)), 'filled')
title('Slope est'), xlim([0 3]), ylim([0 2]); set(gca, 'XTick', 1:2, 'XTickLabel', {'Sim.' 'Seq.'}, 'XTickLabelRotation', -45,'TickDir','out')
box  off;

subplot(1,4,3)
errorbar(1:2, [mean(est_params(1,:,3)) mean(est_params(2,:,3))], [std(est_params(1,:,3)) std(est_params(2,:,3))]/sqrt(size(est_params,2)), 'ok', ...
    'MarkerSize', 10, 'MarkerFaceColor', [0 0 0], 'CapSize', 0, 'LineWidth', 3)
hold all,
h = scatter(ones(size(est_params,2),1), est_params(1,:,3), 100, jet(size(est_params,2)), 'filled');
scatter(2*ones(size(est_params,2),1), est_params(2,:,3), 100, jet(size(est_params,2)), 'filled')
plot([0 3], [0 0], 'k:')
title('Wi est'), xlim([0 3]), ylim([-0.3 0.5]); set(gca, 'XTick', 1:2, 'XTickLabel', {'Sim.' 'Seq.'}, 'XTickLabelRotation', -45,'TickDir','out')
box  off;

subplot(1,4,4)
errorbar(1:2, [mean(ind_r2(1,:,2)) mean(ind_r2(2,:,2))], [std(ind_r2(1,:,2)) std(ind_r2(2,:,2))]/sqrt(size(ind_r2,2)), 'ok', ...
    'MarkerSize', 10, 'MarkerFaceColor', [0 0 0], 'CapSize', 0, 'LineWidth', 3)
hold all,
h = scatter(ones(size(ind_r2,2),1), ind_r2(1,:,2), 100, jet(size(ind_r2,2)), 'filled');
scatter(2*ones(size(ind_r2,2),1), ind_r2(2,:,2), 100, jet(size(ind_r2,2)), 'filled')
plot([0 3], [0 0], 'k:')
title('R2'), xlim([0 3]), ylim([0.5 1]); set(gca, 'XTick', 1:2, 'XTickLabel', {'Sim.' 'Seq.'}, 'XTickLabelRotation', -45,'TickDir','out')
box  off;

suptitle('Normalization Model Fitting')
% T-test
% Wi est - surround induced normalization parameter
% two sampled t est, testing for significance between perception and
% working memory conditions based off of hypothesis of an equal mean
[h_NP,p_NP,ci_NP,stats_NP] = ttest2(est_params(1,:,3),est_params(2,:,3));


%% Super Subject Plots
% if superPlots
%     %         % Total Location Error for baseline condition
%     %         for nCondition = 1:numel(dataFields)
%     %             figure('Position',[50 50 scrsz(3) scrsz(3)/2], 'color', [1 1 1], 'name', [dataFields{nCondition} ' Location Error'])
%     %             for nContrast = 1:numel(centerContrast)
%     %                 currInd = Conditions == nCondition & Contrasts == centerContrast(nContrast);
%     %                 edges_loc = -180:180;
%     %                 [N_bl_loc(nContrast,:),x_loc] = histcounts(allLocationError(currInd),edges_loc);
%     %                 N_bl_loc(nContrast,:) = N_bl_loc(nContrast,:)./max(N_bl_loc(nContrast,:));
%     %                 xvalues_loc = (x_loc(1:end-1)+x_loc(2:end))/2;
%     %
%     %                 subplot(1,5,nContrast)
%     %                 bar(xvalues_loc, N_bl_loc(nContrast,:));
%     %                 hold on
%     %                 title([num2str(round(100*centerContrast(nContrast))) '%'])
%     %                 xlabel('Location Error'); ylabel('Frequency')
%     %             end
%     %         end
%     %         %Location Error histogram split up by order, for baseline condition
%     %         figure
%     %         nbins = 100;
%     %         histfit(allLocationError(Conditions == 3 & Orders == 1),nbins)
%     %         hold on
%     %         histfit(allLocationError(Conditions == 3 & Orders == 0),nbins)
%     %         title("'Super Subject' Location Error by Order")
%     %         legend({'Loc-Con','Con-Loc'})
%     %         xlabel('Location Error'); ylabel('Frequency')
%     
%     baselineLocationErr1 = mean(abs(allLocationError(Conditions == 3)));
%     baselineLocationErr_STE1 = std(abs(allLocationError(Conditions == 3)))/sqrt(length(allLocationError(Conditions == 3)));
%     wmLocationErr1 = mean(abs(allLocationError(Conditions == 2)));
%     wmLocationErr_STE1 = std(abs(allLocationError(Conditions == 2)))/sqrt(length(allLocationError(Conditions == 2)));
%     percLocationErr1 = mean(abs(allLocationError(Conditions == 1)));
%     percLocationErr_STE1 = std(abs(allLocationError(Conditions == 1)))/sqrt(length(allLocationError(Conditions == 1)));
%     
%     for nCon = 1:length(centerContrast)
%         baselineLocationErr2(nCon,:) = mean(abs(allLocationError(Conditions == 3 & Contrasts == centerContrast(nCon))));
%         baselineLocationErr_STE2(nCon,:)  = std(abs(allLocationError(Conditions == 3 & Contrasts == centerContrast(nCon))))/...
%             sqrt(length(allLocationError(Conditions == 3 & Contrasts == centerContrast(nCon))));
%         wmLocationErr2(nCon,:)  = mean(abs(allLocationError(Conditions == 2 & Contrasts == centerContrast(nCon))));
%         wmLocationErr_STE2(nCon,:)  = std(abs(allLocationError(Conditions == 2 & Contrasts == centerContrast(nCon))))/...
%             sqrt(length(allLocationError(Conditions == 2 & Contrasts == centerContrast(nCon))));
%         percLocationErr2(nCon,:)  = mean(abs(allLocationError(Conditions == 1 & Contrasts == centerContrast(nCon))));
%         percLocationErr_STE2(nCon,:)  = std(abs(allLocationError(Conditions == 1 & Contrasts == centerContrast(nCon))))/...
%             sqrt(length(allLocationError(Conditions == 1)));
%     end
%     
%     figure('name','Super Avg. Loc. Err. by condition')
%     bar([baselineLocationErr1 percLocationErr1 wmLocationErr1])
%     set(gca,'TickDir','out')
%     hold on
%     errorbar(1:3,[baselineLocationErr1 percLocationErr1 wmLocationErr1], [baselineLocationErr_STE1 percLocationErr_STE1 wmLocationErr_STE1],'k','LineStyle','None' )
%     xticklabels({'Baseline' 'Perception' 'WM'})
%     box off
%     
%     figure('name','Super Avg. Loc. Err. by contrast & condition')
%     y = [baselineLocationErr2 percLocationErr2 wmLocationErr2];
%     err = [baselineLocationErr_STE2 percLocationErr_STE2 wmLocationErr_STE2];
%     bar(y)
%     hold on
%     ngroups = size(y, 1);
%     nbars = size(y, 2);
%     % Calculating the width for each bar group to plot error
%     groupwidth = min(0.8, nbars/(nbars + 1.5));
%     for i = 1:nbars
%         x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%         errorbar(x, y(:,i), err(:,i), 'k','LineStyle','none');
%     end
%     xlabel('Target Contrasts (%)'); ylabel('Location Error (deg)')
%     xticks(1:5);xticklabels({'10%' '17%' '27%' '45%' '75%'})
%     set(gca,'TickDir','out'); box off;
%     legend({'Baseline' 'Perception' 'WM'})
% end