%%% vmmf_subject_analysis.m
%% Prepare data and working directories
clear all; close all; clc;
scrsz = get(groot,'ScreenSize');

subjPlots = 0;
groupPlots = 0;
superPlots = 0;
normalizationModelPlots = 1;

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
locationErrorMatSuper = [];
conditionMatSuper = [];
locationAveragesByCond = zeros(3,length(subjects));

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
    
    % Save out average location error Per Condition Per Subject
    for i = 1:3
        locationAveragesByCond(i,currSubj) = mean(mean(abs(locationErrorMat(conditionMat==i))));
    end
    
    % Save out organized data structure
    subjectProfile.OrganizedData{currSubj} = ContrastData;
 
    % Save out mean location error
    subjectProfile.LocationError(currSubj) = mean(mean(abs(locationErrorMat)));
   
  
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

%% Probe Effect Analysis
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
%% Distribution Analysis
% % Fit Gaussian distribution to all participant data
% % Y=Amplitude*exp(-0.5*((X-Mean)/SD)^2)
% 
% for nCondition = 1:numel(dataFields)
%     for nContrast = 1:numel(centerContrast)
%         currInd = Conditions == nCondition & Contrasts == centerContrast(nContrast);
%         currData = contrastEstimates(currInd);
%         edges_con = 0:0.05:1;
%         [N_bl_con(nContrast,:),x_con] = histcounts(currData,edges_con);
%         N_bl_con(nContrast,:) = N_bl_con(nContrast,:)./max(N_bl_con(nContrast,:));
%         xvalues_con = (x_con(1:end-1)+x_con(2:end))/2;
%     end
%     
%     % fit all subjects data simultaneous
    % mean, amplitude, sigma
    startValues = [centerContrast' repmat(0.1, [1 numel(centerContrast)])];
   % options = optimset('MaxFunEvals', 5000.*numel(startValues), 'MaxIter', 5000.*numel(startValues));
%     
%     [est_params_tmp, r2_bl] = fminsearch('mygauss_allContrasts', startValues, options, N_bl_con, xvalues_con);
%     tmp = reshape(est_params_tmp, [5 2])';
%     est_params_bl = [tmp(1,:); ones(1,5); tmp(2,:)];
%     
%     % Generate gaussian fits from estimated params
%     
%     figure('Position',[50 50 scrsz(3) scrsz(3)/2], 'color', [1 1 1], 'name', [dataFields{nCondition} ' Contrast Estimates'])
%     for nContrast = 1:numel(centerContrast)
%         x_fit = 0.01:0.01:1;
%         y_est_bl = est_params_bl(2,nContrast)*exp(-(x_fit-est_params_bl(1,nContrast)).^2/(2*(est_params_bl(3,nContrast)^2)));
%         xvalues_con = (x_con(1:end-1)+x_con(2:end))/2;
%         subplot(1,5,nContrast)
%         bar(xvalues_con, N_bl_con(nContrast,:));
%         hold on
%         plot(x_fit, y_est_bl, 'r', 'LineWidth', 2)
%         box off; title([num2str(round(100*centerContrast(nContrast))) '%']); xlabel('Contrast (%)'); ylabel('Normalized count of perceived contrast');
%         axis square
%         %plot parameters alone
%     end
% end
% 
% %% Contrast Model Fitting
experiments = {'Perception Condition','Visual Working Memory'};
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
        baselineMat(subjCount,:) = subjectProfile.ContrastEstimate(subjCount,:,baselineIndex);
        if e == 1 %Perception Experiment %
            variableMat(subjCount,:) = subjectProfile.ContrastEstimate(subjCount,:,perceptionIndex);
        elseif e == 2  % Working Memory Experiment %
            variableMat(subjCount,:) = subjectProfile.ContrastEstimate(subjCount,:,workingMemoryIndex);
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
            Data = {centerContrast, baselineMat(subjCount,:), variableMat(subjCount,:), surroundContrast};
            [est_params(e,subjCount,:), r2(e,subjCount)] = fminsearch('fitNormalizationModel_contrastMatch', startValuesModelFitting, options, Data);


            % Fit Y data based off of estimated parameters
            Y_base = (centerContrast.^est_params(e,subjCount,2)) ./ ...
                ((est_params(e,subjCount,1).^est_params(e,subjCount,2)) + (centerContrast.^est_params(e,subjCount,2)));

            Y_var = (centerContrast.^est_params(e,subjCount,2)) ./ ...
                ((est_params(e,subjCount,1).^est_params(e,subjCount,2)) + (centerContrast.^est_params(e,subjCount,2)) + (est_params(e,subjCount,3))*(surroundContrast.^est_params(e,subjCount,2)));

            % R^2 values for baseline and variable condition
            r2_BL = 1 - (sum((baselineMat(subjCount,:) - Y_base).^2) / sum((baselineMat(subjCount,:) - mean(baselineMat(subjCount,:))).^2));
            r2_V = 1 - (sum((variableMat(subjCount,:) - Y_var).^2) / sum((variableMat(subjCount,:) - mean(variableMat(subjCount,:))).^2));
            indvR2(e,subjCount,:) = [r2_BL r2_V];

            if normalizationModelPlots
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
                    loglog(centerContrast,baselineMat(subjCount,:),'Linewidth',1.25)
                elseif e == 2
                    loglog(centerContrast,variableMat(subjCount,:),'Linewidth',1.25)
                    loglog(centerContrast,baselineMat(subjCount,:),'Linewidth',1.25)
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
suppressionIndex.Collapsed(1:numel(subjectProfile.SubjectName),:,2) = ((subjectProfile.ContrastEstimate(:,:,2)) - (subjectProfile.ContrastEstimate(:,:,3))) ...
    ./((subjectProfile.ContrastEstimate(:,:,2)) + (subjectProfile.ContrastEstimate(:,:,3)));
suppressionIndex.Order0(1:numel(subjectProfile.SubjectName),:,2) = ((subjectProfile.ContrastEstimateOrder0(:,:,2)) - (subjectProfile.ContrastEstimateOrder0(:,:,3))) ...
    ./((subjectProfile.ContrastEstimateOrder0(:,:,2)) + (subjectProfile.ContrastEstimateOrder0(:,:,3)));
suppressionIndex.Order1(1:numel(subjectProfile.SubjectName),:,2) = ((subjectProfile.ContrastEstimateOrder1(:,:,2)) - (subjectProfile.ContrastEstimateOrder1(:,:,3))) ...
    ./((subjectProfile.ContrastEstimateOrder1(:,:,2)) + (subjectProfile.ContrastEstimateOrder1(:,:,3)));
suppressionIndex.Collapsed_PS(1:numel(subjectProfile.SubjectName),:,2) = ((ContrastEstimate_PS(:,:,2)) - (ContrastEstimate_PS(:,:,3))) ...
    ./((ContrastEstimate_PS(:,:,2)) + (ContrastEstimate_PS(:,:,3)));


%% T-test 

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

% Wi est - surround induced normalization parameter
% two sampled t est, testing for significance between perception and
% working memory conditions based off of hypothesis of an equal mean
[h_NP,p_NP,ci_NP,stats_NP] = ttest2(est_params(1,:,3),est_params(2,:,3));




%% Super Subject Plots
if superPlots
%         % Total Location Error for baseline condition
%         for nCondition = 1:numel(dataFields)
%             figure('Position',[50 50 scrsz(3) scrsz(3)/2], 'color', [1 1 1], 'name', [dataFields{nCondition} ' Location Error'])
%             for nContrast = 1:numel(centerContrast)
%                 currInd = Conditions == nCondition & Contrasts == centerContrast(nContrast);
%                 edges_loc = -180:180;
%                 [N_bl_loc(nContrast,:),x_loc] = histcounts(allLocationError(currInd),edges_loc);
%                 N_bl_loc(nContrast,:) = N_bl_loc(nContrast,:)./max(N_bl_loc(nContrast,:));
%                 xvalues_loc = (x_loc(1:end-1)+x_loc(2:end))/2;
%     
%                 subplot(1,5,nContrast)
%                 bar(xvalues_loc, N_bl_loc(nContrast,:));
%                 hold on
%                 title([num2str(round(100*centerContrast(nContrast))) '%'])
%                 xlabel('Location Error'); ylabel('Frequency')
%             end
%         end
%         %Location Error histogram split up by order, for baseline condition
%         figure
%         nbins = 100;
%         histfit(allLocationError(Conditions == 3 & Orders == 1),nbins)
%         hold on
%         histfit(allLocationError(Conditions == 3 & Orders == 0),nbins)
%         title("'Super Subject' Location Error by Order")
%         legend({'Loc-Con','Con-Loc'})
%         xlabel('Location Error'); ylabel('Frequency')

    baselineLocationErr = mean(abs(allLocationError(Conditions == 3))); 
    baselineLocationErr_STE = std(abs(allLocationError(Conditions == 3)))/sqrt(length(allLocationError(Conditions == 3)));
    wmLocationErr = mean(abs(allLocationError(Conditions == 2)));
    wmLocationErr_STE = std(abs(allLocationError(Conditions == 2)))/sqrt(length(allLocationError(Conditions == 2)));
    percLocationErr = mean(abs(allLocationError(Conditions == 1)));
    percLocationErr_STE = std(abs(allLocationError(Conditions == 1)))/sqrt(length(allLocationError(Conditions == 1)));
    
    figure
    bar([baselineLocationErr percLocationErr wmLocationErr])
    set(gca,'TickDir','out')
    hold on
    errorbar(1:3,[baselineLocationErr percLocationErr wmLocationErr], [baselineLocationErr_STE percLocationErr_STE wmLocationErr_STE],'k','LineStyle','None' )
    xticklabels({'Baseline' 'Perception' 'WM'})
    box off
end

%% Group Plots
if groupPlots
    % Total Contrast Estimates
    plotContrasts = 100*round(centerContrast,2);
    figure
    cvp_group = loglog(plotContrasts, 100.*squeeze(mean(subjectProfile.ContrastEstimate)));
    hold on
    errorbar(plotContrasts,100.*mean(subjectProfile.ContrastEstimate(:,:,1)),100.*(std(subjectProfile.ContrastEstimate(:,:,1))/sqrt(size(subjectProfile.ContrastEstimate(:,:,1),1))),'b','LineStyle','None') % perception
    errorbar(plotContrasts,100.*mean(subjectProfile.ContrastEstimate(:,:,2)),100.*(std(subjectProfile.ContrastEstimate(:,:,2))/sqrt(size(subjectProfile.ContrastEstimate(:,:,2),1))),'r','LineStyle','None') % wm
    errorbar(plotContrasts,100.*mean(subjectProfile.ContrastEstimate(:,:,3)),100.*(std(subjectProfile.ContrastEstimate(:,:,3))/sqrt(size(subjectProfile.ContrastEstimate(:,:,3),1))),'k','LineStyle','None') % baseline
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
    box off
    hold off
    
    % Location Error by Condition
    figure
    bar(mean(locationAveragesByCond,2))
    set(gca,'TickDir','out')
    hold on
    errorbar(1:3,mean(locationAveragesByCond,2), std(locationAveragesByCond,[],2)/sqrt(size(locationAveragesByCond,2)),'k','LineStyle','None' )
    xticklabels({'Perception' 'WM' 'Baseline'})
    box off
    
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
    
    % Wi est - both perception and working memory conditions --> Surround induced normalization %
    normparams = mean(est_params(:,:,3),2);
        % with ALL data
        bar([0 1.0], mean(est_params(:,:,3),2), 0.3); % page 3 is surround induced noramlization parameter
        hold all
        handles = get(gca, 'Children');
        set(handles(1), 'FaceColor', [0 0 1], 'EdgeColor', 'none'); 
        errorbar([0 1], mean(est_params(:,:,3),2)',std(est_params(:,:,3),[], 2)'/sqrt(subjCount), 'k.')
        plot(repmat([0 1.0], [length(subjects) 1]), est_params(:,:,3)', 'ok');
        box off; xlim([-0.5 1.5]); set(gca, 'Xtick', [0 1], 'XtickLabel', {'Perc' 'vWM'})
        title(sprintf('Surround Induced Normalization\n Parameter'))
        
           % with ALL data
        bar([0 1.0], mean(est_params(:,:,3),2), 0.3); % page 3 is surround induced noramlization parameter
        hold all
        handles = get(gca, 'Children');
        set(handles(1), 'FaceColor', [0 0 1], 'EdgeColor', 'none'); 
        errorbar([0 1], mean(est_params(:,:,3),2)',std(est_params(:,:,3),[], 2)'/sqrt(subjCount), 'k.')
        plot(repmat([0 1.0], [length(subjects) 1]), est_params(:,:,3)', 'ok');
        box off; xlim([-0.5 1.5]); set(gca, 'Xtick', [0 1], 'XtickLabel', {'Perc' 'vWM'})
        title(sprintf('Surround Induced Normalization\n Parameter'))
        
%         % excluding subject 5 - temporarily
%         allbutvec = [1:4 6:10];
%         for i = 1:length(allbutvec)
%             est_params2(:,i,:) = est_params(:,i,:);
%         end
%         bar([0 1.0], mean(est_params2(:,:,3),2), 0.3); % page 3 is surround induced noramlization parameter
%         hold all
%         handles = get(gca, 'Children');
%         set(handles(1), 'FaceColor', [0 0 1], 'EdgeColor', 'none'); 
%         errorbar([0 1], mean(est_params2(:,:,3),2)',std(est_params2(:,:,3),[], 2)'/sqrt(subjCount), 'k.')
%         plot(repmat([0 1.0], [length(subjects) 1]), est_params2(:,:,3)', 'ok');
%         box off; xlim([-0.5 1.5]); set(gca, 'Xtick', [0 1], 'XtickLabel', {'Perc' 'vWM'})
%         title(sprintf('Surround Induced Normalization\n Parameter - No 5'))
%         ylim([-0.5 0.5])

    %Suppression Index%
    %Collapsed over both orders
    hold off
    figure('Color', [1 1 1]);
    set(gcf, 'Name', sprintf('Suppression Index: Perception vs. Working Memory'));
    hold all;
    %Perception
    errorbar(centerContrast', mean(suppressionIndex.Collapsed(:,:,1),1), std(suppressionIndex.Collapsed(:,:,1),1)...
        /sqrt(length(subjects)), 'r','LineWidth',3);
    %Working Memory
    errorbar(centerContrast', mean(suppressionIndex.Collapsed(:,:,2),1), std(suppressionIndex.Collapsed(:,:,2),1)...
        /sqrt(length(subjects)), 'b','LineWidth',3);
    plot(repmat(centerContrast', [length(subjects) 1]), (suppressionIndex.Collapsed(:,:,1)),'r.','MarkerSize',10);
    plot(repmat(centerContrast', [length(subjects) 1]), (suppressionIndex.Collapsed(:,:,2)),'b.','MarkerSize',10);
    hold on
    plot([0.09 0.8],[0 0],':k')
    ylabel('Suppression Index (surround-nosurround)/(surround+nosurround)');
    xlabel('Contrast (%)');
    xlim([0.09 0.8]);
    ylim([-0.45 0.45]);
    xticks(C_Test)
    axis square;
    legend({'Perception Condition','Working Memory Condition'})
    set(gca,'XScale','log')
    
   
    
    
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
    ylim([-0.5 0.5]);
    title('Order 1: Location, then Contrast');
    axis square;
    legend({'Perception Condition','Working Memory Condition'})
end

%% Group Plots
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
    
    
%     % Total Contrast Estimates for Selective probes - exlcudes 0.75 and 0.1
%     % probe trials
%     plotContrasts = 100*round(centerContrast,2);
%     figure;
%     cvp_group = loglog(plotContrasts, 100.*squeeze(mean(ContrastEstimate_PS)));
%     hold on
%     errorbar(repmat(plotContrasts,1,3),100.*squeeze(mean(ContrastEstimate_PS)),100*squeeze(std(ContrastEstimate_PS)./sqrt(size(ContrastEstimate_PS,1))),'k','LineStyle','None')
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
%     title(['Center vs. Perceived Contrast, without 0.75 or 0.1 Probe Trials'])
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
    set(gca,'XScale','log')
    
%     % Split between orders%
%     figure('Color', [1 1 1]);
%     
%     subplot(1,2,1)
%     set(gcf, 'Name', sprintf('Suppression Index: Perception vs. Working Memory'));
%     hold all;
%     %Perception
%     errorbar(centerContrast', nanmean(suppressionIndex.Order0(:,:,1),1), nanstd(suppressionIndex.Order0(:,:,1),1)...
%         /sqrt(numel(subjectProfile.SubjectName)), 'r','LineWidth',3)
%     %Working Memory
%     errorbar(centerContrast', nanmean(suppressionIndex.Order0(:,:,2),1), nanstd(suppressionIndex.Order0(:,:,2),1)...
%         /sqrt(numel(subjectProfile.SubjectName)), 'b','LineWidth',3); hold all;
%     plot(repmat(centerContrast', [numel(subjectProfile.SubjectName) 1]), (suppressionIndex.Order0(:,:,1)),'r.','MarkerSize',10);
%     plot(repmat(centerContrast', [numel(subjectProfile.SubjectName) 1]), (suppressionIndex.Order0(:,:,2)),'b.','MarkerSize',10);
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
%         /sqrt(numel(subjectProfile.SubjectName)), 'r','LineWidth',3);
%     %Working Memory
%     errorbar(centerContrast', mean(suppressionIndex.Order1(:,:,2),1), std(suppressionIndex.Order1(:,:,2),1)...
%         /sqrt(numel(subjectProfile.SubjectName)), 'b','LineWidth',3);
%     plot(repmat(centerContrast', [numel(subjectProfile.SubjectName) 1]), (suppressionIndex.Order1(:,:,1)),'r.','MarkerSize',10);
%     plot(repmat(centerContrast', [numel(subjectProfile.SubjectName) 1]), (suppressionIndex.Order1(:,:,2)),'b.','MarkerSize',10);
%     hold all
%     plot([0.09 0.8],[0 0],':k')
%     ylabel('Suppression Index (surround-nosurround)/(surround+nosurround)');
%     xlabel('Contrast (%)');
%     xlim([0.09 0.8]);
%     ylim([-0.3 0.3]);
%     title('Order 1: Location, then Contrast');
%     axis square;
%     legend({'Perception Condition','Working Memory Condition'})
    
% %    Suppression index: excluding trials with 0.1 or 0.75 contrasts %
%     hold off
%     figure('Color', [1 1 1]);
%     set(gcf, 'Name', sprintf('Suppression Index: Perception vs. Working Memory, without 0.1 or 0.75 contrast probes'));
%     hold all;
%     %Perception
%     errorbar(centerContrast', mean(suppressionIndex.Collapsed_PS(:,:,1),1), std(suppressionIndex.Collapsed_PS(:,:,1),1)...
%         /sqrt(size(visualmemory_subjectsRan,2)), 'r','LineWidth',3);
%     %Working Memory
%     errorbar(centerContrast', mean(suppressionIndex.Collapsed_PS(:,:,2),1), std(suppressionIndex.Collapsed_PS(:,:,2),1)...
%         /sqrt(size(visualmemory_subjectsRan,2)), 'b','LineWidth',3);
%     plot(repmat(centerContrast', [numel(subjectProfile.SubjectName) 1]), (suppressionIndex.Collapsed_PS(:,:,1)),'r.','MarkerSize',10);
%     plot(repmat(centerContrast', [numel(subjectProfile.SubjectName) 1]), (suppressionIndex.Collapsed_PS(:,:,2)),'b.','MarkerSize',10);
%     hold on
%     plot([0.09 0.8],[0 0],':k')
%     ylabel('Suppression Index (surround-nosurround)/(surround+nosurround)');
%     xlabel('Contrast (%)');
%     xlim([0.09 0.8]);
%     ylim([-0.3 0.3]);
%     axis square;
%     legend({'Perception Condition','Working Memory Condition'})
end



%% sanity check
if groupPlots
figure;
    for i = 1:length(subjects)
        subplot(2,5,i)
        loglog(centerContrast,subjectProfile.ContrastEstimate(i,:,1))
        hold all
        loglog(centerContrast,subjectProfile.ContrastEstimate(i,:,2))
        loglog(centerContrast,subjectProfile.ContrastEstimate(i,:,3))
        legend({'Perc','WM','BL'})
        xlim([0.1 0.75]); ylim([0.1 0.75])
        plot([0.1 0.75],[0.1 0.75],'k--')
    end
end