
% trialEvents = TheData(runs).p.TrialEvents;
% differenceContrast = TheData(1).data.DifferenceContrast;

close all
clear all

% expDir = '/Users/ywatanabe/Dropbox/NormalizationVisualMemory/Experiment1';
expDir = '/Users/juliaschwartz/Desktop/visualmemorymf';
p = genpath(expDir);
addpath(p);
dataDir = '/Users/juliaschwartz/Desktop/visualmemorymf/data_master';
cd(dataDir)
indvFiguresOn = 1;

experiment = 'Exp 1';
load('visualmemory_subjectsRan.mat')
load('data_visualmemorymf_overallData.mat')

numSubjects = length(visualmemory_subjectsRan);
numContrast = 5;
numExperiments = numel(experiment);
TotalSuppressionIndexColinear = NaN(numExperiments,numSubjects, numContrast);
TotalSuppressionIndexOrthogonal = NaN(numExperiments,numSubjects, numContrast); % Not used here %
subjAge = [19,19,24,24]; % Figure out and correct order %


% Model fit setup
C_Test = 10.^(linspace(log10(0.1),log10(0.75), 5)); %center contrast
C_Surround = 1; %surround contrast
C50 = .6; 
Wi = 0; 
Wi_coll = 0; % Collinear orientation %
Wi_orth = 0; % ignore, we didnt do orthoganal orientation %
b = 0;
n = 2; %nonlinear transducer (steepness of the function)

for e = 1:numel(experiment)
    
    % Subject Loop %
    if indvFiguresOn
        figure(e), if e  == 1, set(gcf, 'Name', 'Perception', 'Color', [1 1 1]), else set(gcf, 'Name', 'Visual working memory', 'Color', [1 1 1]), end
    end
    subjects = (1:4);
    % Load subject Data %
    files = struct2cell(dir(dataDir))';
    [numFiles, ~] = size(files);
    possibleFileNames = cell(length(visualmemory_subjectsRan),1);
    for i = 1:length(visualmemory_subjectsRan)
        filename = strcat('data_visualmemorymf_exp_',visualmemory_subjectsRan{i},'.mat');
        possibleFileNames{i,1} = filename;
    end
    for i = 1:numel(possibleFileNames)
        if exist(possibleFileNames{i,1},'file') == 1 ||exist(possibleFileNames{i,1},'file') == 2
            load(possibleFileNames{i,1})
            totalData{i} = theData;
            fprintf('\nLoading %s',possibleFileNames{i,1})
        end
    end

    meancoll = []; % Mean Collinear %
    meanorth = []; % Mean Orthoganal %
    meanbase = []; % Mean Baseline %
    trialDurations = [];
    runDurations = [];
        
        for runs = 1:length(theData)
            %accuracies = TheData(runs).data.DifferenceContrast;
            accuracies = theData(runs).data.EstimatedContrast;
            contrast = theData(runs).p.trialEvents(:,3); %Actual contrast
            %orientation = theData(runs).p.TrialEvents(:,1); -- We dont
            %have orientation
            trialDurations = [trialDurations; theData(runs).data.ResponseTime_Contrast]; %Either response time contrast or response time location, both??
            %runDurations = [runDurations; theData(runs).t.EndTime];
%             collinear = orientation == 1;
%             orthogonal = orientation == 2;
%             base = orientation == 3;
            
            %number of elements in this array
            numContrasts = numel(theData(runs).p.numContrasts);
            
%             for c = 1:numContrasts
%                 %index of which trials have a specific contrast in them
%                 %finding accuracies with collinear orientations and 5 different
%                 %contrasts
%                 collcontrasts(:,c) = (orientation == 1) & (contrast == TheData(runs).p.testContrasts(c));
%                 collaccuracies(:,c) = (accuracies(collcontrasts(:,c)));
%                 
%                 %finding accuracies with orthogonal orientation and 5 different
%                 %contrasts
%                 orthcontrasts(:,c) = (orientation == 2) & (contrast == TheData(runs).p.testContrasts(c));
%                 orthaccuracies(:,c) = (accuracies(orthcontrasts(:,c)));
%                 
%                 %finding accuracies with no surround and 5 different contrasts
%                 basecontrasts(:,c) = (orientation == 3) & (contrast == TheData(runs).p.testContrasts(c));
%                 baseaccuracies(:,c) = (accuracies(basecontrasts(:,c)));
%                 
%                 
%             end
%             
            % collect all the responses over multiple runs
%             meancoll = [meancoll; collaccuracies];
%             %meanorth = [meanorth; orthaccuracies];
%             meanbase = [meanbase; baseaccuracies];
%             
%             clear collaccuracies orthaccuracies baseaccuracies collcontrasts orthcontrasts basecontrasts
%             
        end
        
        
        %responseDurations(e,subjCount) = mean(trialDurations);
        
%         %mean over runs for baseline
%         runsbasemean(subjCount,:) = mean(meanbase); %Here do we want 5 values?
%         %mean over runs for collinear
%         runscollmean(subjCount,:) = mean(meancoll);
%         %mean over runs for orthogonal
%         runsorthmean(subjCount,:) = mean(meanorth);
        
%         runsbasemean = overallData.baselineForP;%Make sure its only the corresponding baseline data to the condition: wm/p
%         runscollmean = overallData.perceptionmean; %change with workingmemmean

%           runsbasemean = overallData.baselineforP;
          runsbasemean = overallData.baselineForWM;
          runscollmean = overallData.workingmemmat;
          runscollmean = overallData.perceptionmat;
       for subjCount = 1:numel(visualmemory_subjectsRan) 
        % fit individual data with Normalization model
        options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000);
        
        % Use formula from Xing&Heeger 2001 to use baseline data to constrain alpha and n, and look for the inhibitory weight 
        % to explain the suppressive influence of the surround 
        indv_r2_coll = zeros(2,numel(subjects),2);
        startValues = [C50 n Wi_coll Wi_orth];
        
        % Fit Collinear and Orthogonal surround conditions together
        
        Data = {C_Test, runsbasemean(subjCount,:), runscollmean(subjCount,:), C_Surround};
        [est_params(e,subjCount,:), r2(e,subjCount)] = fminsearch('fitNormalizationModel_contrastMatch', startValues, options, Data);

        Y_base = (C_Test.^est_params(e,subjCount,2)) ./ ... 
            ((est_params(e,subjCount,1).^est_params(e,subjCount,2)) + (C_Test.^est_params(e,subjCount,2)));
        
        Y_coll = (C_Test.^est_params(e,subjCount,2)) ./ ... 
            ((est_params(e,subjCount,1).^est_params(e,subjCount,2)) + (C_Test.^est_params(e,subjCount,2)) + (est_params(e,subjCount,3))*(C_Surround.^est_params(e,subjCount,2)));
      
        Y_orth = (C_Test.^est_params(e,subjCount,2)) ./ ... 
            ((est_params(e,subjCount,1).^est_params(e,subjCount,2)) + (C_Test.^est_params(e,subjCount,2)) + (est_params(e,subjCount,4))*(C_Surround.^est_params(e,subjCount,2)));
      
        r2_m = 1 - (sum((runsbasemean(subjCount,:) - Y_base).^2) / sum((runsbasemean(subjCount,:) - mean(runsbasemean(subjCount,:))).^2));
        r2_c = 1 - (sum((runscollmean(subjCount,:) - Y_coll).^2) / sum((runscollmean(subjCount,:) - mean(runscollmean(subjCount,:))).^2));
%         r2_o = 1 - (sum((runsorthmean(subjCount,:) - Y_orth).^2) / sum((runsorthmean(subjCount,:) - mean(runsorthmean(subjCount,:))).^2));
        
%         indvR2(e,subjCount,:) = [r2_m r2_c r2_o]; no othoganal
indvR2(e,subjCount,:) = [r2_m r2_c];

        
        % plots
        %Difference in perceived contrast compared to test contrast
        if indvFiguresOn
            
            C_fit = 0.1:0.01:0.8;
            
            Y_coll = (C_fit.^est_params(e,subjCount,2)) ./ ... 
            ((est_params(e,subjCount,1).^est_params(e,subjCount,2)) + (C_fit.^est_params(e,subjCount,2)) + (est_params(e,subjCount,3)*(C_Surround.^est_params(e,subjCount,2))));
            
            Y_orth = (C_fit.^est_params(e,subjCount,2)) ./ ... 
            ((est_params(e,subjCount,1).^est_params(e,subjCount,2)) + (C_fit.^est_params(e,subjCount,2)) + (est_params(e,subjCount,4)*(C_Surround.^est_params(e,subjCount,2))));
        
                      
            subplot(2,round(numel(visualmemory_subjectsRan)/2), subjCount)
            
            % display fits
            loglog(C_fit, Y_coll, 'r') 
            hold all;
%             plot(C_fit, Y_orth, 'b')
%             plot(C_fit, Y_base, 'k')
            % display data with errorbars
%             errorbar(C_Test, runscollmean(subjCount,:), std(meancoll),
%             'or'); %%% ADD MEANCOLL BACK IN
%             errorbar(theData(runs).p.testContrasts, runsorthmean(subjCount,:), std(meanorth) , 'ob');            
%             errorbar(TheData(runs).p.testContrasts, runsbasemean(subjCount,:), std(meanbase), '.k');
            ylim([0.05 0.8]); xlim([0.09 0.8]); box off
            % ylabel('Difference in perceived contrast')
%             ylabel('Perceived contrast'); xlabel('Contrast');
            plot([0.09 0.8], [0.09, 0.8], 'k-');
            title(['Subject ' num2str(subjCount)]);
            axis square;
            if subjCount == numel(subjects)
                legend('Collinear')
            end
                            
        end
        
        runTime(e,subjCount) = mean(runDurations);
        
        
    end
end
    %data pooled over subjects
    subjectsbasemean{e} = (runsbasemean);
    subjectscollmean{e} = (runscollmean);
%     subjectsorthmean{e} = (runsorthmean);    
    
    %Suppression Index over subjects    
    TotalSuppressionIndexColinear(e,1:numel(visualmemory_subjectsRan),:) = ((runscollmean) - (runsbasemean))./((runscollmean) + (runsbasemean));
%     TotalSuppressionIndexOrthogonal(e,1:numel(visualmemory_subjectsRan),:) = ((runsorthmean) - (runsbasemean))./((runsorthmean) + (runsbasemean));
    
    clear runsbasemean runscollmean runsorthmean


cd ..


figure('Color', [1 1 1]);
% subplot(2,2,1)
hold all;

errorbar(C_Test, nanmean(squeeze(TotalSuppressionIndexColinear(1,:,:))), ...
    nanstd(squeeze(TotalSuppressionIndexColinear(1,:,:)))/sqrt(numel(visualmemory_subjectsRan)), 'ro-');
errorbar(C_Test, nanmean(squeeze(TotalSuppressionIndexColinear(2,:,:))), ...
    nanstd(squeeze(TotalSuppressionIndexColinear(2,:,:)))/sqrt(numel(visualmemory_subjectsRan)), 'r:^');
errorbar(C_Test, nanmean(squeeze(TotalSuppressionIndexOrthogonal(1,:,:))), ...
    nanstd(squeeze(TotalSuppressionIndexOrthogonal(1,:,:)))/sqrt(numel(visualmemory_subjectsRan)), 'bo-');
errorbar(C_Test, nanmean(squeeze(TotalSuppressionIndexOrthogonal(2,:,:))), ...
    nanstd(squeeze(TotalSuppressionIndexOrthogonal(2,:,:)))/sqrt(numel(visualmemory_subjectsRan)), 'b:^');
% plot(repmat(TheData(runs).p.testContrasts, [12 1]), (squeeze(TotalSuppressionIndexColinear(1,:,:))), 'ro');
% plot(repmat(TheData(runs).p.testContrasts, [12 1]), (squeeze(TotalSuppressionIndexColinear(2,:,:))), 'ro');
% plot(repmat(TheData(runs).p.testContrasts, [12 1]), (squeeze(TotalSuppressionIndexOrthogonal(1,:,:))), 'bo');
% plot(repmat(TheData(runs).p.testContrasts, [12 1]), (squeeze(TotalSuppressionIndexOrthogonal(2,:,:))), 'bo');

legend({'Coll Surround Suppression' 'Coll vWM Suppression' ...
    'Orth Surround Suppression' 'Orth vWM Suppression'})
ylabel('Suppression Index ( (surround-nosurround)/(surround+nosurround) )')
xlabel('Contrast (%)')
set(gca, 'XTick', C_Test, 'XTickLabel', round(C_Test*100), 'XScale', 'log')
xlim([0.09 0.8]); ylim([-0.15 0.15]); axis square

SubjectDifferenceCollinear = squeeze(TotalSuppressionIndexColinear(1,:,:))- ...
    squeeze(TotalSuppressionIndexColinear(2,:,:));
% SubjectDifferenceOrthogonal = squeeze(TotalSuppressionIndexOrthogonal(1,:,:)) - ...
%     squeeze(TotalSuppressionIndexOrthogonal(2,:,:));



%%% Plot model Parameters
% model fits
figure('Color', [1 1 1])
subplot(1,2,1)
bar([indvR2(1,:,1)' indvR2(1,:,2)']')
set(gca, 'xtickLabel', {'Baseline' 'Collinear' 'Orthogonal'})
ylabel('R2'), box off; title('Perception');
subplot(1,2,2)
bar([indvR2(2,:,1)' indvR2(2,:,2)' ]')
set(gca, 'xtickLabel', {'Baseline' 'Collinear' 'Orthogonal'})
ylabel('R2'), legend(num2str([1:12]')); box off; title('Visual working memory');


figure('Color', [1 1 1])
subplot(1,4,1)
% c50 -
bar(1:2, mean(squeeze(est_params(:,:,1)),2))
hold all
handles = get(gca, 'Children');
set(handles(1), 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none'); 
xlim([0 3]); box off
errorbar([1 2], [mean(squeeze(est_params(1,:,1))) mean(squeeze(est_params(2,:,1)))], ...
    [std(squeeze(est_params(1,:,1)))/sqrt(subjCount) std(squeeze(est_params(2,:,1)))/sqrt(subjCount)], 'k.')
plot(repmat([1 2], [12 1]), (squeeze(est_params(:,:,1)))', 'ok')
set(gca, 'Xtick', [1 2], 'XtickLabel', {'Perc' 'vWM'})
title('C50 - no surround condition')

subplot(1,4,2)
% c50 -
bar(1:2, mean(squeeze(est_params(:,:,2)),2))
hold all
handles = get(gca, 'Children');
set(handles(1), 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none'); 
xlim([0 3]); box off
errorbar([1 2], [mean(squeeze(est_params(1,:,2))) mean(squeeze(est_params(2,:,2)))], ...
    [std(squeeze(est_params(1,:,2)))/sqrt(subjCount) std(squeeze(est_params(2,:,2)))/sqrt(subjCount)], 'k.')
plot(repmat([1 2], [12 1]), (squeeze(est_params(:,:,2)))', 'ok')
set(gca, 'Xtick', [1 2], 'XtickLabel', {'Perc' 'vWM'})
title('n - no surround condition')

subplot(1,4,3)
bar([0 1.3], mean(est_params(:,:,3),2), 0.3);
hold all
bar([0.5 1.8], mean(est_params(:,:,4),2), 0.3)
handles = get(gca, 'Children');
set(handles(1), 'FaceColor', [0 0 1], 'EdgeColor', 'none'); set(handles(2), 'FaceColor', [1 0 0], 'EdgeColor', 'none')
errorbar([0 1.3 0.5 1.8], [mean(est_params(:,:,3),2)' mean(est_params(:,:,4),2)'], ...
    [std(est_params(:,:,3),[], 2)'/sqrt(subjCount) std(est_params(:,:,4),[],2)'/sqrt(subjCount)], 'k.')
plot(repmat([0 1.3 0.5 1.8], [12 1]), [est_params(:,:,3)' est_params(:,:,4)'], 'ok');
box off; xlim([-0.5 2.5]); set(gca, 'Xtick', [0.25 1.55], 'XtickLabel', {'Perc' 'vWM'})
title('Parameter surround conditions')

subplot(1,4,4)
h1 = scatter(est_params(1,:,3), est_params(2,:,3));
hold all
set(h1, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none')
h2 = scatter(est_params(1,:,4), est_params(2,:,4));
set(h2, 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', 'none')
ylim([-0.25 0.25]); xlim([-0.1 0.4]); axis square
plot([0 0], [-0.4 0.4], ':k', [-0.4 0.4], [0 0 ], ':k')
ylabel('Visual working memory'), xlabel('Perception'); title('Gamma estimate')
legend({'Collinear' 'Orthogonal'})


cd(expDir)

%% Final figure paper: plot perceived contrast over observers, with estimated fit
figure('Color', [1 1 1], 'Name', 'Perceived Contrast')
for n = 1:2 % perception and vWM
    subplot(1,3,n)
    
    Y_coll(n,:) = (C_fit.^mean(est_params(n,:,2))) ./ ...
        ((mean(est_params(n,:,1)).^mean(est_params(n,:,2))) + (C_fit.^mean(est_params(n,:,2))) + (mean(est_params(n,:,3))*(C_Surround.^mean(est_params(n,:,2)))));
    
    y_data_coll = mean(subjectscollmean{n});
    
    Y_orth(n,:) = (C_fit.^mean(est_params(n,:,2))) ./ ...
        ((mean(est_params(n,:,1)).^mean(est_params(n,:,2))) + (C_fit.^mean(est_params(n,:,2))) + (mean(est_params(n,:,4))*(C_Surround.^mean(est_params(n,:,2)))));

    y_data_orth = mean(subjectsorthmean{n});

    Y_data_base = mean(subjectsbasemean{n});
        
    % plot
    loglog(C_Test, y_data_coll, '.r')
    hold all 
    loglog(C_Test, y_data_orth, '.b')
    loglog(C_Test, Y_data_base, '.k')

    errorbar(C_Test, mean(subjectscollmean{n}), 1*(std(subjectscollmean{n})/sqrt(12)), '-or')
    errorbar(C_Test, mean(subjectsorthmean{n}), 1*(std(subjectsorthmean{n})/sqrt(12)), '-ob')
    errorbar(C_Test, mean(subjectsbasemean{n}), 1*(std(subjectsbasemean{n})/sqrt(12)), '-ok')
    
    plot([0.1 0.8], [0.1 0.8], ':k')
    xlim([0.1 0.8]), ylim([0.1 0.8])
    if  n ==1, title('Perception'), legend({'Collinear' 'Orthogonal' 'No surround'})
    else title('Visual working memory'), end
    ylabel('Perceived contrast'); xlabel('Center Contrast'); axis square; box off
end
subplot(1,3,3)
bar([0 1.3], mean(est_params(:,:,3),2), 0.3);
hold all
bar([0.5 1.8], mean(est_params(:,:,4),2), 0.3)
handles = get(gca, 'Children');
set(handles(1), 'FaceColor', [0 0 1], 'EdgeColor', 'none'); set(handles(2), 'FaceColor', [1 0 0], 'EdgeColor', 'none')
errorbar([0 1.3 0.5 1.8], [mean(est_params(:,:,3),2)' mean(est_params(:,:,4),2)'], ...
    [std(est_params(:,:,3),[], 2)'/sqrt(subjCount) std(est_params(:,:,4),[],2)'/sqrt(subjCount)], 'k.')
plot(repmat([0 1.3 0.5 1.8], [12 1]), [est_params(:,:,3)' est_params(:,:,4)'], 'ok');
box off; xlim([-0.5 2.5]); set(gca, 'Xtick', [0.25 1.55], 'XtickLabel', {'Perc' 'vWM'})
title('Parameter surround conditions'); axis square
legend({'Collinear' 'Orthogonal'})

%% Save out data as txt file
fileName = 'exp1_normalizationStrength.txt';
dataLabels = {'Perception_coll', 'Perception_orth', 'Visual_memory_coll', 'Visual_memory_orth'};

if ~exist(fileName,'file')
    T = table(squeeze(est_params(1,:,3))', squeeze(est_params(1,:,4))', squeeze(est_params(2,:,3))', squeeze(est_params(2,:,4))','VariableNames', dataLabels);
    writetable(T, fileName);
end

%% STATISTICS

% Test whether there is a difference in baseline parameters between the different surround configurations

% Should do a 2x2 repeated measures anova.. rule out there is no difference between perception and vwm either 
% Alpha - perception
[~,p_ap,ci_ap,st_ap] = ttest(est_params(1,:,1), est_params(2,:,1))
% n - perception
[~,p_np,ci_np,st_np] = ttest(est_params(1,:,2), est_params(2,:,2))


% Test whether there is a difference in the inhibitory weigth

% % Should do a 2x2 repeated measures anova.. rule out there is no difference between perception and vwm either 
% % inhibitory weigth - collinear
% [~,p_coll,ci_coll,st_coll] = ttest(est_params(1,:,3), est_params(2,:,3))
% % inhibitory weigth - orthogonal
% [~,p_orth,ci_orth,st_orth] = ttest(est_params(1,:,4), est_params(2,:,4))
% 
% % inhibitory weigth - orientation tuned
% [~,p_perc,ci_perc,st_perc] = ttest(est_params(1,:,3), est_params(1,:,4))
% [~,p_vwm,ci_vwm,st_vwm] = ttest(est_params(2,:,3), est_params(2,:,4))

% Test all weigths against 0
[~,p_ttest, ci_ttest, st_ttest] = ttest([est_params(1,:,3)' est_params(1,:,4)' est_params(2,:,3)' est_params(2,:,4)'], [], 'tail', 'right') 

% Do power analysis permutation[perc_coll perc_orth vwm_coll vwm_orth]
[power, cohens_d] = powerAnalysis_exp1([squeeze(est_params(1,:,3:4)) squeeze(est_params(2,:,3:4))], 10000);

noPerm_cohens_d = st_ttest.tstat./sqrt(12)


% save out baseline data in order to compare to Exp 2
save BaselineModelParameters_Exp1.mat est_params


