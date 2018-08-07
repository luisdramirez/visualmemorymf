close all
clear all

%% SETUP %%

%expDir = '/Users/juliaschwartz/Desktop/visualmemorymf'; %Lab computer
expDir = '/Users/julia/Desktop/Ling Lab/Experiments/visualmemorymf'; %Laptop
p = genpath(expDir);
addpath(p);
%dataDir = '/Users/juliaschwartz/Desktop/visualmemorymf/data_master'; %Lab computer
dataDir = '/Users/julia/Desktop/Ling Lab/Experiments/visualmemorymf/data_master'; %Laptop
cd(dataDir)
indvFiguresOn = 1; %Display plots
load('visualmemory_subjectsRan.mat')
load('data_visualmemorymf_overallData.mat')
experiments = {'Perception Condition','Visual Working Memory'};

% Load Subject files
totalData = cell(1:numel(visualmemory_subjectsRan));
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

numContrasts = 5;
numExperiments = numel(experiments); % 2 experiments now, perception is 1 and working memory is 2...
numSubjects = round(numel(visualmemory_subjectsRan)); % 4 subjects
TotalSuppressionIndexVariable = NaN(numExperiments,numSubjects, numContrasts);
subjAge = [19,19,26,26];

% Model fit setup
C_Test = 10.^(linspace(log10(0.1),log10(0.75), 5)); %Equivalent to center contrast variable.
C_Surround = 1;
C50 = .6;
Wi = 0;
Wi_var = 0;
b = 0;
n = 2;

%% EXPERIMENT LOOP %%
for e = 1:numel(experiments)
    subjects = (1:numel(visualmemory_subjectsRan)); %number of subjects
    subjCount = 0;
    
    %Figure Title
    if indvFiguresOn
        figure(e)
        if e  == 1
            set(gcf, 'Name', 'Perception', 'Color', [1 1 1])
        elseif e == 2
            set(gcf, 'Name', 'Visual Working Memory', 'Color', [1 1 1])
        end
    end
    
    %% SUBJECT LOOP %%
    for subject = subjects
        subjCount = subjCount + 1; % Subject Counter

        % our trial duration is an estimate, calculations for actual trial
        % duration.
        % responseDurations(e,subjCount) = mean(trialDurations);
        
        % Perception or vWM matrices for average baseline/variable
        % responses per subject.
        if e == 1
            % Perception Experiment
            baselineMat(subjCount,:) = overallData.baselineForP(subjCount,:);
            variableMat(subjCount,:) = overallData.perceptionmat(subjCount,:);
        elseif e == 2
            % Working Memory Experiment
            baselineMat(subjCount,:) = overallData.baselineForWM(subjCount,:);
            variableMat(subjCount,:) = overallData.workingmemmat(subjCount,:);
        end
          
        % fit individual data with Normalization model
        options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000);
        
        % Use formula from Xing&Heeger 2001 to use baseline data to constrain alpha and n, and look for the inhibitory weight 
        % to explain the suppressive influence of the surround 
        % Wi: suppression weight
        % C50: inflection point
        % n: nonlinear transducer, determining steepness
        
        indv_r2_coll = zeros(2,numel(subjects),2);
        startValues = [C50 n Wi_var];
       
        %% FIT: Variable and Baseline Conditions per experiment %%
        Data = {C_Test, baselineMat(subjCount,:), variableMat(subjCount,:), C_Surround};
        [est_params(e,subjCount,:), r2(e,subjCount)] = fminsearch('fitNormalizationModel_contrastMatch', startValues, options, Data);
        % Estimated parameters:
        % rows(experiment),col(subject),page(contrast?)
        
        % Fit Y data based off of estimated parameters
        Y_base = (C_Test.^est_params(e,subjCount,2)) ./ ... 
            ((est_params(e,subjCount,1).^est_params(e,subjCount,2)) + (C_Test.^est_params(e,subjCount,2)));
        
        Y_var1 = (C_Test.^est_params(e,subjCount,2)) ./ ... 
            ((est_params(e,subjCount,1).^est_params(e,subjCount,2)) + (C_Test.^est_params(e,subjCount,2)) + (est_params(e,subjCount,3))*(C_Surround.^est_params(e,subjCount,2)));
        
        % R^2 values for baseline and variable condition
        r2_BL = 1 - (sum((baselineMat(subjCount,:) - Y_base).^2) / sum((baselineMat(subjCount,:) - mean(baselineMat(subjCount,:))).^2));
        r2_V = 1 - (sum((variableMat(subjCount,:) - Y_var1).^2) / sum((variableMat(subjCount,:) - mean(variableMat(subjCount,:))).^2));

        indvR2(e,subjCount,:) = [r2_BL r2_V];
        
        %% PLOTTING %%
        
        %Difference in perceived contrast compared to center contrast
        if indvFiguresOn
            C_fit = 0.1:0.01:0.8;
            Y_var = (C_fit.^est_params(e,subjCount,2)) ./ ... 
            ((est_params(e,subjCount,1).^est_params(e,subjCount,2)) + (C_fit.^est_params(e,subjCount,2)) + (est_params(e,subjCount,3)*(C_Surround.^est_params(e,subjCount,2)))); 

            %Display Fits 
            subplot(2,round(numel(visualmemory_subjectsRan)/2), subjCount)
             
            if e == 1
                errorbar(C_Test,Y_var1,overallData.perror(subjCount,:),'ob');
            elseif e == 2
                errorbar(C_Test,Y_var1,overallData.wmerror(subjCount,:),'or');
            end
            hold all;
            loglog(C_fit, Y_var, 'r')
            %plot(C_fit, Y_base, 'k')
            
            %plot(C_Test,variableMat(subjCount,:))
            % DISPLAY DATA WITH ERROR BARS, COME BACK AND EXPORT FULL DATA
            % dont have meancoll - need to go back into data and analysis
            % and export !!!!!
            % errorbar(C_Test, runscollmean(subjCount,:), std(meancoll), 'or');
            % errorbar(TheData(runs).p.testContrasts, runsorthmean(subjCount,:), std(meanorth) , 'ob');            
            % errorbar(TheData(runs).p.testContrasts, runsbasemean(subjCount,:), std(meanbase), '.k');
            hold on
            ylim([0.05 0.8]); xlim([0.09 0.8]); box off
            ylabel('Perceived Contrast')
            xlabel('Center Contrast');
            plot([0.09 0.8], [0.09, 0.8], 'k-');
            title(['Subject ' num2str(subjCount)]);
            axis square;
            legend('Variable Condition','Baseline Condition') 
        end
    end % END OF SUBJECT LOOP %
    
    % Avg. Data pooled over subjects %
    if e == 1
        subjectsbasemean{e} = (overallData.baselineForPMean);
        subjectsvarmean{e} = (overallData.perceptionmean);    
    elseif e == 2
        subjectsbasemean{e} = (overallData.baselineForWMMean);
        subjectsvarmean{e} = (overallData.workingmemmean); 
    end
    
    %Suppression Index over subjects    
    TotalSuppressionIndexVariable(e,1:numel(visualmemory_subjectsRan),:) = ((variableMat) - (baselineMat))./((variableMat) + (baselineMat)); 
    
    totalbaselineMat(:,:,e) = baselineMat;
    totalvarMat(:,:,e) = variableMat;
end
%% END OF EXPERIMENT LOOP %%

%Plot contrast
figure('Color', [1 1 1]);
hold all;
for e = 1:numel(experiments)
    if e == 1, color = 'or'; else , color = ' ob'; end %red is perception, blue is wm
    plot(repmat(C_Test, [4 1]), (squeeze(TotalSuppressionIndexVariable(e,:,:))),color);
    ylabel('Suppression Index ( (surround-nosurround)/(surround+nosurround) )')
    xlabel('Contrast (%)')
    set(gca, 'XTick', C_Test, 'XTickLabel', round(C_Test*100), 'XScale', 'log')
    xlim([0.09 0.8]); ylim([-0.3 0.3]); axis square
    plot([0.09 0.8], [0 0], ':k')
end

for e = 1:numel(experiments),errorbar(C_Test, nanmean(squeeze(TotalSuppressionIndexVariable(e,:,:))), ...
        nanstd(squeeze(TotalSuppressionIndexVariable(e,:,:)))/sqrt(numel(visualmemory_subjectsRan)), 'ko-'); end
legend({'Perception Error','WM Error','Perception','WM'})

%%% Plot model Parameters
% R^2
figure('Color', [1 1 1])
title('R^2 Value Per subject Condition')
subplot(1,2,1)
bar([indvR2(1,:,1)' indvR2(1,:,2)']') %first page is baseline, first row is perception condition
set(gca, 'xtickLabel', {'Baseline' 'Perception'})
ylabel('R2'), box off; title('Perception');
subplot(1,2,2)
bar([indvR2(2,:,1)' indvR2(2,:,2)']')
set(gca, 'xtickLabel', {'Baseline' 'Working Memory'})
ylabel('R2'), legend(num2str((1:4)')); box off; title('Visual working memory');

% c50 - both perception and working memory conditions
figure('Color', [1 1 1])
subplot(1,4,1)
bar(1:2, mean(squeeze(est_params(:,:,1)),2)') %page 1 is c50
hold all
handles = get(gca, 'Children');
set(handles(1), 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none'); 
xlim([0 3]); box off
errorbar([1 2], [mean(squeeze(est_params(1,:,1))) mean(squeeze(est_params(2,:,1)))], ...
    [std(squeeze(est_params(1,:,1)))/sqrt(subjCount) std(squeeze(est_params(2,:,1)))/sqrt(subjCount)], 'k.')
plot(repmat([1 2], [4 1]), (squeeze(est_params(:,:,1)))', 'ok')
set(gca, 'Xtick', [1 2], 'XtickLabel', {'Perc' 'vWM'})
title('C50')

subplot(1,4,2)
% n - both perception and working memory conditions
bar(1:2, mean(squeeze(est_params(:,:,2)),2)) %page 2 is n
hold all
handles = get(gca, 'Children');
set(handles(1), 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none'); 
xlim([0 3]); box off
errorbar([1 2], [mean(squeeze(est_params(1,:,2))) mean(squeeze(est_params(2,:,2)))], ...
    [std(squeeze(est_params(1,:,2)))/sqrt(subjCount) std(squeeze(est_params(2,:,2)))/sqrt(subjCount)], 'k.')
plot(repmat([1 2], [4 1]), (squeeze(est_params(:,:,2)))', 'ok')
set(gca, 'Xtick', [1 2], 'XtickLabel', {'Perc' 'vWM'})
title('n')

subplot(1,4,3)
% Wi est - both perception and working memory conditions
bar([0 1.3], mean(est_params(:,:,3),2), 0.3);
hold all
handles = get(gca, 'Children');
set(handles(1), 'FaceColor', [0 0 1], 'EdgeColor', 'none'); 
errorbar([0 1.3], [mean(est_params(:,:,3),2)'], ...
    [std(est_params(:,:,3),[], 2)'/sqrt(subjCount)], 'k.')
plot(repmat([0 1.3 ], [4 1]), [est_params(:,:,3)'], 'ok');
box off; xlim([-0.5 1.5]); set(gca, 'Xtick', [0.25 1.55], 'XtickLabel', {'Perc' 'vWM'})
title('Parameter surround conditions')
%the third and fourth page were collinear (our variable) and orthogonal,
%fourth deleted


% h1, h2
subplot(1,4,4)
h1 = scatter(est_params(1,:,3), est_params(2,:,3));
hold all
set(h1, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none')
%h2 = scatter(est_params(1,:,4), est_params(2,:,4));
%set(h2, 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', 'none')
ylim([-0.25 0.25]); xlim([-0.1 0.4]); axis square
plot([0 0], [-0.4 0.4], ':k', [-0.4 0.4], [0 0 ], ':k')
ylabel('Visual working memory'), xlabel('Perception'); title('Gamma estimate')
legend('Variable')



%% Final figure paper: plot perceived contrast over observers, with estimated fit
figure('Color', [1 1 1], 'Name', 'Perceived Contrast')
for  n = 1:2 % EITHER perception / vwm
    subplot(1,2,n)
    Y_var(n,:) = (C_fit.^mean(est_params(n,:,2))) ./ ...
        ((mean(est_params(n,:,1)).^mean(est_params(n,:,2))) + (C_fit.^mean(est_params(n,:,2))) + (mean(est_params(n,:,3))*(C_Surround.^mean(est_params(n,:,2)))));
    y_data_var = mean(subjectsvarmean{n});
    y_data_base = mean(subjectsbasemean{n}); 
    % plot
    loglog(C_Test, y_data_var, '.r')
    hold all 
    loglog(C_Test, y_data_base, '.k')
    errorbar(C_Test, subjectsvarmean{n}, 1*(std(totalvarMat(:,:,n))/sqrt(numSubjects)), '-or')
    errorbar(C_Test, subjectsbasemean{n}, 1*(std(totalbaselineMat(:,:,n))/sqrt(numSubjects)), '-ok')
    plot([0.1 0.8], [0.1 0.8], ':k')
    xlim([0.1 0.8]), ylim([0.1 0.8])
    if  n ==1, title('Perception'), legend({'Collinear', 'No surround'})
    else ,title('Visual working memory'), end
    ylabel('Perceived contrast'); xlabel('Center Contrast'); axis square; box off
end

%% STATISTICS %%

% Test whether there is a difference in baseline parameters between the different surround configurations
% -------------------- %

% c50 - between perception and working memory
[h_c50,p_c50,ci_c50,st_c50] = ttest(est_params(1,:,1), est_params(2,:,1));
% n - between perception and working memory
[h_n,p_n,ci_n,st_n] = ttest(est_params(1,:,2), est_params(2,:,2));

% Test whether there is a difference in the inhibitory weight
% inhibitory weight - between perception and working memory
[h_wi,p_wi,ci_wi,st_wi] = ttest(est_params(1,:,2), est_params(2,:,2));

% Should do a 2x2 repeated measures anova.. rule out there is no difference between perception and vwm either 

% Test all weights against 0
[h_ttest,p_ttest, ci_ttest, st_ttest] = ttest([est_params(1,:,3)' est_params(2,:,3)'], [], 'tail', 'right');

% % Do power analysis permutation[perc_coll perc_orth vwm_coll vwm_orth]
% [power, cohens_d] = powerAnalysis_exp1([squeeze(est_params(1,:,3)) squeeze(est_params(2,:,3))], 10000);
% 
% noPerm_cohens_d = st_ttest.tstat./sqrt(numSubjects);
% 


