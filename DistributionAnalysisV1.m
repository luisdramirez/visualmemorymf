clear
close all
% expDir = '/Users/ywatanabe/Dropbox/NormalizationVisualMemory/Experiment1';
expDir = '/Users/thelinglab/Dropbox/NormalizationVisualMemory/Experiment1';
cd([expDir '/Data'])

plotOn = 1;
experiments = {'surrSuppression_2000_Delay*' 'vWM_surrSuppression_2000_Delay*'};

SubjectNumbers = [1:8 10:13];

numSubjects = round(length(dir('*.mat'))/2);
numContrast = 5;
numExperiments = numel(experiments);

%experiment loop
for e = 1:numel(experiments);
    
    list = dir(experiments{e});
    subjCount = 0;
    
    
    %subject loop
    
    for subject = SubjectNumbers
        
        subjCount = subjCount +1;
        
        if subject < 10 && str2num(list(subjCount).name(end-4)) == subject                             
            load(list(subjCount).name);
        elseif subject >= 10 && str2num(list(subjCount).name(end-5:end-4)) == subject
            load(list(subjCount).name);   
        else % subject name does not match the order of the list
            for subj = 1:numel(list) % look in the list to find the right subject
                if subject < 10 
                    if str2num(list(subj).name(end-4)) == subject
                        load(list(subj).name);
                    end
                elseif subject >= 10 
                    if str2num(list(subj).name(end-5:end-4)) == subject
                        load(list(subj).name);
                    end
                end
            end
        end
        
            
        %figure is outside the contrast loop because want each subplot in one
        %figure
        figure ('name', ['Exp. ' num2str(e) ' Subject ' num2str(subject)])
        hold all
        % title(['Distribution Plot Subject ' num2str(subject)]);
        
        
        %number of elements in this array
        numContrasts = numel(TheData(1).p.testContrasts);
        
        for c = 1:numContrasts
            % to set meancoll to blank first so that variable contains something
            responsecoll = [];
            responseorth = [];
            responsebase = [];
            
            for runs = 1:length(TheData);
                
                responses = TheData(runs).data.EstimatedContrast;
                
                contrast = TheData(runs).p.TrialEvents(:,2);
                orientation = TheData(runs).p.TrialEvents(:,1);
                
                collinear = orientation == 1;
                orthogonal = orientation == 2;
                base = orientation == 3;
                
                %finding responses with collinear orientations and 5 different
                % contrasts
                collcontrasts = (orientation == 1) & (contrast == TheData(runs).p.testContrasts(c));
                collresponses = (responses(collcontrasts));
                
                %finding responses with orthogonal orientation and 5 different
                %contrasts
                orthcontrasts = (orientation == 2) & (contrast == TheData(runs).p.testContrasts(c));
                orthresponses = (responses(orthcontrasts));
                
                %finding responses with no surround and 5 different contrasts
                basecontrasts = (orientation == 3) & (contrast == TheData(runs).p.testContrasts(c));
                baseresponses = (responses(basecontrasts));
                
                % collect all the responses over multiple runs
                responsecoll = [responsecoll; collresponses];
                responseorth = [responseorth; orthresponses];
                responsebase = [responsebase; baseresponses];
                
                clear collresponses orthresponses baseresponses collcontrasts orthcontrasts basecontrasts
            end
            
            %collect responses for each contrast into columns (needed because needs
            %to be cleared after every run)
            contrastcoll(:,c) = responsecoll;
            contrastorth(:,c) = responseorth;
            contrastbase(:,c) = responsebase;
            
                        
            % Plot of 3 x (Collinear, Orthogonal, Baseline) by 5 x (Different
            % Contrasts) - done for each subject twice (one for SS and one for vWM)
            
            [N_coll, ~] = histcounts(contrastcoll(:,c), 'BinWidth', 0.05, 'BinLimits',[0,1]);
            [N_orth, ~] = histcounts(contrastorth(:,c), 'BinWidth', 0.05, 'BinLimits',[0,1]);
            [N_base, x] = histcounts(contrastbase(:,c), 'BinWidth', 0.05, 'BinLimits',[0,1]);

            % mean, sigma, amplitude
            startValues = [(TheData(runs).p.testContrasts(c)), (0.1), 10];
            options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000);
            
            [est_params_coll, r2_coll(subjCount,c)] = fminsearch('mygauss', startValues, options, N_coll, x(1:end-1));
            [est_params_orth, r2_orth(subjCount,c)] = fminsearch('mygauss', startValues, options, N_orth, x(1:end-1));
            [est_params_base, r2_base(subjCount,c)] = fminsearch('mygauss', startValues, options, N_base, x(1:end-1));
         
            % substitute baseline - not allowed to freely vary
            est_params_coll(4) = 0;
            est_params_orth(4) = 0;
            est_params_base(4) = 0;
            
            params_coll(subjCount,c,:) = est_params_coll;
            params_orth(subjCount,c,:) = est_params_orth;
            params_base(subjCount,c,:) = est_params_base;
            
            if plotOn == 1
                
                x_fit = 0:0.01:1;
                
                subplot (3,5,c)
                estData_coll = est_params_coll(4) + est_params_coll(3)*exp(-(x_fit-est_params_coll(1)).^2/(est_params_coll(2)^2));
%                 histogram(contrastcoll(:,c),20)
                plot(x(1:end-1), N_coll, 'xr')
                hold on,
                plot(x_fit, estData_coll, 'b')
                xlabel([num2str(round(TheData(1).p.testContrasts(c)*100)) '%']);
                ylabel('Collinear');
                


                subplot(3,5,c+5)
                estData_orth = est_params_orth(4) + est_params_orth(3)*exp(-(x_fit-est_params_orth(1)).^2/(est_params_orth(2)^2));
                plot(x(1:end-1), N_orth, 'xr')
                hold on,
                plot(x_fit, estData_orth, 'b')
%                 histogram(contrastorth(:,c),20)
                xlabel([num2str(round(TheData(1).p.testContrasts(c)*100)) '%' ]);
                ylabel('Orthogonal');
                



                subplot(3,5,c+10)
                estData_base = est_params_base(4) + est_params_base(3)*exp(-(x_fit-est_params_base(1)).^2/(est_params_base(2)^2));
                plot(x(1:end-1), N_base, 'xr')
                hold on,
                plot(x_fit, estData_base, 'b')
%                 histogram(contrastbase(:,c),20)
                xlabel([num2str(round(TheData(1).p.testContrasts(c)*100)) '%']);
                ylabel('Baseline');
                
            end
            
        end
        
        clear contrastcoll contrastorth contrastbase
    end
    
%    TotalSuppressionColinear(e,1:numel(list),:) = ((repmat(TheData(runs).p.testContrasts, [numel(list),1,1]) - runscollmean) - (repmat(TheData(runs).p.testContrasts, [numel(list),1,1]) - runsbasemean)) ;%./  ((repmat(TheData(runs).p.testContrasts, [numel(list),1,1]) - runscollmean) + (repmat(TheData(runs).p.testContrasts, [numel(list),1,1]) - runsbasemean));
%    TotalSuppressionOrthogonal(e,1:numel(list),:) = ((repmat(TheData(runs).p.testContrasts, [numel(list),1,1]) - runsorthmean) - (repmat(TheData(runs).p.testContrasts, [numel(list),1]) - runsbasemean)) ;
    
end


for con = 1:5

coll(:,con) = squeeze(params_coll(:,con,2)); %all subjects, first contrast, STDev
orth(:,con) = squeeze(params_orth(:,con,2));
bases(:,con) = squeeze(params_base(:,con,2));

[h, p, ci, stats] = ttest(coll(:,con), orth(:,con))

[h, p, ci, stats] = ttest(bases(:,con), coll(:,con))

[h, p, ci, stats] = ttest(bases(:,con), orth(:,con))


end 


figure('Color', [1 1 1])
subplot(1,3,1)
bar([-r2_base(:,1) -r2_base(:,2) -r2_base(:,3) -r2_base(:,4) -r2_base(:,5)]')
set(gca, 'xtickLabel', {'10%' '16%' '27%' '45%' '75%'})
ylabel('R2'), box off; title('Baseline R2');

subplot(1,3,2)
bar([-r2_coll(:,1) -r2_coll(:,2) -r2_coll(:,3) -r2_coll(:,4) -r2_coll(:,5)]')
set(gca, 'xtickLabel', {'10%' '16%' '27%' '45%' '75%'})
ylabel('R2'), box off; title('Collinear R2');

subplot(1,3,3)
bar([-r2_orth(:,1) -r2_orth(:,2) -r2_orth(:,3) -r2_orth(:,4) -r2_orth(:,5)]')
set(gca, 'xtickLabel', {'10%' '16%' '27%' '45%' '75%'})
ylabel('R2'), legend(num2str([1:12]')); box off; title('Orthogonal R2');

% 1 = coll vs. orth 2 = coll vs. base 3 = orth vs. base
% [h(:,:), p, ci, stats] = ttest(coll(:,con), orth(:,con));
% [h(:,:), p, ci, stats] = ttest(coll(:,con), bases(:,con));
% [h(:,:), p, ci, stats] = ttest(bases(:,con), orth(:,con));


% figure('Color', [1 1 1]);
% hold all;
% 
% errorbar(TheData(runs).p.testContrasts, nanmean(squeeze(TotalSuppressionIndexColinear(1,:,:))), ...
%     nanstd(squeeze(TotalSuppressionIndexColinear(1,:,:)))/sqrt(numel(list)), 'ro-');
% errorbar(TheData(runs).p.testContrasts, nanmean(squeeze(TotalSuppressionIndexColinear(2,:,:))), ...
%     nanstd(squeeze(TotalSuppressionIndexColinear(2,:,:)))/sqrt(numel(list)), 'r:^');
% errorbar(TheData(runs).p.testContrasts, nanmean(squeeze(TotalSuppressionIndexOrthogonal(1,:,:))), ...
%     nanstd(squeeze(TotalSuppressionIndexOrthogonal(1,:,:)))/sqrt(numel(list)), 'bo-');
% errorbar(TheData(runs).p.testContrasts, nanmean(squeeze(TotalSuppressionIndexOrthogonal(2,:,:))), ...
%     nanstd(squeeze(TotalSuppressionIndexOrthogonal(2,:,:)))/sqrt(numel(list)), 'b:^');
% 
% legend({'Coll Surround Suppression' 'Coll vWM Suppression' ...
%     'Orth Surround Suppression' 'Orth vWM Suppression'})
% ylabel('Suppression Index')
% xlabel('Contrast (%)')
% set(gca, 'XTick', TheData(runs).p.testContrasts, 'XTickLabel', round(TheData(runs).p.testContrasts*100), 'XScale', 'log')
% xlim([0.08 1])
