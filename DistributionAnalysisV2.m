clear
close all
% expDir = '/Users/ywatanabe/Dropbox/NormalizationVisualMemory/Experiment1';
expDir = '~/Dropbox/NormalizationVisualMemory/Experiment1';
cd([expDir '/Data'])

addpath(genpath(expDir));

plotOn = 1;
experiments = {'surrSuppression_2000_Delay*' 'vWM_surrSuppression_2000_Delay*'};
indvSubjectAnalysis = 0;

SubjectNumbers = [1:8 10:13];

numSubjects = length(SubjectNumbers);
numContrast = 5;
testContrast = 10.^linspace(log10(0.1),log10(0.75),numContrast);
numExperiments = numel(experiments);


if ~exist('perceivedContrast_distributionAnalysis.mat', 'file')
    % stimulus configurations x experimental condition
    allData = cell(3,2);
    subjData_perc = cell(3,numSubjects);
    subjData_vwm = cell(3,numSubjects);
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
            if plotOn == 1
                figure ('name', ['Exp. ' num2str(e) ' Subject ' num2str(subject)])
                hold all
            end
            
            %number of elements in this array
            numContrasts = numel(TheData(1).p.testContrasts);
            testContrast = TheData(1).p.testContrasts;
            
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
                
            end
            
            % append all observers data, to create one big observer
            allData{1,e} = [allData{1,e}; contrastcoll];
            allData{2,e} = [allData{2,e}; contrastorth];
            allData{3,e} = [allData{3,e}; contrastbase];
            
            if e == 1
                subjData_perc{1,subjCount} = contrastcoll;
                subjData_perc{2,subjCount} = contrastorth;
                subjData_perc{3,subjCount} = contrastbase;
            else
                subjData_vwm{1,subjCount} = contrastcoll;
                subjData_vwm{2,subjCount} = contrastorth;
                subjData_vwm{3,subjCount} = contrastbase;
            end
            clear contrastcoll contrastorth contrastbase
        end
        
    end
    
    % save out data
    save perceivedContrast_distributionAnalysis.mat allData subjData_perc subjData_vwm
else
    load('perceivedContrast_distributionAnalysis.mat')
end
cd(expDir);

%% Fit gaussian distribution to data of all observers
cd([expDir '/Figures/Distributions_fits'])
% options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000);

scrsz = get(groot,'ScreenSize');
figure(40), set(gcf, 'color', [1 1 1], 'Position',[50 50 scrsz(3) scrsz(3)/2]);
figure(41), set(gcf, 'color', [1 1 1]);

indx = 0:numContrast:3*numContrast;
x_fit = 0:0.01:1;
est_params_perc = cell(1,3);
est_params_vwm = cell(1,3);
% loop through the 3 surround conditions
for conf = 1:3
    
    for con = 1:numContrast
        %                 edges = 10.^(linspace(log10(0.05),log10(1), 21));
        edges = 0:0.05:1;
        [N_perc(con,:), x] = histcounts(allData{conf,1}(:,con), edges);
        N_perc(con,:) = N_perc(con,:)./max(N_perc(con,:));
        [N_vwm(con,:), ~] = histcounts(allData{conf,2}(:,con), edges);
        N_vwm(con,:) = N_vwm(con,:)./max(N_vwm(con,:));

        xvalues = (x(1:end-1)+x(2:end))/2;
        
    end
    
    % fit all subjects data simultaneous
    % mean, amplitude, sigma
    %             startValues = [testContrast repmat(median(max(N_perc,[],2)), [1 numContrast]) repmat(0.1, [1 numContrast])];
    startValues = [testContrast repmat(0.1, [1 numContrast])];
    options = optimset('MaxFunEvals', 5000.*numel(startValues), 'MaxIter', 5000.*numel(startValues));
    
    [est_params_tmp, r2_coll] = fminsearch('mygauss_allContrasts', startValues, options, N_perc, xvalues);
    tmp = reshape(est_params_tmp, [5 2])';
    est_params_perc{conf} = [tmp(1,:); ones(1,5); tmp(2,:)];
    
    % Get estimates for vwm when both mean and width are allowed to freely vary
    [est_params_tmp, r2_coll] = fminsearch('mygauss_allContrasts', startValues, options, N_vwm, xvalues);
    est_params_vwm_free{conf} = [est_params_tmp(1:5); ones(1,5); est_params_tmp(6:end)];
    
    % When the width is fixed by taking the estimates from the perception
    % condtion - does this explain the data as good as allow it to vary?
    startValues_nested = testContrast;
    fixedValues = {N_vwm, est_params_perc{conf}(3,:)};
    [est_params_tmp, r2_coll] = fminsearch('mygauss_allContrasts_nested', startValues_nested, options, fixedValues, xvalues);
    est_params_vwm_fixed{conf} = [est_params_tmp(1:5); ones(1,5); est_params_perc{conf}(3,:)];
    
    
    for con = 1:numContrast

        y_est_perc = est_params_perc{conf}(2,con)*exp(-(xvalues-est_params_perc{conf}(1,con)).^2/(2*(est_params_perc{conf}(3,con)^2)));
        y_est_vwm_free = est_params_vwm_free{conf}(2,con)*exp(-(xvalues-est_params_vwm_free{conf}(1,con)).^2/(2*((est_params_vwm_free{conf}(3,con))^2)));
        y_est_vwm_fixed = est_params_vwm_fixed{conf}(2,con)*exp(-(xvalues-est_params_vwm_fixed{conf}(1,con)).^2/(2*((est_params_vwm_fixed{conf}(3,con))^2)));
        
        allr2_perc(conf,con) = 1 - (sum((N_perc(con,:) - y_est_perc).^2) / sum((N_perc(con,:) - mean(N_perc(con,:))).^2));
        
        allr2_free_vwm(conf,con) = 1 - (sum((N_vwm(con,:) - y_est_vwm_free).^2) / sum((N_vwm(con,:) - mean(N_vwm(con,:))).^2));
        sse_free_vwm(con) = sum((N_vwm(con,:) - y_est_vwm_free).^2);
        
        allr2_fixed_vwm(conf,con) = 1 - (sum((N_vwm(con,:) - y_est_vwm_fixed).^2) / sum((N_vwm(con,:) - mean(N_vwm(con,:))).^2));
        sse_fixed_vwm(con) = sum((N_vwm(con,:) - y_est_vwm_fixed).^2);
        
        k = 1;
        prm = 1;
        n = length(N_vwm(con,:));
        v1 = k;
        v2 = n-(k+prm+1);
        
        Fratio_width(conf,con) = ((sse_fixed_vwm(con) - sse_free_vwm(con)) / v1) / (sse_free_vwm(con) / (v2));
       
        significance(conf,con) = 1-fcdf(Fratio_width(conf,con), v1, v2);
        p_corrected = pval_adjust(significance, 'bonferroni');
        
        % save individual figures
        x_fit = 0.01:0.01:1;
        y_est_perc = est_params_perc{conf}(2,con)*exp(-(x_fit-est_params_perc{conf}(1,con)).^2/(2*(est_params_perc{conf}(3,con)^2)));
        y_est_vwm_fixed = est_params_vwm_fixed{conf}(2,con)*exp(-(x_fit-est_params_vwm_fixed{conf}(1,con)).^2/(2*((est_params_vwm_fixed{conf}(3,con))^2)));
        y_est_vwm_free = est_params_vwm_free{conf}(2,con)*exp(-(x_fit-est_params_vwm_free{conf}(1,con)).^2/(2*((est_params_vwm_free{conf}(3,con))^2)));
        
        if conf == 1
            titleName = 'Collinear_contrast';
        elseif conf == 2
            titleName = 'Orthogonal_contrast';
        else
            titleName = 'NoSurround_contrast';
        end
        xvalues = (x(1:end-1)+x(2:end))/2;
        figure('Position',[50 50 scrsz(3) scrsz(3)/2], 'color', [1 1 1], 'name', [titleName num2str(con)])
        subplot(1,2,1)
        bar(xvalues, N_perc(con,:));
        hold all
        plot(x_fit, y_est_perc, 'r', 'LineWidth', 2)
        box off; title('Perception'); xlabel('contrast'); ylabel('Normalized count of perceived contrast')
        axis square
        subplot(1,2,2)
        bar(xvalues, N_vwm(con,:))
        hold all
        plot(x_fit, y_est_vwm_fixed, 'r', 'LineWidth', 2)
        plot(x_fit, y_est_vwm_free, 'm', 'LineWidth', 2)
        box off; title('Visual Memory')
        axis square
        
        
%         print(gcf, '-dpsc2', [titleName num2str(con) '.eps'])
        
        figure(40)
        subplot(3,5,indx(conf)+con)
        hold all
        plot(x_fit, y_est_perc, 'r');
        plot(x_fit, y_est_vwm_fixed, 'b');
        plot(x_fit, y_est_vwm_free, '.-k');
        plot(sqrt(x(1:end-1).*x(2:end)), N_perc(con,:), 'xr')
        plot(sqrt(x(1:end-1).*x(2:end)), N_vwm(con,:), 'xb')
        ylim([0 1]);
        
    end
    
    figure(41)
    subplot(2,3,conf)
    bar(1:5,[est_params_vwm_fixed{conf}(1,:)' est_params_vwm_free{conf}(1,:)' ])
    if conf == 1
        ylabel('Mean distributions')
    end
    box off
    
    subplot(2,3,3+conf)
    bar(1:5,[est_params_vwm_fixed{conf}(3,:)' (est_params_vwm_free{conf}(3,:))'])
    if conf == 1
        ylabel('Width distributions')
    end
    xlabel('Contrast levels');
    box off
   
    
end

% summary figures for all estimated parameters
figure('Color', [1 1 1], 'Name', 'Summary parameter estimates');
configurations = {'Collinear' 'Orthogonal' 'No surround'};
for n = 1:3
    % r2
    subplot(3,3,n)
    bar(1:5, [allr2_perc(n,:); allr2_free_vwm(n,:); allr2_fixed_vwm(n,:)]')
    set(gca, 'xtickLabel', {'10%' '17%' '27%' '45%' '75%'}); xlabel('Contrast (%)')
    box off; ylabel('r2'); title([configurations{n} ' r2'])
    axis square;

    subplot(3,3,n+3)
    bar(1:5, [est_params_perc{n}(1,:); est_params_vwm_free{n}(1,:); est_params_vwm_fixed{n}(1,:)]')
    set(gca, 'xtickLabel', {'10%' '17%' '27%' '45%' '75%'}); xlabel('Contrast (%)')
    box off; ylabel('Perceived contrast (%)'); title([configurations{n} ' mean']); 
    axis square;
    
    subplot(3,3,n+6)
    bar(1:5, [est_params_perc{n}(3,:); est_params_vwm_free{n}(3,:); est_params_vwm_fixed{n}(3,:)]')
    set(gca, 'xtickLabel', {'10%' '17%' '27%' '45%' '75%'}); xlabel('Contrast (%)')
    box off; ylabel('Perceived contrast (%)'); title([configurations{n} ' width']); 
    axis square;
    
end
subplot(3,3,4)
legend({'Simultaneous' 'Sequential_ free' 'Sequential_ fixed'}, 'Location', 'NorthWest'); 

