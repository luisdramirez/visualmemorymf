%%% vmmf_subject_analysis_v2.m

%% Setup wd & load data
clear all; close all; clc;
scrsz = get(groot,'ScreenSize');
expDir=pwd;
dataDir='data_master';
logtransform = 0;
indPlots = 0; superPlots = 1;

subjects = [1 2 3 4];

allData = struct([]);

% Load in data
cd(dataDir)
for nSubj = 1:numel(subjects)
    fileName = ['data_vmmf_v2_' num2str(subjects(nSubj)) '.mat'];
    if exist(fileName,'file') ~= 0
        load(fileName);
    end
    allData{nSubj} = theData;
    clear theData
end
cd(expDir)

%% organization

targetContrasts = allData{1}(1).p.centerContrasts;
probeOffsets =  [0 round(10.^linspace(0,log10(allData{1}(1).p.probeLocWidth),allData{1}(1).p.numOffsetLoc-1))]';

runCount = 1;
for nSubj = 1:numel(subjects)
    for nRun = 1:length(allData{nSubj})
        probeOffsetTemp = allData{nSubj}(nRun).p.trialEvents(:,4) - allData{nSubj}(nRun).p.trialEvents(:,2);
        probeOffsetTemp(probeOffsetTemp > 180) = probeOffsetTemp(probeOffsetTemp > 180)-360;
        probeOffsetTemp(probeOffsetTemp < -180) = 360+probeOffsetTemp(probeOffsetTemp < -180);
        probeOffsetTemp = abs(probeOffsetTemp);
        for nCon = 1:length(targetContrasts)
            for nOffset = 1:length(probeOffsets)
                currIndx = allData{nSubj}(nRun).p.trialEvents(:,1) == targetContrasts(nCon) & probeOffsetTemp == probeOffsets(nOffset);
                estContrasts(:,nOffset,nCon,nRun,nSubj) = allData{nSubj}(nRun).data.EstimatedContrast(currIndx);
                subjectData.(['Subject_' num2str(subjects(nSubj))]).(['Offset_' num2str(nOffset)]).(['Contrast_' num2str(nCon)])(:,nRun)= estContrasts(:,nOffset,nCon,nRun,nSubj);
                superSubj.(['Offset_' num2str(nOffset)]).(['Contrast_' num2str(nCon)])(:,runCount)= estContrasts(:,nOffset,nCon,nRun,nSubj);
            end
        end
        runCount=runCount+1;
    end
end

%% fitting

global numContrasts
numContrasts =  numel(targetContrasts);

%%% Individual Subject Fitting
estMeans = nan(numel(targetContrasts),numel(probeOffsets),numel(subjects));
estWidths = nan(numel(targetContrasts),numel(probeOffsets),numel(subjects));
for nSubj = 1:length(subjects)
    for nOffset = 1:length(probeOffsets)
        for nCon = 1:length(targetContrasts)
            currData = subjectData.(['Subject_' num2str(subjects(nSubj))]).(['Offset_' num2str(nOffset)]).(['Contrast_' num2str(nCon)])(:);
            edges = 0:0.05:1;
            [N(nCon,:),x] = histcounts(currData,edges);
            N(nCon,:) = N(nCon,:)./max(N(nCon,:));
            xvalues =(x(1:end-1)+x(2:end))/2;
        end
        
        startValues = [(targetContrasts) (repmat(0.1, [1 numel(targetContrasts)]))];
        options = optimset('MaxFunEvals', 5000.*numel(startValues), 'MaxIter', 5000.*numel(startValues));
        
        [est_params_tmp, r2] = fminsearch('mygauss_allContrasts', startValues, options, N, xvalues);
        tmp = reshape(est_params_tmp, [numContrasts 2])';
        est_params = [tmp(1,:); ones(1,numContrasts); tmp(2,:)];
        
        estMeans(:,nOffset,nSubj) = est_params(1,:);
        estWidths(:,nOffset,nSubj) = est_params(3,:);
        
        if indPlots
            figure('Position',[50 50 scrsz(3) scrsz(3)/2], 'color', [1 1 1], 'name', ['Offset_' num2str(probeOffsets(nOffset))])
            for nCon = 1:numel(targetContrasts)
                x_fit = 0.01:0.01:1;
                y_est = est_params(2,nCon)*exp(-(x_fit-est_params(1,nCon)).^2/(2*(est_params(3,nCon)^2)));
                xvalues = (x(1:end-1)+x(2:end))/2;
                subplot(1,numContrasts,nCon)
                bar(xvalues, N(nCon,:));
                hold on
                plot(x_fit, y_est, 'r', 'LineWidth', 2)
                box off; title(['S' num2str(nSubj) ' | ' 'Offset=' num2str(probeOffsets(nOffset)) 'deg | ' num2str(round(100*targetContrasts(nCon))) '% Contrast']);
                xlabel('Contrast (%)'); ylabel('Normalized count of perceived contrast');
                axis square
            end
        end
    end
end

%%% Super Subject Fitting
figure('Position',[50 50 scrsz(3) scrsz(3)/2], 'color', [1 1 1])
super.estMeans = nan(numel(targetContrasts),numel(probeOffsets));
super.estWidths = nan(numel(targetContrasts),numel(probeOffsets));
plotCount = 1;
for nOffset = 1:length(probeOffsets)
    for nCon = 1:length(targetContrasts)
        currData = superSubj.(['Offset_' num2str(nOffset)]).(['Contrast_' num2str(nCon)])(:);
        edges = 0:0.05:1;
        [N(nCon,:),x] = histcounts(currData,edges);
        N(nCon,:) = N(nCon,:)./max(N(nCon,:));
        xvalues =(x(1:end-1)+x(2:end))/2;
    end
    
    startValues = [(targetContrasts) (repmat(0.1, [1 numel(targetContrasts)]))];
    options = optimset('MaxFunEvals', 5000.*numel(startValues), 'MaxIter', 5000.*numel(startValues));
    
    [est_params_tmp, r2] = fminsearch('mygauss_allContrasts', startValues, options, N, xvalues);
    tmp = reshape(est_params_tmp, [numContrasts 2])';
    est_params = [tmp(1,:); ones(1,numContrasts); tmp(2,:)];
    
    super.estMeans(:,nOffset) = est_params(1,:);
    super.estWidths(:,nOffset) = est_params(3,:);
    
        for nCon = 1:numel(targetContrasts)
            x_fit = 0.01:0.01:1;
            y_est = est_params(2,nCon)*exp(-(x_fit-est_params(1,nCon)).^2/(2*(est_params(3,nCon)^2)));
            xvalues = (x(1:end-1)+x(2:end))/2;
            subplot(numel(probeOffsets),numContrasts,plotCount)
            bar(xvalues, N(nCon,:));
            hold on
            plot(x_fit, y_est, 'r', 'LineWidth', 2)
            box off; title(['Offset=' num2str(probeOffsets(nOffset)) 'deg | ' num2str(round(100*targetContrasts(nCon))) '% Contrast']);
            axis square
            plotCount=plotCount+1;
        end      
end
xlabel('Contrast (%)'); ylabel('Normalized count of perceived contrast');

%% Summaries

tmpContrasts = round(100.*targetContrasts);
plotLabels.contrasts  = {num2str(tmpContrasts(1)) num2str(tmpContrasts(2)) ...
    num2str(tmpContrasts(3)) num2str(tmpContrasts(4))}; % strings for axes labels
plotLabels.offsets = {num2str(probeOffsets(1)) num2str(probeOffsets(2)) num2str(probeOffsets(3)) ...
    num2str(probeOffsets(4)) num2str(probeOffsets(5))}; % strings for legend

%%% Mean Contrast Estimate Fitting Bar Plot (Subject Average)
y = 100.*mean(estMeans,3);
err = 100.*std(estMeans,[],3)/sqrt(size(estMeans,3));

figure('name','Mean Contrast Estimate (Subj. Avg.)')
for i = 1:nOffset
    errorbar(targetContrasts*100 , y(:,i), err(:,i))
    set(gca, 'YScale', 'log', 'XScale', 'log')
    hold on
end
plot([16.55 75],[16.55 75],'k--')
set(gca, 'YScale', 'log', 'XScale', 'log')
xticklabels(16.5,27,45,75);
xlabel('Target Contrast (%)')
ylabel('Percieved Contrast (%)')
legend({0,1,4,13,45})
hold off

figure('name','Mean Contrast Estimate (Subject Average)')
bar(y)
box off, hold on, xlabel('% Contrast Level'), ylabel('% Contrast (Mean)')
xticklabels(plotLabels.contrasts);
ngroups = size(y, 1);
nbars = size(y, 2);

% Calculating the width for each bar group to plot error
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), err(:,i), 'k','LineStyle','none');
end
legend(plotLabels.offsets,'Location','NorthWest')
hold off

%%% Width of Contrast Estimate Fitting Bar Plot (Subject Average)
y = 100.*mean(estWidths,3);
err = 100.*std(estWidths,[],3)/sqrt(size(estWidths,3));

figure('name','Width Contrast Estimates (Subject Average)')
bar(y)
box off, hold on, xlabel('% Contrast Level'), ylabel('% Contrast (Width)')
xticklabels(plotLabels.contrasts)
ngroups = size(y, 1);
nbars = size(y, 2);
% Calculating the width for each bar group to plot error
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), err(:,i), 'k','LineStyle','none');
end
hold off

%%% Actual vs. Perceived Contrast Plot (Subject Average)
y1=100.*mean(estMeans(:,1,:),3);
y2=100.*mean(estMeans(:,5,:),3);
x=100.*targetContrasts;
figure('name','Actual v Percevied Contrast (Subject Average)'), loglog(x,y1), hold on, loglog(x,y2), loglog(x,x,'--k');
box off, legend({'0' '45'},'Location','NorthWest'), xlabel('Actual Contrast (%)'), ylabel('Perceived Contrast (%)')

%%% Mean Contrast Estimate Fitting Bar Plot (Super Average)
y=100.*super.estMeans;
figure('name','Mean Contrast Estimate (Super Average)')
bar(y)
box off, hold on, xlabel('% Contrast Level'), ylabel('% Contrast (Mean)')
xticklabels(plotLabels.contrasts);

%%% Width of Contrast Estimate Fitting Bar Plot (Super Average)
y=100.*super.estWidths;
figure('name','Width Contrast Estimates (Super Average)')
bar(y)
box off, hold on, xlabel('% Contrast Level'), ylabel('% Contrast (Width)')
xticklabels(plotLabels.contrasts)

%%% Actual vs. Perceived Contrast Plot (Super Average)
y1=100.*estMeans(:,1);
y2=100.*estMeans(:,5);
x=100.*targetContrasts;
figure('name','Actual v Percevied Contrast (Super Average)'), loglog(x,y1), hold on, loglog(x,y2), loglog(x,x,'--k');
box off, legend({'0' '45'},'Location','NorthWest'), xlabel('Actual Contrast (%)'), ylabel('Perceived Contrast (%)')
