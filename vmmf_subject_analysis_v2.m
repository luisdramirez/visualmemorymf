%%% vmmf_subject_analysis_v2.m

%% Setup wd & load data
clear all; close all; clc;
scrsz = get(groot,'ScreenSize');
expDir=pwd;
dataDir='data_master';

subjects = 1;

allData = struct([]);

% Load in data
cd(dataDir)
for nSubj = subjects
    fileName = ['data_vmmf_v2_' num2str(nSubj) '.mat'];
    if exist(fileName,'file') ~= 0
        load(fileName);
    end
    allData{nSubj} = theData;
    clear theData
end
cd(expDir)
%% Data organization

targetContrasts = allData{1}(1).p.centerContrasts;
probeOffsets =  [0 round(10.^linspace(0,log10(allData{1}(1).p.probeLocWidth),allData{1}(1).p.numOffsetLoc-1))]';

for nSubj = subjects
    for nRun = 1:length(allData{nSubj}) % Number of runs
        probeOffsetTemp = allData{nSubj}(nRun).p.trialEvents(:,4) - allData{nSubj}(nRun).p.trialEvents(:,2);
        probeOffsetTemp(probeOffsetTemp > 180) = probeOffsetTemp(probeOffsetTemp > 180)-360;
        probeOffsetTemp(probeOffsetTemp < -180) = 360+probeOffsetTemp(probeOffsetTemp < -180);
        probeOffsetTemp = abs(probeOffsetTemp);
        
        for nCon = 1:length(targetContrasts)
            for nOffset = 1:length(probeOffsets)
                
                currIndx = allData{nSubj}(nRun).p.trialEvents(:,1) == targetContrasts(nCon) & probeOffsetTemp == probeOffsets(nOffset);
                estContrasts(:,nOffset,nCon,nRun,nSubj) = allData{nSubj}(nRun).data.EstimatedContrast(currIndx);
                subjectData.(['Subject_' num2str(subjects(nSubj))]).(['Offset_' num2str(nOffset)]).(['Contrast_' num2str(nCon)])(:,nRun)= estContrasts(:,nOffset,nCon,nRun,nSubj);
                
                %                 %store means in separate part
                %                 meanAllSubjs(nOffset,nCon,nRun,nSubj) = mean(estConProbe);
                %                 semAllSubjs(nOffset,nCon,nRun,nSubj) = std(estConProbe)/sqrt(numel(estConProbe));
                %                 %Plots that graph a distribution for each contrast (Not
                %                 %considering probe offset)
                %                 figure(nSubj); nbins = 20;
                %                 set(gcf,'Name',sprintf('Subject %i', nSubj))
                %                 subplot(length(targetContrasts),1,nCon)
                %                 hist(allData{nSubj}(nRun).data.EstimatedContrast(conindx),nbins);
                %                 title(sprintf('Contrast Level: %.2f',targetContrasts(nCon)))
                %                 allData{nSubj}(nRun).ContrastEstimates(nSubj,nCon) = mean(allData{nSubj}(nRun).data.EstimatedContrast(conindx));
                %                 xlabel('Contrast');
                %                 ylabel('Instances');
            end
        end
    end
end

%% fitting

for nSubj = 1:length(subjects)
    for nOffset = 1:length(probeOffsets)
        for nCon = 1:length(targetContrasts)
            currData = subjectData.(['Subject_' num2str(subjects(nSubj))]).(['Offset_' num2str(nOffset)]).(['Contrast_' num2str(nCon)])(:);
            edges = 0:0.05:1;
            [N(nCon,:),x] = histcounts(currData,edges);
            N(nCon,:) = N(nCon,:)./max(N(nCon,:));
            xvalues = (x(1:end-1)+x(2:end))/2;
        end
        
        global numContrasts
        numContrasts =  numel(targetContrasts);
        startValues = [targetContrasts repmat(0.1, [1 numel(targetContrasts)])];
        options = optimset('MaxFunEvals', 5000.*numel(startValues), 'MaxIter', 5000.*numel(startValues));
        
        [est_params_tmp, r2] = fminsearch('mygauss_allContrasts', startValues, options, N, xvalues);
        tmp = reshape(est_params_tmp, [numContrasts 2])';
        est_params = [tmp(1,:); ones(1,numContrasts); tmp(2,:)];
        
        figure('Position',[50 50 scrsz(3) scrsz(3)/2], 'color', [1 1 1], 'name', ['Offset_' num2str(probeOffsets(nOffset))])
        for nCon = 1:numel(targetContrasts)
            x_fit = 0.01:0.01:1;
            y_est = est_params(2,nCon)*exp(-(x_fit-est_params(1,nCon)).^2/(2*(est_params(3,nCon)^2)));
            xvalues = (x(1:end-1)+x(2:end))/2;
            subplot(1,numContrasts,nCon)
            bar(xvalues, N(nCon,:));
            hold on
            plot(x_fit, y_est, 'r', 'LineWidth', 2)
            box off; title(['Offset=' num2str(probeOffsets(nOffset)) ' | ' num2str(round(100*targetContrasts(nCon))) '% Contrast']);
            xlabel('Contrast (%)'); ylabel('Normalized count of perceived contrast');
            axis square
        end
    end
end

%% plotting

