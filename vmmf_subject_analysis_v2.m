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

% for level = 1:length(offsetLevels)
%     
%     figure;
%     set(gcf,'Name',['Subject ' num2str(subjects(nSubj)) ', Offset ' num2str(currOffset)])
%     %loop through contrasts
%     for nCon = 1:length(targetContrasts)
% 
%         
%         subplot(1,length(targetContrasts),nCon)
%         histogram(estConProbe,20)
%         title(['Contrast ' num2str(targetContrasts(nCon))])
%         hold all
%         xlim([0.1 0.75])
%         ylim([0 7])
%         
% 
%     end
% end

%% plotting

