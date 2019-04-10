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
        for nCon = 1:length(targetContrasts) 
            conindx = allData{nSubj}(nRun).p.trialEvents(:,1) == targetContrasts(nCon);
            
            %Plots that graph a distribution for each contrast (Not
            %considering probe offset)
            figure(nSubj); nbins = 20;
            set(gcf,'Name',sprintf('Subject %i', nSubj))
            subplot(length(targetContrasts),1,nCon)
            hist(allData{nSubj}(nRun).data.EstimatedContrast(conindx),nbins);
            title(sprintf('Contrast Level: %.2f',targetContrasts(nCon)))
            allData{nSubj}(nRun).ContrastEstimates(nSubj,nCon) = mean(allData{nSubj}(nRun).data.EstimatedContrast(conindx));
            xlabel('Contrast');
            ylabel('Instances');

            estContrast(:,nCon,nRun,nSubj) = allData{nSubj}(nRun).data.EstimatedContrast(conindx);
            probeOffsetTemp = allData{nSubj}(nRun).p.trialEvents(find(conindx == 1),4) - allData{nSubj}(nRun).p.trialEvents(find(conindx == 1),2);
            probeOffsetTemp(probeOffsetTemp > 180) = probeOffsetTemp(probeOffsetTemp > 180)-360;
            probeOffsetTemp(probeOffsetTemp < -180) = 360+probeOffsetTemp(probeOffsetTemp < -180);
            probeOffset(:,nCon,nRun,nSubj) = probeOffsetTemp;

        end

        offsetLevels = unique(abs(probeOffset));
%         figure;
%         set(gcf,'Name',['Subject_' num2str(subjects(nSubj))])
        for level = 1:length(offsetLevels)
             offset = offsetLevels(level);
             figure;
             set(gcf,'Name',['Subject ' num2str(subjects(nSubj)) ', Offset ' num2str(offset)])
             %loop through contrasts
             for nCon = 1:length(targetContrasts)
                 indxOffsetLvl = find(abs(probeOffset(:,nCon)) == offset);
                 estConProbe = estContrast(indxOffsetLvl,nCon);
                    
                 subplot(1,length(targetContrasts),nCon)
                 histogram(estConProbe,20)
                 title(['Contrast ' num2str(targetContrasts(nCon))])
                 hold all
                 xlim([0.1 0.75])
                 ylim([0 7])
                 
                 %store the entirety of the data in totalSubjectData
                 %totalSubjectData.(['Subject_' num2str(subjects(nSubj))]).(['Probe_' num2str(level)]).(['Contrast_' num2str(nCon)]) = [];
                 totalSubjectData.(['Subject_' num2str(subjects(nSubj))]).(['Probe_' num2str(level)]).(['Contrast_' num2str(nCon)]) = estConProbe;
                    
                 %store means in separate part
                 meanAllSubjs(level,nCon,nRun,nSubj) = mean(estConProbe);
                 semAllSubjs(level,nCon,nRun,nSubj) = std(estConProbe)/sqrt(numel(estConProbe));
             end   
        end
    end
end
