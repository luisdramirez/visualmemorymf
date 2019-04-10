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

targetContrasts = allData{1}.p.centerContrasts;

for nSubj = subjects
    %for nRun = 1:length(allData.TheData{nSubj}.data) % Number of runs
        for nCon = 1:length(targetContrasts) 
            conindx = allData{nSubj}.p.trialEvents(:,1) == targetContrasts(nCon);
            
            %Plots that graph a distribution for each contrast (Not
            %considering probe offset)
            figure(nSubj);
            set(gcf,'Name',' Subject %i')
            subplot(length(targetContrasts),1,nCon)
            hist(allData{nSubj}.data.EstimatedContrast(conindx),20);
            title(sprintf('Contrast Level: %.2f',targetContrasts(nCon)))
            allData{nSubj}.ContrastEstimates(nSubj,nCon) = mean(allData{nSubj}.data.EstimatedContrast(conindx));
            xlabel('Contrast');
            ylabel('Instances');

            estContrast(:,nCon,nSubj) = allData{nSubj}.data.EstimatedContrast(conindx);
            probeOffsetTemp = allData{1,nSubj}.p.trialEvents(find(conindx == 1),4) - allData{1,nSubj}.p.trialEvents(find(conindx == 1),2);
            probeOffsetTemp(probeOffsetTemp > 180) = probeOffsetTemp(probeOffsetTemp > 180)-360;
            probeOffsetTemp(probeOffsetTemp < -180) = 360+probeOffsetTemp(probeOffsetTemp < -180);
            probeOffset(:,nCon,nSubj) = probeOffsetTemp;

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
                 estConProbe = estContrast(indxOffsetLvl);
                    
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
                 meanAllSubjs(level,nCon,nSubj) = mean(estConProbe);
                 semAllSubjs(level,nCon,nSubj) = std(estConProbe)/sqrt(numel(estConProbe));
             end   
        end
    %end
end
