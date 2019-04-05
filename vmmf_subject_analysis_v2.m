%%% vmmf_subject_analysis_v2.m

%% Setup wd & load data
clear all; close all; clc;
scrsz = get(groot,'ScreenSize');
expDir=pwd;
dataDir='data_master';
subjects = 1;
subjPlots = 0;
groupPlots = 0;

data = struct('TheData',[],'ContrastEstimates',[],'ProbeOffset',[]);

% Load in data
cd(dataDir)

for nSubj = subjects
    fileName = ['data_vmmf_v2_' num2str(nSubj) '.mat'];
    if exist(fileName,'file') ~= 0
        load(fileName);
    end
    data.TheData{nSubj} = theData;
    clear theData
end
data.ContrastEstimates = nan(length(subjects),length(data.TheData{1, 1}.p.centerContrasts)); % initialize average contrast estimates

cd(expDir)
%% Data organization

targetContrasts = data.TheData{1}.p.centerContrasts;


for nSubj = subjects
    for runNum = 1:size(size(data.TheData{1, 1},1) % Number of runs
        for nCon = 1:length(targetContrasts) 
            conindx = data.TheData{nSubj}.p.trialEvents(:,1) == targetContrasts(nCon);
            if subjPlots 
                %Plots that graph a distribution for each contrast (Not
                %considering probe offset)
                figure(nSubj);
                set(gcf,'Name',' Subject %i')
                subplot(length(targetContrasts),1,nCon)
                hist(data.TheData{nSubj}.data.EstimatedContrast(conindx),20);
                title(sprintf('Contrast Level: %.2f',targetContrasts(nCon)))
                data.ContrastEstimates(nSubj,nCon) = mean(data.TheData{nSubj}.data.EstimatedContrast(conindx));
                xlabel('Contrast');
                ylabel('Instances');
            end
            estContrast(:,nCon,nSubj) = data.TheData{nSubj}.data.EstimatedContrast(conindx);
            probeOffsetTemp = data.TheData{nSubj}.p.trialEvents(find(conindx == 1),4) - data.TheData{nSubj}.p.trialEvents(find(conindx == 1),2);
            probeOffsetTemp(probeOffsetTemp > 180) = probeOffsetTemp(probeOffsetTemp > 180)-360;
            probeOffsetTemp(probeOffsetTemp < -180) = 360+probeOffsetTemp(probeOffsetTemp < -180);
            probeOffset(:,nCon,nSubj) = probeOffsetTemp;

        end

        offsetLevels = unique(abs(probeOffset));
        for level = 1:length(offsetLevels)
             offset = offsetLevels(level);
             %loop through contrasts
             for nCon = 1:length(targetContrasts)
                 indxOffsetLvl = find(abs(probeOffset(:,nCon)) == offset);
                 estConProbe = estContrast(indxOffsetLvl);
                 %store this in two seperate places now - the subject loop and
                 %the main matrix
                 totalSubjectData.(['Subject_' num2str(subjects(1),'%d')]) = [];


             end
                 %data.TheData{nSubj}.data.EstimatedContrast(conindx)
                 %probe ABs, find(probeAbs == 1) , estContrast(find)
        end
    end
    
end
