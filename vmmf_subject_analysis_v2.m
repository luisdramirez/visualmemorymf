%%% vmmf_subject_analysis_v2.m

%% Setup wd & load data
clear all; close all; clc;
scrsz = get(groot,'ScreenSize');
expDir=pwd;
dataDir='data_master';
subjects = 1;

data = struct('TheData',[],'ContrastEstimates',[],'ProbeOffset',[]);
data.ContrastEstimate = nan(length(subjects),5); % initialize average contrast estimates

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

cd(expDir)
%% Data organization

targetContrasts = data.TheData{1}.p.centerContrasts;

for nSubj = subjects
    for nCon = 1:length(targetContrasts) 
        conindx = data.TheData{nSubj}.p.trialEvents(:,1) == targetContrasts(nCon);
        data.ContrastEstimates(nSubj,nCon) = mean(data.TheData{nSubj}.data.EstimatedContrast(conindx));
%         probeOffset(:,nCon,nSubj) = data.TheData{nSubj}.p.trialEvents(:,4) - data.TheData{nSubj}.p.trialEvents(:,2);
%         probeOffset(probeOffset > 180) = probeOffset(probeOffset(:,nCon,nSubj) > 180)-360;
%         probeOffset(probeOffset < -180) = 360-probeOffset(probeOffset(:,nCon,nSubj) < -180); 
    end
end
