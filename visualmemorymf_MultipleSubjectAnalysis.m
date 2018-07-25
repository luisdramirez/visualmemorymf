%% MULTIPLE SUBJECT ANALYSIS %% 
% This script is designed to run based off of data calculated in the single
% subject analysis. It is designed to take all the average data collected
% from each subject who has ran (and ideally completeld 4) visualmemorymf
% trials. It will compara this data amongst subjects to produce concise
% results.

%% SETUP %%
clear;
close all;
expDir = '/Users/juliaschwartz/Desktop/visualmemorymf';
dataDir = '/Users/juliaschwartz/Desktop/visualmemorymf/data_master';
cd(dataDir)
load('visualmemory_condition_order')
load('visualmemory_subjectsRan')  

cd(dataDir)
files = struct2cell(dir(dataDir))';

[numFiles, ~] = size(files);
possibleFileNames = cell(length(visualmemory_subjectsRan),1);
for i = 1:length(visualmemory_subjectsRan)
    filename = strcat('data_visualmemorymf_exp_',visualmemory_subjectsRan{i},'.mat');
    possibleFileNames{i,1} = filename;
end

plotVar = 1; %if equal to 0, doesnt plot
printVar = 1; %if equal to 0, doesnt print

%preallocate a cell that will load theData structures from each participant
%into one cell
master_subjectData = cell(length(possibleFileNames),2); %files by 2 columns

% if any of files.name = possibleFilesNames then load the file and put
% into a cell array
for currfilenum = 1:numFiles
    dataFile = files{currfilenum,1};
    for i = 1:length(possibleFileNames)
        if strcmp(dataFile,possibleFileNames{i,1}) == 1
            load(dataFile)
            %index into the subjectsRan and compare to p.subject to find
            %out who is who if needed
            master_subjectData{i,1} = theData;
            master_subjectData{i,2} = subject;
            fprintf('\n%s and %s are compatible. Loading',dataFile, possibleFileNames{i,1})
        else
            fprintf('\n%s and %s not compatible',dataFile, possibleFileNames{i,1})
        end
    end
end
[subjectsLong,~] = size(master_subjectData);
p = theData(1).p;
centerContrast = (10.^linspace(log10(p.minContrast),log10(p.maxContrast),p.numContrasts));

%% ANALYSIS %%

% Preallocate matrices for subject data %
master.avgContEstimation = nan(length(visualmemory_subjectsRan),master_subjectData{1,1}(1).p.numContrasts,3);
master.avgDiffLocation = nan(length(visualmemory_subjectsRan),master_subjectData{1,1}(1).p.numContrasts,3);
master.avgDiffContrast = nan(length(visualmemory_subjectsRan),master_subjectData{1,1}(1).p.numContrasts,3);

% Loop through the subject data and add to matrices
for subj = 1:subjectsLong
    currentSubject = master_subjectData{subj,2};%There shouldnt be any Nans here if all people have been run 4 times total. For testing purposes take out nans
    for contNum = 1:master_subjectData{1,1}(1).p.numContrasts
        
        subject = master_subjectData{subj,2};
         
        master.avgContEstimation(subj,contNum,1) = subject.meanEstContPerception(contNum);
        master.avgContEstimation(subj,contNum,2) = subject.meanEstContWorkingMemory(contNum);
        master.avgContEstimation(subj,contNum,3) = subject.meanEstContBaseline(contNum);
        
        master.avgDiffContrast(subj,contNum,1) = subject.meanDiffContPerception(contNum);
        master.avgDiffContrast(subj,contNum,2) = subject.meanDiffContWorkingMemory(contNum);
        master.avgDiffContrast(subj,contNum,3) = subject.meanDiffContBaseline(contNum);
        
        master.avgDiffLocation(subj,contNum,1) = subject.meanEstContPerception(contNum);
        master.avgDiffLocation(subj,contNum,2) = subject.meanEstContWorkingMemory(contNum);
        master.avgDiffLocation(subj,contNum,3) = subject.meanEstContBaseline(contNum);
    end
end
%cant take nans out because then the size of the 3D matrix changes
    

%Take mean down the column of each page of each field of "master" to find
%averages over subjects. Then plot these as done in single analysis.
%Final average over all subjects will be in the last row, subjectsLong+1
for i = 1:subjectsLong
    if isnan(master.avgContEstimation(i,:,1)) == 0
        if exist('notNanPer','var') == 0
            notNanPer = (i);
        else
            notNanPer = horzcat([notNanPer,i]);
        end
    end
    if isnan(master.avgContEstimation(i,:,2)) == 0
        if exist('notNanWM','var') == 0
            notNanWM = (i);
        else
            notNanWM = horzcat([notNanWM,i]);
        end
    end
    if isnan(master.avgContEstimation(i,:,3)) == 0
        if exist('notNanBL','var') == 0
            notNanBL = (i);
        else
            notNanBL = horzcat([notNanBL,i]);
        end
    end
 end
for i = 1:length(notNanPer)
   if exist('perceptionmat','var') == 0
        perceptionmat = master.avgContEstimation(notNanPer(i),:,1);
    else
    perceptionmat = vertcat(perceptionmat,master.avgContEstimation(notNanPer(i),:,1));
    end
end
    master.avgContEstimation(subjectsLong+1,:,1) = mean(perceptionmat,1);
for i = 1:length(notNanWM)
    if exist('workingmemmat','var') == 0
        workingmemmat = master.avgContEstimation(notNanWM(i),:,2);
    else
    workingmemmat = vertcat(workingmemmat,master.avgContEstimation(notNanWM(i),:,2));
    end
end
    master.avgContEstimation(subjectsLong+1,:,2) = mean(workingmemmat,1);
for i = 1:length(notNanBL)
   if exist('baselinemat','var') == 0
        baselinemat = master.avgContEstimation(notNanBL(i),:,3);
    else
    baselinemat = vertcat(baselinemat,master.avgContEstimation(notNanBL(i),:,3));
    end
end
    master.avgContEstimation(subjectsLong+1,:,3) = mean(baselinemat,1);
    master.avgDiffContrast(subjectsLong+1,:,:) = mean(master.avgDiffContrast,1);
    master.avgDiffLocation(subjectsLong+1,:,:) = mean(master.avgDiffLocation,1);

if plotVar ~= 0
    % Total Avg. Estimated Contrast
        figure(1)
        set(gcf, 'Name', sprintf('Collective Perceived (Estimated) Contrast Versus Center Contrast'));
        loglog(centerContrast,mean(baselinemat,1),'-o') %Baseline
        hold on
        loglog(centerContrast,mean(perceptionmat,1),'-o') %Perception
        hold on
        loglog(centerContrast,mean(workingmemmat,1),'-o') %Working Memory
        hold on 
        loglog([0.1 0.8],[0.1 0.8],'k--') %log scale line
        xlabel('Center Contrast')
        ylabel('Perceived Contrast')
        xticks([0.1 0.8]); yticks([0.1 0.8]);
        xticklabels({'10','80'});yticklabels({'10','80'});
        % Legend incorporating nans
        if sum(isnan(master.avgContEstimation(subjectsLong+1,:,1))) ~= 0 && sum(isnan(master.avgContEstimation(subjectsLong+1,:,2))) ~= 0
            legend('Baseline','Log Scale') %perception and working memory are nans
        elseif sum(isnan(master.avgContEstimation(subjectsLong+1,:,1))) ~= 0 && sum(isnan(master.avgContEstimation(subjectsLong+1,:,2))) == 0
            legend('Baseline','Working Memory','Log Scale') % perception is nan
        elseif sum(isnan(master.avgContEstimation(subjectsLong+1,:,1))) == 0 && sum(isnan(master.avgContEstimation(subjectsLong+1,:,2))) ~= 0
            legend('Baseline','Perception','Log Scale') % working memory is nan
        else
            legend('Baseline','Perception','Working Memory','Log Scale')
        end 
        hold off 
        [howmanypercep,~] = size(perceptionmat);
        [howmanywm,~] = size(workingmemmat);
        [howmanybl,~] = size(baselinemat);
    % Histogram Plotting
        figure(2)
        set(gcf,'Name',sprintf('Histograms of Estimated Contrast at Each Contrast/Condition'));
        for i = 1:theData(1).p.numContrasts
            subplot(3,theData(1).p.numContrasts,i)
            hist(perceptionmat(:,i),howmanypercep) %perception
            hold on
            line([master.avgContEstimation(subjectsLong+1,i,1) master.avgContEstimation(subjectsLong+1,i,1)],ylim,'Linewidth',1.75,'Color','g');
            ylabel('PERCEPTION')
            xlabel('Contrast Level')
            xlim([0 1])
            hold off
            subplot(3,theData(1).p.numContrasts,i+5)
            hist(workingmemmat(:,i),howmanywm) %Working mem
            hold on
            line([master.avgContEstimation(subjectsLong+1,i,2) master.avgContEstimation(subjectsLong+1,i,2)],ylim,'Linewidth',1.75,'Color','g');
            hold off
            ylabel('WORKING MEMORY')
            xlabel('Contrast Level')
            xlim([0 1 ])
            subplot(3,theData(1).p.numContrasts,i+10)
            hist(baselinemat(:,i),howmanybl) %baseline
            hold on
            line([master.avgContEstimation(subjectsLong+1,i,3) master.avgContEstimation(subjectsLong+1,i,3)],ylim,'Linewidth',1.75,'Color','g');
            hold off
            ylabel('BASELINE')
            xlabel('Contrast Level')
            xlim([0 1])
        end   
        
        for i = 1:subjectsLong
            figure(3)
            set(gcf, 'Name', sprintf('Perception: Estimated versus Center Contrast'));
            subplot(2,subjectsLong/2,i)
            loglog(centerContrast,master.avgContEstimation(i,:,1),'-o') % Perception
            hold on
            loglog([0.1 0.8],[0.1 0.8],'k--') %log scale line
            xlabel('Center Contrast')
            ylabel('Estimated Contrast')
            hold off
            figure(4)
            set(gcf, 'Name', sprintf('Working Memory: Estimated versus Center Contrast'));
            subplot(2,subjectsLong/2,i)
            loglog(centerContrast,master.avgContEstimation(i,:,2),'-o') % Working Memory
            hold on
            loglog([0.1 0.8],[0.1 0.8],'k--') %log scale line
            xlabel('Center Contrast')
            ylabel('Estimated Contrast')
            hold off
        end
end
 


% T test between the different contrast levels between perception, working
% memory, and baseline average responses.
% JUST DO A T TEST THAT CORRESPONS WITH THE CONDITION IT WAS ON  (THE
% COMPARABLE VALUES FROM EACH PERSON)

% PERCEPTION & WORKING MEMORY STATISTICAL SIGNIFICANCE:
% Will only run if the same number of working memory and perception trials
% have been ran. (same size arrays necessary)
if sum(howmanypercep == howmanywm) == 1
    % Preallocate arrays to contain ttest information
    h_P_WM = NaN(1,theData(1).p.numContrasts);
    p_P_WM = NaN(1,theData(1).p.numContrasts);
    % Loops through number of contrasts - P & WM
    for i = 1:theData(1).p.numContrasts
    [h_P_WM(i),p_P_WM(i)] = ttest(perceptionmat(1:howmanypercep,i),workingmemmat(1:howmanywm,i));
    end
end
    %Loops through number of contrasts - BL & P
        %%% ARRAYS MUST BE SAME SIZE FOR T TEST, BUT BASELINE WILL HAVE
        %%% TWICE AS MANY AS THE OTHER CONDITIONS
    %Loops through number of contrasts - BL & WM
    
    
 
% error bars in the contrast graph (?)
 
