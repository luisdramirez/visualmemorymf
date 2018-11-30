%% MULTIPLE SUBJECT ANALYSIS %% 
% This script is designed to run based off of data calculated in the single
% subject analysis. It is designed to take all the average data collected
% from each subject who has ran (and ideally completed 4) visualmemorymf
% trials. It will compare this data amongst subjects to produce concise
% results.

%% SETUP %%
clear;
% close all;

experiment = 'exp';
%expDir = '/Users/juliaschwartz/Desktop/visualmemorymf'; %Lab computer
expDir =pwd; %Laptop
%dataDir = '/Users/juliaschwartz/Desktop/visualmemorymf/data_master'; %Lab computer
dataDir = 'data_master'; %Laptop
files = struct2cell(dir(dataDir))';

cd(dataDir)
load('visualmemory_condition_order')
load('visualmemory_subjectsRan') 


% removeSubjs = {'005'};
% for nSubj = 1:numel(removeSubjs)
%     visualmemory_subjectsRan{find(strcmp(visualmemory_subjectsRan{1,:}, removeSubjs{nSubj})),:} = [];
% end

[numFiles, ~] = size(files);
possibleFileNames = cell(size(visualmemory_subjectsRan,2),1);
%preallocate a cell that will load theData structures from each participant
%into one cell
master_subjectData = cell(length(possibleFileNames),2); %files by 2 columns
for i = 1:length(visualmemory_subjectsRan)
    filename = strcat('analyzed_visualmemorymf_exp_',visualmemory_subjectsRan{1,i},'.mat');
    if exist(filename,'file') ~= 0
        load(filename); % Loads HC, test, and Regular trials
        runNumbers = 1:length(theData);
        [fields, numRuns] = size(theData);
    end
    possibleFileNames{i,1} = filename;
    master_subjectData{i,1} = theData;
    master_subjectData{i,2} = subject;
end

plotVar = 1; %if equal to 0, doesnt plot
printVar = 1; %if equal to 0, doesnt print

cd(expDir)
fprintf('\n\n')
[subjectsLong,~] = size(master_subjectData);
p = theData(1).p;
centerContrast = (10.^linspace(log10(p.minContrast),log10(p.maxContrast),p.numContrasts));

%% ANALYSIS %%
perceptionIndex = 1;
workingmemIndex = 2;
baselineIndex = 3;

% Preallocate matrices for subject data %
master.avgContEstimation = nan(length(visualmemory_subjectsRan),master_subjectData{1,1}(1).p.numContrasts,3);
master.avgDiffLocation = nan(length(visualmemory_subjectsRan),master_subjectData{1,1}(1).p.numContrasts,3);
master.avgDiffContrast = nan(length(visualmemory_subjectsRan),master_subjectData{1,1}(1).p.numContrasts,3);

% Loop through the subject data and add to matrices
for subj = 1:subjectsLong
    currentSubject = master_subjectData{subj,2};%There shouldnt be any Nans here if all people have been run 4 times total. For testing purposes take out nans
    for contNum = 1:master_subjectData{1,1}(1).p.numContrasts
        
        subject = master_subjectData{subj,2};
         
        master.avgContEstimation(subj,contNum,perceptionIndex) = subject.meanEstContPerception(contNum);
        master.avgContEstimation(subj,contNum,workingmemIndex) = subject.meanEstContWorkingMemory(contNum);
        master.avgContEstimation(subj,contNum,baselineIndex) = subject.meanEstContBaseline(contNum);
        
        master.avgDiffContrast(subj,contNum,perceptionIndex) = subject.meanDiffContPerception(contNum);
        master.avgDiffContrast(subj,contNum,workingmemIndex) = subject.meanDiffContWorkingMemory(contNum);
        master.avgDiffContrast(subj,contNum,baselineIndex) = subject.meanDiffContBaseline(contNum);
        
        master.avgDiffLocation(subj,contNum,perceptionIndex) = subject.meanEstContPerception(contNum);
        master.avgDiffLocation(subj,contNum,workingmemIndex) = subject.meanEstContWorkingMemory(contNum);
        master.avgDiffLocation(subj,contNum,baselineIndex) = subject.meanEstContBaseline(contNum);
    end
end

%Take mean down the column of each page of each field of "master" to find
%averages over subjects. Then plot these as done in single analysis.
%Final average over all subjects will be in the last row, subjectsLong+1
for i = 1:subjectsLong
    if isnan(master.avgContEstimation(i,:,perceptionIndex)) == 0
        if exist('notNanPer','var') == 0
            notNanPer = (i);
        else
            notNanPer = horzcat([notNanPer,i]);
        end
    end
    if isnan(master.avgContEstimation(i,:,workingmemIndex)) == 0
        if exist('notNanWM','var') == 0
            notNanWM = (i);
        else
            notNanWM = horzcat([notNanWM,i]);
        end
    end
    if isnan(master.avgContEstimation(i,:,baselineIndex)) == 0
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
    perceptionmat = vertcat(perceptionmat,master.avgContEstimation(notNanPer(i),:,perceptionIndex));
    end
end
    master.avgContEstimation(subjectsLong+1,:,1) = mean(perceptionmat,1);
for i = 1:length(notNanWM)
    if exist('workingmemmat','var') == 0
        workingmemmat = master.avgContEstimation(notNanWM(i),:,workingmemIndex);
    else
    workingmemmat = vertcat(workingmemmat,master.avgContEstimation(notNanWM(i),:,workingmemIndex));
    end
end
    master.avgContEstimation(subjectsLong+1,:,2) = mean(workingmemmat,1);
for i = 1:length(notNanBL)
   if exist('baselinemat','var') == 0
        baselinemat = master.avgContEstimation(notNanBL(i),:,baselineIndex);
    else
    baselinemat = vertcat(baselinemat,master.avgContEstimation(notNanBL(i),:,baselineIndex));
    end
end
    master.avgContEstimation(subjectsLong+1,:,3) = mean(baselinemat,1);
    master.avgDiffContrast(subjectsLong+1,:,:) = mean(master.avgDiffContrast,1);
    master.avgDiffLocation(subjectsLong+1,:,:) = mean(master.avgDiffLocation,1);

if plotVar ~= 0
    % Total Avg. Estimated Contrast
        figure(1)
        set(gcf, 'Name', sprintf('Collective Perceived (Estimated) Contrast Versus Center Contrast'));
        % Baseline
        howmanybl = length(notNanBL);
        berror = (std(master.avgContEstimation(1:subjectsLong-1,:,3),'omitnan')/sqrt(howmanybl));
        bly = mean(baselinemat,1);
        errorbar(centerContrast,bly,berror,'LineWidth',1.25);   
        hold on
        
        % Perception
        hold on
        howmanyp = length(notNanPer);
        perror = (std(master.avgContEstimation(1:subjectsLong-1,:,1),'omitnan')/sqrt(howmanyp));
        py = mean(perceptionmat,1);
        errorbar(centerContrast,py,perror,'LineWidth',1.25);   
        hold on
        
        % Working Memory
        howmanywm = length(notNanWM);
        wmerror = (std(master.avgContEstimation(1:subjectsLong-1,:,2),'omitnan')/sqrt(howmanywm));
        wmy = mean(workingmemmat,1);
        errorbar(centerContrast,wmy,wmerror,'LineWidth',1.25);   
        hold on
        set(gca,'YScale','log','XScale','log');
        
        %log scale line
        hold on 
        loglog([0.1 0.8],[0.1 0.8],'k--') 
        xlabel('Center Contrast')
        ylabel('Perceived Contrast')
        xticks([0.1 0.8]); yticks([0.1 0.8]);
        xticklabels({'10','80'});yticklabels({'10','80'});
        % Legend incorporating nans
        if sum(isnan(master.avgContEstimation(subjectsLong+1,:,1))) ~= 0 && sum(isnan(master.avgContEstimation(subjectsLong+1,:,2))) ~= 0
            legend({'Baseline','Log Scale'},'Location','northwest') %perception and working memory are nans
        elseif sum(isnan(master.avgContEstimation(subjectsLong+1,:,1))) ~= 0 && sum(isnan(master.avgContEstimation(subjectsLong+1,:,2))) == 0
            legend({'Baseline','Working Memory','Log Scale'},'Location','northwest') % perception is nan
        elseif sum(isnan(master.avgContEstimation(subjectsLong+1,:,1))) == 0 && sum(isnan(master.avgContEstimation(subjectsLong+1,:,2))) ~= 0
            legend({'Baseline','Perception','Log Scale'},'Location','northwest') % working memory is nan
        else
            legend('Baseline','Perception','Working Memory','Log Scale')
        end 
        hold off 
        [howmanypercep,~] = size(perceptionmat);
        [howmanywm,~] = size(workingmemmat);
        [howmanybl,~] = size(baselinemat);
        
    %Bar graph plot (replace hist)    
        figure(2)
        set(gcf, 'Name', sprintf('Avg. Estimation per Contrast & Standard Error'));
        for i = 1:theData(1).p.numContrasts
            subplot(2,3,i)
            bar(3,mean(baselinemat(:,i)))
            hold on
            bar(1,mean(perceptionmat(:,i)))
            hold on
            bar(2,mean(workingmemmat(:,i)))
            xticks([1 2 3])
            xticklabels({'Perc','vWM','BL'})
            title(sprintf('Contrast of %.3f',centerContrast(i)))
            xlabel('Condition')
            ylabel('Contrast')
            hold on
            line([0 4],[centerContrast(i) centerContrast(i)],'Color','black','LineStyle','--')
            ylim([0 1])
            xlim([0.5 3.5])
            hold on
            for j = 1:howmanypercep
                plot(1,perceptionmat(j,i),'ko')
            end
            hold on
            for j = 1:howmanywm
                plot(2,workingmemmat(j,i),'ko')
            end
            hold on
            for j = 1:howmanybl
                plot(3,baselinemat(j,i),'ko')
            end
            hold on
            %Standard Error Bars
            line([3.05 3.05],[mean(baselinemat(:,i))-berror(i) mean(baselinemat(:,i))+berror(i)],'LineWidth',2,'Color','black')
            line([1.05 1.05],[mean(perceptionmat(:,i))-perror(i) mean(perceptionmat(:,i))+perror(i)],'LineWidth',2,'Color','black')
            line([2.05 2.05],[mean(workingmemmat(:,i))-wmerror(i) mean(workingmemmat(:,i))+wmerror(i)],'LineWidth',2,'Color','black')
            %ylim
            ylim([0 (max([max(perceptionmat(:,i)) max(baselinemat(:,i)) max(workingmemmat(:,i))]))+0.05])
        end
        
        for i = 1:subjectsLong
           % change if there is an even number of participants and fix
           % later on
           subjectsColNumber = round(subjectsLong);
%             evens = [2 4 6 8 10];
%             odds = [ 1 3 5 7 9];
%             subjectsColumn = 0;
%             for i = length(evens)
%                 if isequal(subjectsLong, odds(i)) == 1
%                     subjectsColumn = subjectsLong + 1;
%                     if subjectsColumn == 1
%                         subjectsColNumber = subjectsColumn;
%                     end
%                 elseif isequal(subjectsLong, evens(i)) == 1
%                     subjectsColNumber = subjectsLong ;
%                     if subjectsColumn == 1
%                       subjectsColNumber = subjectsColumn;
%                     end
%                 end
%             end
                   
            figure(3)
            set(gcf, 'Name', sprintf('Perception: Estimated versus Center Contrast'));
            
            subplot(2,subjectsColNumber/2,i)
            [howmanypercep,~] = size(master_subjectData{i,2}.perceptionmat);
            if howmanypercep == 1
                loglog(centerContrast,master_subjectData{i,2}.perceptionmat,'-o','LineWidth',1.25)
                hold on
            else
                perror = (std(master_subjectData{i,2}.perceptionmat)/sqrt(howmanypercep));
                py = mean(master_subjectData{i,2}.perceptionmat,1);
                errorbar(centerContrast,py,perror,'LineWidth',1.25,'Color','red')
            end
            overallData.perror(i,:) = perror;
            hold on
            loglog([0.1 0.8],[0.1 0.8],'k--') %log scale line
            set(gca,'YScale','log','XScale','log')
            xlim([0.1 0.8])
            ylim([0.1 0.8])
            ylabel('Estimated Contrast')
            set(gca,'YScale','log','XScale','log')
            title(sprintf('Subject %i',i))
            hold off
        end
        
        for i = 1:subjectsLong
            figure(4)
            set(gcf, 'Name', sprintf('Working Memory: Estimated versus Center Contrast'));
            subplot(2,subjectsColNumber/2,i)
           
            [howmanywm,~] = size(master_subjectData{i,2}.workingmemmat);
            if howmanywm == 1
                loglog(centerContrast,master_subjectData{i,2}.workingmemmat,'-o','LineWidth',1.25,'Color','red')
                hold on
            else
                wmerror = (std(master_subjectData{i,2}.workingmemmat)/sqrt(howmanywm));
                wmy = mean(master_subjectData{i,2}.workingmemmat,1);
                errorbar(centerContrast,wmy,wmerror,'LineWidth',1.25,'Color','red')
            end
            overallData.wmerror(i,:) = wmerror;
            hold on
            loglog([0.1 0.8],[0.1 0.8],'k--') %log scale line
            set(gca,'YScale','log','XScale','log')
            xlim([0 0.8])
            ylim([0 0.8])
            xlabel('Center Contrast')
            ylabel('Estimated Contrast')
            title(sprintf('Subject %i',i))
            hold off
        end
end

%Location Averages over 10 degree long bins
for subj = 1:length(master_subjectData)
    binAvg(:,subj) = master_subjectData{subj,2}.binAvgs(:,1);
    binContAvg(:,subj) = master_subjectData{subj,2}.binContAvgs(:,1);
end
binMeans = mean(binAvg,2,'omitnan');
binContMeans = mean(binContAvg,2,'omitNan');
if plotVar ~= 0
    figure('Color',[1 1 1])
    set(gcf,'Name','Average Location Difference Away Per 10 Degree Bin')
    plot(1:length(binMeans),binMeans','r')
    xlabel('Location on 360 Degree Circle')
    ylabel('Difference in Estimated versus Actual Degree Location (Abs)')
    xticks([0 9 18 27 36]);
    xticklabels({'East','North','West','South','East'})
end
binMeans(:,2) = (1:length(binMeans));  
binContMeans(:,2) = (1:length(binContMeans));
sortedBinMeans = sortrows(binMeans,1);
sortedBinContMeans = sortrows(binContMeans,1);

X = repmat(1,1,36);
labels = cell(1,length(binMeans));
for i = 1:length(binMeans)
    labels{1,i} = num2str(i*10);
end
numBins = length(X);
fig = figure;
ax = axes('Parent',fig);
locPie = pie(ax,ones(1,numBins),labels);
set(gcf,'name','Avg Location Difference Per 10 Degree Bin (yellow lowest) (red highest)')
midcol = fliplr((0:0.0278:1))';
autumnmat = zeros(36,3);
autumnmat(:,1) = 1;
autumnmat(:,2) = midcol;
autumnmat(:,4:5) = binMeans;
autumnmat = sortrows(autumnmat,4);
% sortedBinMeans col 2 is the order for the autumnmat
colorOrder = sortedBinMeans(:,2)'; %lowest first, so yellow, highest last, reed
% rgbmatrix = [ 1+(X(:) < 0).*X(:), 1-(X(:) > 0).*X(:),1-abs(X(:))];
for i = 1:length(binMeans)
    pieColorMap = autumnmat(i,1:3);
    set(locPie(i*2-1),'FaceColor',pieColorMap);
end
camroll(-90);

fig1 = figure;
ax1 = axes('Parent',fig1);
contPie = pie(ax1,ones(1,numBins),labels);
set(gcf,'name','Avg Contrast Error Per 10 Degree Bin (yellow lowest) (red highest)')
autumnmat = zeros(36,3);
autumnmat1(:,1) = repmat(1,36,1);
autumnmat1(:,2) = midcol;
autumnmat1(:,4:5) = binContMeans;
autumnmat1 = sortrows(autumnmat1,4);
for i = 1:length(binContMeans)
    pieColorMap1 = autumnmat1(i,1:3);
    set(contPie(i*2-1),'FaceColor',pieColorMap1);
end
camroll(-90);

perceptionLocDiffMean = zeros(subjectsLong,length(centerContrast));
workingmemLocDiffMean = zeros(subjectsLong,length(centerContrast));
baselineLocDiffMean = zeros(subjectsLong,length(centerContrast));
for subj = 1:length(master_subjectData)
    perceptionLocDiffMean(subj,:) = master_subjectData{subj,2}.perceptionLocDiffMean;
    workingmemLocDiffMean(subj,:) = master_subjectData{subj,2}.workingmemLocDiffMean;
    baselineLocDiffMean(subj,:) = master_subjectData{subj,2}.baselineLocDiffMean;
end
percepLocError = std(perceptionLocDiffMean)/sqrt(length(master_subjectData));
wmLocError = std(workingmemLocDiffMean)/sqrt(length(master_subjectData));
blLocError = std(baselineLocDiffMean)/sqrt(length(master_subjectData));


%% COMPARING DIFFERENT REPORT ORDERS %%
% report orders: index into visualmemory_subjectsRan, second row of a
% subjects column in the cell.
    %the first 4 subjects were asked in condition a: location, then
    %contrast. Also subject 10.
    %the next five subjects were asked in condition b: contrast, then
    %location.
 % we want to compare between the different reporting orders.
 for subj = 1:size(visualmemory_subjectsRan,2)
     firstOrder = visualmemory_subjectsRan{2,subj};
     secondOrder = visualmemory_subjectsRan{3,subj}; %save first and second orders into their subject struct
     master_subjectData{subj, 2}.firstOrder = firstOrder;
     if master_subjectData{1, 2}.numRuns  > 4
        master_subjectData{subj, 2}.secondOrder = secondOrder;
     end
     %THIS IS WHERE YOU WANT TO ADD THE OTHER ORDERS DATA TO THE MATRIX
     %if loop that will add specific subjects into a seperate master
     %subject dad matrix based off of their order number. then from there
     %you can split up the data.
     if firstOrder == 'a'
         subjectDataA{subj,1} = master_subjectData(subj,1); %should be theData in column 1, subject in column 2.
         subjectDataA{subj,2} = master_subjectData(subj,2);
     elseif firstOrder == 'b'
         subjectDataB{subj,1} = master_subjectData(subj,1); %should be theData in column 1, subject in column 2.
         subjectDataB{subj,2} = master_subjectData(subj,2);
     else
         error('subjects ran file doesnt contain order numbers')
     end
     
     if size(master_subjectData{subj,1},2) > 4
         if secondOrder == 'a'
             subjectDataA{subj,1} = master_subjectData(subj,1); % for SECOND order, should be theData in column 1, subject in column 2.
             subjectDataA{subj,2} = master_subjectData(subj,2);
         elseif secondOrder == 'b'
             subjectDataB{subj,1} = master_subjectData(subj,1); % for SECOND order,should be theData in column 1, subject in column 2.
             subjectDataB{subj,2} = master_subjectData(subj,2);
         else
             error('subjects ran file doesnt contain order numbers')
         end
     end
 end
     

%% A & B order comparison

for subj = 1:size(subjectDataA,1)
    if isempty(subjectDataA{subj,2}) == 0
        currentSubj = subjectDataA{subj,2};
            if currentSubj{1,1}.firstOrder == 'a'
                range = 1:4;
            elseif currentSubj{1,1}.secondOrder == 'a'
                range = 5:currentSubj{1, 1}.numRuns;
            end
        %Contrast Data
        meanContAPerception(subj,:) = nanmean(currentSubj{1,1}.avgEstContrast(range,:,1),1); %first page of avg contrast
        meanContAWM(subj,:) = nanmean(currentSubj{1,1}.avgEstContrast(range,:,2),1); %second page
        meanContABL(subj,:) = nanmean(currentSubj{1,1}.avgEstContrast(range,:,3),1); % third page
 
        %Location Data
        meanLocDiffAPerception(subj,:) = nanmean(currentSubj{1,1}.avgDiffLoc(range,:,1),1); %first page of avg contrast
        meanLocDiffAWM(subj,:) = nanmean(currentSubj{1,1}.avgDiffLoc(range,:,2),1); %second page
        meanLocDiffABL(subj,:) = nanmean(currentSubj{1,1}.avgDiffLoc(range,:,3),1); % third page
    end
 end

for subj = 1:size(subjectDataB,1)
    if isempty(subjectDataB{subj,2}) == 0
        currentSubj = subjectDataB{subj,2};
            if currentSubj{1,1}.firstOrder == 'b'
                range = 1:4;
            elseif currentSubj{1,1}.secondOrder == 'b'
                range = 5:currentSubj{1, 1}.numRuns;
            end
        %Contrast Data
        meanContBPerception(subj,:) = nanmean(currentSubj{1,1}.avgEstContrast(range,:,1),1); %first page of avg contrast
        meanContBWM(subj,:) = nanmean(currentSubj{1,1}.avgEstContrast(range,:,2),1); %second page
        meanContBBL(subj,:) = nanmean(currentSubj{1,1}.avgEstContrast(range,:,3),1); % third page
 
        %Location Data
        meanLocDiffBPerception(subj,:) = nanmean(currentSubj{1,1}.avgDiffLoc(range,:,1),1); %first page of avg contrast
        meanLocDiffBWM(subj,:) = nanmean(currentSubj{1,1}.avgDiffLoc(range,:,2),1); %second page
        meanLocDiffBBL(subj,:) = nanmean(currentSubj{1,1}.avgDiffLoc(range,:,3),1); % third page
    end
end

%Loop through an take 0s out of the matrices
rowswith0 = any(meanContAPerception==0,2);
meanContAPerception = meanContAPerception(~rowswith0, :);
clear rowswith0
rowswith0 = any(meanContAWM==0,2);
meanContAWM = meanContAWM(~rowswith0, :);
clear rowswith0
rowswith0 = any(meanContABL==0,2);
meanContABL = meanContABL(~rowswith0, :);
clear rowswith0
rowswith0 = any(meanLocDiffABL==0,2);
meanLocDiffABL = meanLocDiffABL(~rowswith0, :);
clear rowswith0
rowswith0 = any(meanLocDiffAPerception==0,2);
meanLocDiffAPerception = meanLocDiffAPerception(~rowswith0, :);
clear rowswith0
rowswith0 = any(meanLocDiffAWM==0,2);
meanLocDiffAWM = meanLocDiffAWM(~rowswith0, :);
clear rowswith0

rowswith0 = any(meanContBPerception==0,2);
meanContBPerception = meanContBPerception(~rowswith0, :);
clear rowswith0
rowswith0 = any(meanContBWM==0,2);
meanContBWM = meanContBWM(~rowswith0, :);
clear rowswith0
rowswith0 = any(meanContBBL==0,2);
meanContBBL = meanContBBL(~rowswith0, :);
clear rowswith0
rowswith0 = any(meanLocDiffBBL==0,2);
meanLocDiffBBL = meanLocDiffBBL(~rowswith0, :);
clear rowswith0
rowswith0 = any(meanLocDiffBPerception==0,2);
meanLocDiffBPerception = meanLocDiffBPerception(~rowswith0, :);
clear rowswith0
rowswith0 = any(meanLocDiffBWM==0,2);
meanLocDiffBWM = meanLocDiffBWM(~rowswith0, :);
clear rowswith0

  
%Means of each response
totalaverageContAPerception = nanmean(meanContAPerception);
totalaverageContAWM = nanmean(meanContAWM);
totalaverageContABL = nanmean(meanContABL);
totalavgPlocA = nanmean(meanLocDiffAPerception);
totalavgBLlocA = nanmean(meanLocDiffABL);
totalavgWMlocA = nanmean(meanLocDiffAWM);

totalaverageContBPerception = nanmean(meanContBPerception);
totalaverageContBWM = nanmean(meanContBWM);
totalaverageContBBL = nanmean(meanContBBL); 
totalavgPlocB = nanmean(meanLocDiffBPerception);
totalavgBLlocB = nanmean(meanLocDiffBBL);
totalavgWMlocB = nanmean(meanLocDiffBWM);



stdAP = (nanstd(meanContAPerception)/sqrt(size(meanLocDiffAPerception,1)));
stdAWM = (nanstd(meanContAWM)/sqrt(size(meanLocDiffAWM,1)));
stdABL = (nanstd(meanContABL)/sqrt(size(meanContABL,1)));

stdBP = (nanstd(meanContBPerception)/sqrt(size(meanContBPerception,1)));
stdBWM = (nanstd(meanContBWM)/sqrt(size(meanContBWM,1)));
stdBBL = (nanstd(meanContBBL)/sqrt(size(meanContBBL,1)));


% PLOT THE CENTER CONTRAST VERSUS ESTIMATED CONTRAST FOR BOTH ORDERING
% CONDITIONS
if plotVar ~= 0
    figure
    subplot(1,2,1) %A (location thwn contrast)
    errorbar(centerContrast,totalaverageContAPerception,stdAP,'LineWidth',2)
    hold on
    errorbar(centerContrast,totalaverageContAWM,stdAWM,'LineWidth',2)
    hold on
    errorbar(centerContrast,totalaverageContABL,stdABL,'LineWidth',2)
    hold on
    % change to include errorbars!
    loglog([0.1 0.8],[0.1 0.8],'k--') %log scale line
    set(gca,'YScale','log','XScale','log')
    xlim([0 0.8])
    ylim([0 0.8])
    xlabel('Center Contrast')
    ylabel('Estimated Contrast')
    title('Location then Contrast')
    legend('Perception','Working Memory','Baseline')
    
    subplot(1,2,2) %B, contrast then location
    errorbar(centerContrast,totalaverageContBPerception,stdBP,'LineWidth',2)
    hold on
    errorbar(centerContrast,totalaverageContBWM,stdBWM,'LineWidth',2)
    hold on
    errorbar(centerContrast,totalaverageContBBL,stdBBL,'LineWidth',2)
    hold on
    % change to include errorbars!
    loglog([0.1 0.8],[0.1 0.8],'k--') %log scale line
    set(gca,'YScale','log','XScale','log')
    xlim([0 0.8])
    ylim([0 0.8])
    xlabel('Center Contrast')
    ylabel('Estimated Contrast')
    title('Contrast then Location Trials')
    legend('Perception','Working Memory','Baseline')
    
    %Location difference between conditions
    figure
    subplot(3,2,1)
    %condition a right side, condition b is left side bl, percep, working memory
    bar(totalavgBLlocA,'r'); title('A: Baseline'); xlabel('Contrast Level'); ylabel('Location Difference'); ylim([0 17])
    subplot(3,2,2)
    bar(totalavgBLlocB,'r'); title('B: Baseline'); xlabel('Contrast Level'); ylabel('Location Difference'); ylim([0 17])
    subplot(3,2,3)
    bar(totalavgPlocA,'g'); title('A: Perception'); xlabel('Contrast Level'); ylabel('Location Difference'); ylim([0 17])
    subplot(3,2,4)
    bar(totalavgPlocB,'g'); title('B: Perception'); xlabel('Contrast Level'); ylabel('Location Difference'); ylim([0 17])
    subplot(3,2,5)
    bar(totalavgWMlocA,'b'); title('A: Working Memory'); xlabel('Contrast Level'); ylabel('Location Difference'); ylim([0 17])
    subplot(3,2,6)
    bar(totalavgWMlocB,'b'); title('B: Working Memory'); xlabel('Contrast Level'); ylabel('Location Difference'); ylim([0 17])
    
end
 
%%
figure('Color',[1 1 1])
    set(gcf,'name','Average Location Error/Bin, conditions split')
    locdiffmeanMat = zeros(3,length(centerContrast));
    locdiffmeanMat(1,:) = mean(perceptionLocDiffMean);
    locdiffmeanMat(2,:) = mean(workingmemLocDiffMean);
    locdiffmeanMat(3,:) = mean(baselineLocDiffMean);
    locerrorperception = std(locdiffmeanMat(1,:))/sqrt(length(locdiffmeanMat(1,:)));
    locerrorworkingmem = std(locdiffmeanMat(2,:))/sqrt(length(locdiffmeanMat(2,:)));
    locerrorbaseline = std(locdiffmeanMat(3,:))/sqrt(length(locdiffmeanMat(3,:)));
    locdiffmeanMat = locdiffmeanMat';
    locdiffmeanMatmean = mean(locdiffmeanMat);
 figure
    bar(1:3,locdiffmeanMatmean);
    xticks(fprintf('Perception','vWM','Baseline'));
    hold all
        errorbar(1,locdiffmeanMatmean(1),locerrorperception,'k-','LineWidth',2)
        errorbar(2,locdiffmeanMatmean(2),locerrorworkingmem,'k-','LineWidth',2)
        errorbar(3,locdiffmeanMatmean(3),locerrorbaseline,'k-','LineWidth',2)
    for i = 1:length(centerContrast)
        plot(1-0.05,locdiffmeanMat(i,1),'k.','MarkerSize',20)
        plot(2-0.05,locdiffmeanMat(i,2),'k.','MarkerSize',20)
        plot(3-0.05,locdiffmeanMat(i,3),'k.','MarkerSize',20)
    end

    xlabel('Contrast Level')
    ylabel('Difference in Location Estimate (Degrees)')
    
    
%%
%Subject Data A location difference graphs
for subj = 1:size(subjectDataA,1)
    if isempty(subjectDataA{subj,2}) == 0
        currentsubjectlocA = subjectDataA{subj,2};
        %Perception - order A
        if exist('locdiffAP','var') == 0
            locdiffAP = currentsubjectlocA{1,1}.perceptionLocDiffMean;
        else
            locdiffAP = [locdiffAP; currentsubjectlocA{1,1}.perceptionLocDiffMean];
        end
        %WM
        if exist('locdiffAWM','var') == 0
            locdiffAWM = currentsubjectlocA{1,1}.workingmemLocDiffMean;
        else
            locdiffAWM = [locdiffAWM; currentsubjectlocA{1,1}.workingmemLocDiffMean];
        end
        %baseline
        if exist('locdiffABL','var') == 0
            locdiffABL = currentsubjectlocA{1,1}.baselineLocDiffMean;
        else
            locdiffABL = [locdiffABL; currentsubjectlocA{1,1}.baselineLocDiffMean];
        end
    end
end
%Subject Data B location difference graphs
for subj = 1:size(subjectDataB,1)
    if isempty(subjectDataB{subj,2}) == 0
        currentsubjectlocB = subjectDataB{subj,2};
        %Perception - order B
        if exist('locdiffBP','var') == 0
            locdiffBP = currentsubjectlocB{1,1}.perceptionLocDiffMean;
        else
            locdiffBP = [locdiffBP; currentsubjectlocB{1,1}.perceptionLocDiffMean];
        end
        %WM
        if exist('locdiffBWM','var') == 0
            locdiffBWM = currentsubjectlocB{1,1}.workingmemLocDiffMean;
        else
            locdiffBWM = [locdiffBWM; currentsubjectlocB{1,1}.workingmemLocDiffMean];
        end
        %baseline
        if exist('locdiffBBL','var') == 0
            locdiffBBL = currentsubjectlocB{1,1}.baselineLocDiffMean;
        else
            locdiffBBL = [locdiffBBL; currentsubjectlocB{1,1}.baselineLocDiffMean];
        end
    end
end
MeanlocdiffAP = mean(locdiffAP);
MeanlocdiffAWM = mean(locdiffAWM);
MeanlocdiffABL = mean(locdiffABL);
MeanlocdiffBP = mean(locdiffBP);
MeanlocdiffBWM = mean(locdiffBWM);
MeanlocdiffBBL = mean(locdiffBBL);

stdAPloc = (std(locdiffAP)/sqrt(size(locdiffAP,1)));
stdAWMloc = (std(locdiffAWM)/sqrt(size(locdiffAWM,1)));
stdABLloc = (std(locdiffABL)/sqrt(size(locdiffABL,1)));
stdBPloc = (std(locdiffBP)/sqrt(size(locdiffBP,1)));
stdBWMloc = (std(locdiffBWM)/sqrt(size(locdiffBWM,1)));
stdBBLloc = (std(locdiffBBL)/sqrt(size(locdiffBBL,1)));

figure
subplot(1,2,1)
locdiffmeanMatA(1,:) = MeanlocdiffAP;
locdiffmeanMatA(2,:) = MeanlocdiffAWM;
locdiffmeanMatA(3,:) = MeanlocdiffABL;
locdiffmeanMatA = locdiffmeanMatA';
bar(1:5,locdiffmeanMatA);
hold all
title('Condition A: Location Estimation Error')
ylim([0 27])
for cont = 1:length(centerContrast)
        errorbar(cont,locdiffmeanMatA(cont,2),stdAWMloc(cont),'k')
        hold all 
        errorbar(cont-0.225,locdiffmeanMatA(cont,1),stdAPloc(cont),'k')
        errorbar(cont+0.225,locdiffmeanMatA(cont,3),stdABLloc(cont),'k')
        plot(cont-0.225,locdiffAP(:,cont),'ko')
        plot(cont,locdiffAWM(:,cont),'ko')
        plot(cont,locdiffABL(:,cont),'ko')
end
xticks([10 17 27 45 75])
xlabel('Contrast Level')
ylabel('Difference in Location Estimate (Degrees)')
legend({'Perception','Working Memory','Baseline'})
xticks([10 17 27 45 75]);


subplot(1,2,2)
locdiffmeanMatB(1,:) = MeanlocdiffBP;
locdiffmeanMatB(2,:) = MeanlocdiffBWM;
locdiffmeanMatB(3,:) = MeanlocdiffBBL;
locdiffmeanMatB = locdiffmeanMatB';
bar(1:5,locdiffmeanMatB);
hold all
title('Condition B: Location Estimation Error')
ylim([0 27])
for cont = 1:length(centerContrast)
        errorbar(cont,locdiffmeanMatB(cont,2),stdBWMloc(cont),'k-')
        errorbar(cont-0.225,locdiffmeanMatB(cont,1),stdBPloc(cont),'k-')
        errorbar(cont+0.225,locdiffmeanMatB(cont,3),stdBBLloc(cont),'k-')
        plot(cont-0.225,locdiffBP(:,cont),'ko')
        plot(cont,locdiffBWM(:,cont),'ko')
        plot(cont,locdiffBBL(:,cont),'ko')
end
xlabel('Contrast Level')
ylabel('Difference in Location Estimate (Degrees)')
legend({'Perception','Working Memory','Baseline'})
xticks([10 17 27 45 75]);
set(gcf,'Name','Location Error: Condition A vs. Condition B')


%% Statitsical Significance %%
% T test between the different contrast levels between perception, working
% memory, and baseline average responses.

% P & WM - ttest
if sum(howmanypercep == howmanywm) == 1
    % Preallocate arrays to contain ttest information
    h_P_WM = NaN(1,theData(1).p.numContrasts);
    p_P_WM = NaN(1,theData(1).p.numContrasts);
    % Loops through number of contrasts - P & WM
    for i = 1:theData(1).p.numContrasts
    [h_P_WM(i),p_P_WM(i)] = ttest(perceptionmat(1:howmanypercep,i),workingmemmat(1:howmanywm,i));
    end
end
  
% BL & WM - ttest
% Preallocate arrays to contain ttest information
    h_BL_WM = zeros(1,theData(1).p.numContrasts);
    p_BL_WM = zeros(1,theData(1).p.numContrasts);
    avgdBL_WM = zeros(subjectsLong,theData(1).p.numContrasts);
    fourArray = zeros(1,subjectsLong);
for i = 1:subjectsLong
    if master_subjectData{i,2}.numRuns == 4
        fourArray(i) = 1;
    end
end
for i = 1:subjectsLong
    if sum(fourArray == 1) == subjectsLong
        avgdBL_WM(i,:) = mean(master_subjectData{i,2}.baselinemat(master_subjectData{i,1}(1).p.trialSchedule  == 2,:));
    end
end
for i = 1:theData(1).p.numContrasts
    [h_BL_WM(i),p_BL_WM(i)] = ttest(avgdBL_WM(:,i)',workingmemmat(:,i)');
end

% BL & P - ttest
h_BL_P = zeros(1,theData(1).p.numContrasts);
p_BL_P = zeros(1,theData(1).p.numContrasts);
avgdBL_P = zeros(subjectsLong,theData(1).p.numContrasts);
for i = 1:subjectsLong
    avgdBL_P(i,:) = mean(master_subjectData{i,2}.baselinemat(master_subjectData{i,1}(1).p.trialSchedule  == 1,:));
end
for i = 1:theData(1).p.numContrasts
    [h_BL_P(i),p_BL_P(i)] = ttest(avgdBL_P(:,i),perceptionmat(:,i));
end

% Statistical Difference between how well people did in reporting contrast,
% versus how well people did in reporting location.

% % subjectContLocCompare 
% for subj = 1:size(visualmemory_subjectsRan,2)
%     if exist('subjectContCompare' ,'var') == 0
%         subjectContCompare = master_subjectData{subj,2}.avgLocationContrastCompareDiff(2,:);
%     else
%         subjectContCompare = [subjectContCompare; master_subjectData{subj,2}.avgLocationContrastCompareDiff(2,:)];
%     end
%     if exist('subjectLocCompare' ,'var') == 0
%         subjectLocCompare = master_subjectData{subj,2}.avgLocationContrastCompareDiff(1,:);
%     else
%         subjectLocCompare = [subjectLocCompare; master_subjectData{subj,2}.avgLocationContrastCompareDiff(1,:)];
%     end
%     meanPerPerson_LocComapre = mean(subjectLocCompare,2);
%     
   %statistical significant 
    
    
    
%     %0-100 scoring system:
%     if exist('comparisonScoreSystem' ,'var') == 0
%         comparisonScoreSystem = master_subjectData{subj,2}.avgLocationContrastCompareDiff(2,:);
%     else
%         comparisonScoreSystem = [subjectContCompare; master_subjectData{subj,2}.avgLocationContrastCompareDiff(2,:)];
%     end
%     
% end       
% come back to this point and seperate based on condition if necessary
% in order to compare: convert each person to a certain score for both
% contrast and location: 0 - 100 ?, and then add these scores together for
% a superscore between the two to assess how well people did on the task.
% figure; set(gcf,'Name','Comparing Location and Contrast Error: How well did a subject do?')
% for subj = 1:size(visualmemory_subjectsRan,2)
%     subplot(2,round(size(visualmemory_subjectsRan,2)/2),subj)
%     plot(subjectContCompare(subj,:),subjectLocCompare(subj,:),'.r','MarkerSize',20)
%     hold all
%     xlabel('Contrast Diff.');ylabel('Location Diff.'); title(sprintf('Subject %i',subj))
%     xlim([0 0.3]); ylim([0 30]);
% end


%% PRINTING %%
if printVar ~= 0
    for i = 1:theData(1).p.numContrasts
    fprintf('\n   BASELINE: At contrast %.3f, Estimated Contrast was %.4f',centerContrast(i),mean(baselinemat(:,i)))
    fprintf('\nWORKING MEM: At contrast %.3f, Estimated Contrast was %.4f',centerContrast(i),mean(workingmemmat(:,i)))
    fprintf('\n PERCEPTION: At contrast %.3f, Estimated Contrast was %.4f\n',centerContrast(i),mean(perceptionmat(:,i)))
    end
end

%% Saving variables %%
cd(dataDir)
overallData.baselineForWM = avgdBL_WM;
overallData.baselineForP = avgdBL_P;
overallData.perceptionmat = perceptionmat;
overallData.workingmemmat = workingmemmat;
overallData.baselineForWMMean = mean(avgdBL_WM);
overallData.baselineForPMean = mean(avgdBL_P);
overallData.perceptionmean = mean(perceptionmat);
overallData.workingmemmean = mean(workingmemmat);

overallData.master_subjectData = master_subjectData;
overallData.subjectDataA = subjectDataA;
overallData.meanContAPerception = meanContAPerception;
overallData.meanContAWM = meanContAWM;
overallData.meanContABL = meanContABL;
overallData.subjectDataB = subjectDataB;
overallData.meanContBPerception = meanContBPerception;
overallData.meanContBWM = meanContBWM;
overallData.meanContBBL = meanContBBL;
save(['data_visualmemorymf_overallData.mat'], 'overallData')
cd(expDir)