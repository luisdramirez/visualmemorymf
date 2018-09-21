%% MULTIPLE SUBJECT ANALYSIS %% 
% This script is designed to run based off of data calculated in the single
% subject analysis. It is designed to take all the average data collected
% from each subject who has ran (and ideally completeld 4) visualmemorymf
% trials. It will compara this data amongst subjects to produce concise
% results.

%% SETUP %%
clear;
close all;
%expDir = '/Users/juliaschwartz/Desktop/visualmemorymf'; %Lab computer
expDir = '/Users/julia/Desktop/Ling Lab/Experiments/visualmemorymf'; %Laptop
%dataDir = '/Users/juliaschwartz/Desktop/visualmemorymf/data_master'; %Lab computer
dataDir = '/Users/julia/Desktop/Ling Lab/Experiments/visualmemorymf/data_master'; %Laptop
cd(dataDir)
load('visualmemory_condition_order')
load('visualmemory_subjectsRan')  

cd(dataDir)
files = struct2cell(dir(dataDir))';

[numFiles, ~] = size(files);
possibleFileNames = cell(size(visualmemory_subjectsRan,2),1);
for i = 1:length(visualmemory_subjectsRan)
    filename = strcat('analyzed_visualmemorymf_exp_',visualmemory_subjectsRan{1,i},'.mat');
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
            fprintf('\nLoading %s',dataFile)
        end
    end
end
fprintf('\n\n')
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
           subjectsColNumber = subjectsLong + 1;
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
percepLocError = std(perceptionLocDiffMean)/length(master_subjectData);
wmLocError = std(workingmemLocDiffMean)/length(master_subjectData);
blLocError = std(baselineLocDiffMean)/length(master_subjectData);


%% COMPARING DIFFERENT REPORT ORDERS %%
% report orders: index into visualmemory_subjectsRan, second row of a
% subjects column in the cell.
    %the first 4 subjects were asked in condition a: location, then
    %contrast. Also subject 10.
    %the next five subjects were asked in condition b: contrast, then
    %location.
 % we want to compare between the different reporting orders.
 for subj = 1:size(visualmemory_subjectsRan,2)
     currentOrder = visualmemory_subjectsRan{2,subj};
     %if loop that will add specific subjects into a seperate master
     %subject dad matrix based off of their order number. then from there
     %you can split up the data.
     if currentOrder == 'a'
         subjectDataA{subj,1} = master_subjectData(subj,1); %should be theData in column 1, subject in column 2.
         subjectDataA{subj,2} = master_subjectData(subj,2);
     elseif currentOrder == 'b'
         subjectDataB{subj,1} = master_subjectData(subj,1); %should be theData in column 1, subject in column 2.
         subjectDataB{subj,2} = master_subjectData(subj,2);
     else
         error('subjects ran file doesnt contain order numbers')
     end
 end
     
% Use these indexes to see wherein the cells actually have elements, use it
% to loop through
%from here make an array that has the index numbers for the rows equal to
%one, then loop through that array
indexA = cellfun(@isempty, subjectDataA) == 0;
indexB = cellfun(@isempty, subjectDataB) == 0;

arrayforA = zeros(1,size(indexA,1));
for i = 1:size(indexA,2)
    if indexA(i,1) == 1 || indexA(i,2) == 1
        arrayforA(i) = i;
    end
end
arrayforA(arrayforA==0) = [];

arrayforB = zeros(1,size(indexB,1));
for i = 1:size(indexB,2) 
    if indexB(i,1) == 1 || indexA(i,2) == 1
        arrayforB(i) = i;
    end
end
arrayforB(arrayforB==0) = [];

%% move on from the array thing and come back!!! 
%for now just go through subject data a and b and if the cell row is empty
%skip that row.

for subj = 1:size(subjectDataA,1)
    if isempty(subjectDataA{subj,2}) == 0
        currentSubj = subjectDataA{subj,2};
        %Perception - order A
        if exist('meanContAPerception','var') == 0
            meanContAPerception = currentSubj{1,1}.meanEstContPerception;
        else
            meanContAPerception = [meanContAPerception; currentSubj{1,1}.meanEstContPerception];
        end
        % Working Memory - order A
        if exist('meanContAWM','var') == 0
            meanContAWM = currentSubj{1,1}.meanEstContWorkingMemory;
        else
            meanContAWM = [meanContAWM; currentSubj{1,1}.meanEstContWorkingMemory];
        end
        % Baseline - order A
        if exist('meanContABL','var') == 0
            meanContABL = currentSubj{1,1}.meanEstContBaseline;
        else
            meanContABL = [meanContABL; currentSubj{1,1}.meanEstContBaseline];
        end
    end
end
%Means of each response
totalaverageContAPerception = mean(meanContAPerception);
totalaverageContAWM = mean(meanContAWM);
totalaverageContABL = mean(meanContABL);

% Subject B
for subj = 1:size(subjectDataB,1)
    if isempty(subjectDataB{subj,2}) == 0
        currentsubjectB = subjectDataB{subj,2};
        %Perception - order B
        if exist('meanContBPerception','var') == 0
            meanContBPerception = currentsubjectB{1,1}.meanEstContPerception;
        else
            meanContBPerception = [meanContBPerception; currentsubjectB{1,1}.meanEstContPerception];
        end
        % Working Memory - order B
        if exist('meanContBWM','var') == 0
            meanContBWM = currentsubjectB{1,1}.meanEstContWorkingMemory;
        else
            meanContBWM = [meanContBWM; currentsubjectB{1,1}.meanEstContWorkingMemory];
        end
        % Baseline - order B
        if exist('meanContBBL','var') == 0
            meanContBBL = currentsubjectB{1,1}.meanEstContBaseline;
        else
            meanContBBL = [meanContBBL; currentsubjectB{1,1}.meanEstContBaseline];
        end
    end
end
       
totalaverageContBPerception = mean(meanContBPerception);
totalaverageContBWM = mean(meanContBWM);
totalaverageContBBL = mean(meanContBBL); 

stdAP = (std(meanContAPerception)/sqrt(size(meanContAPerception,1)));
stdAWM = (std(meanContAWM)/sqrt(size(meanContAWM,1)));
stdABL = (std(meanContABL)/sqrt(size(meanContABL,1)));

stdBP = (std(meanContBPerception)/sqrt(size(meanContBPerception,1)));
stdBWM = (std(meanContBWM)/sqrt(size(meanContBWM,1)));
stdBBL = (std(meanContBBL)/sqrt(size(meanContBBL,1)));


% PLOT THE CENTER CONTRAST VERSUS ESTIMATED CONTRAST FOR BOTH ORDERING
% CONDITIONS
if plotVar ~= 0
    figure
    subplot(1,2,1) %A (location then contrast)
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
    title('Location then Contrast Trials')
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
end
    
 
   
    
    
%     howmanywm = length(notNanWM);
%         wmerror = (std(master.avgContEstimation(1:subjectsLong-1,:,2),'omitnan')/sqrt(howmanywm));
%         wmy = mean(workingmemmat,1);
%         errorbar(centerContrast,wmy,wmerror,'LineWidth',1.25);

   
             
     %from here grab the data from each subject (within the for loop, ie
     %contrast differences, location differences (seperate by level of
     %actual contrast).
     %comparing the two, make graphs where it will display the two on the
     %same window but with different features.
     %statistical sigificance - ttest, and all of the rest.

 
%%
figure('Color',[1 1 1])
    set(gcf,'name','Average Location Error/Bin, conditions split')
    locdiffmeanMat = zeros(3,length(centerContrast));
    locdiffmeanMat(1,:) = mean(perceptionLocDiffMean);
    locdiffmeanMat(2,:) = mean(workingmemLocDiffMean);
    locdiffmeanMat(3,:) = mean(baselineLocDiffMean);
    locdiffmeanMat = locdiffmeanMat';
    bar(1:5,locdiffmeanMat);
    hold all
    for cont = 1:length(centerContrast)
        errorbar(cont,locdiffmeanMat(cont,2),wmLocError(cont),'k-')
        errorbar(cont-0.225,locdiffmeanMat(cont,1),percepLocError(cont),'k-')
        errorbar(cont+0.225,locdiffmeanMat(cont,3),blLocError(cont),'k-')
    end
    xlabel('Contrast Level')
    ylabel('Difference in Location Estimate (Degrees)')
    legend({'Perception','Working Memory','Baseline'})
    xticks([1 2 3 4 5]);
    

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
    if master_subjectData{i,2}.nTrials == 4
        fourArray(i) = 1;
    end
end
for i = 1:subjectsLong
    if sum(fourArray == 1) == 4
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
    

%% PRINTING %%
if printVar ~= 0
    for i = 1:theData(1).p.numContrasts
    fprintf('\n   BASELINE: At contrast %.3f, Estimated Contrast was %.4f',centerContrast(i),mean(baselinemat(:,i)))
    fprintf('\nWORKING MEM: At contrast %.3f, Estimated Contrast was %.4f',centerContrast(i),mean(workingmemmat(:,i)))
    fprintf('\n PERCEPTION: At contrast %.3f, Estimated Contrast was %.4f\n',centerContrast(i),mean(perceptionmat(:,i)))
    end
end

%% Saving variables %%
overallData.baselineForWM = avgdBL_WM;
overallData.baselineForP = avgdBL_P;
overallData.perceptionmat = perceptionmat;
overallData.workingmemmat = workingmemmat;
overallData.baselineForWMMean = mean(avgdBL_WM);
overallData.baselineForPMean = mean(avgdBL_P);
overallData.perceptionmean = mean(perceptionmat);
overallData.workingmemmean = mean(workingmemmat);
save(['data_visualmemorymf_overallData.mat'], 'overallData')