%% Analysis Designed to Compare statistical results between data sets of different runs on the same subject %%
clear all
close all
load('data_visualmemorymf_test_JS.mat') 

% Under assumption that a completed experiment on someone is 4 runs 
[~,trialsran] = size(theData);
if trialsran == 4
    run1Stats = theData(1).stats;
    run2Stats = theData(2).stats;
    run3Stats = theData(3).stats;
    run4Stats = theData(4).stats;
else
    disp('Not a complete data set.') % Add choice to still run data on nopt complete runs
end 

%Organize trials - put perception conditions as first two, then working
%memory as second two in a matrix (?)

%6 different combinations

% Trial 1 vs. Trial 2
[results.h12,results.p12] = ttest(run1Stats.data(:,5),run2Stats.data(:,5));
    if results.h12 == 1
        fprintf('\nThe estimated contrast differences were statistically significant between %s (trial 1) and %s (trial 2) with a calculated probability of %0.3f\n',run1Stats.condition, run2Stats.condition, results.p12);
    else
        fprintf('The estimated contrast between %s (trial 1) and %s (trial 2) is not statistically significant\n', run1Stats.condition, run2Stats.condition)
    end
 if sum(strcmp(run1Stats.condition,'Perception') || strcmp(run2Stats.condition,'Perception')) == 0
     results.trials12 = 0; % =0 is both working memory conditioned trials
 elseif sum(strcmp(run1Stats.condition,'Perception') || strcmp(run2Stats.condition,'Perception')) == 1
     results.trials12 = 1; % =1 is one working memory trial, and one perception trial
 elseif sum(strcmp(run1Stats.condition,'Perception') || strcmp(run2Stats.condition,'Perception')) == 2 
     results.trials12 = 2; % =2 is both perception conditioned trials
 else
     error('Error with condition pairings.');
 end

 % Trial 2 vs. Trial 3
[results.h23,results.p23] = ttest(run2Stats.data(:,5),run3Stats.data(:,5));
    if results.h23 == 1
        fprintf('\nThe estimated contrast differences were statistically significant between %s (trial 2) and %s (trial 3) with a calculated probability of %0.3f\n',run2Stats.condition, run3Stats.condition, results.p23);
    else
        fprintf('The estimated contrast between %s (trial 2) and %s (trial 3) is not statistically significant\n', run2Stats.condition, run3Stats.condition)
    end
 if sum(strcmp(run2Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 0
     results.trials23 = 0; % =0 is both working memory conditioned trials
 elseif sum(strcmp(run2Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 1
     results.trials23 = 1; % =1 is one working memory trial, and one perception trial
 elseif sum(strcmp(run2Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 2 
     results.trials23 = 2; % =2 is both perception conditioned trials
 else
     error('Error with condition pairings.');
 end

 




% trial 1 and 2
if sum(strcmp(run1Stats.condition,'Perception') || strcmp(run2Stats.condition,'Perception')) == 1  %One of them is perception, one is working memory
    [h,p] = ttest(run1Stats.data(:,5),run2Stats.data(:,5));
    if h == 1
        fprintf('The estimated contrast differences were statistically significant between %s (trial 1) and %s (trial 2) with a calculated probability of %0.3f',run1Stats.condition, run2Stats.condition, p);
    else
        disp('The estimated contrast between %s (trial 1) and %s (trial 2) is not statistically significant', run1Stats.condition, run2Stats.condition)
    end
elseif  sum(strcmp(run1Stats.condition,'Perception') || strcmp(run2Stats.condition,'Perception')) == 2 % Both are perception
   [h,p] = ttest(run1Stats.data(:,5),run2Stats.data(:,5));
    if h == 1
        fprintf('The estimated contrast differences were statistically significant between %s (trial 1) and %s (trial 2) with a calculated probability of %0.3f',run1Stats.condition, run2Stats.condition, p);
    else
        disp('The estimated contrast between %s (trial 1) and %s (trial 2) is not statistically significant', run1Stats.condition, run2Stats.condition)
    end
elseif sum(strcmp(run1Stats.condition,'Perception') || strcmp(run2Stats.condition,'Perception')) == 0  % Both are working memory 
    [h,p] = ttest(run1Stats.data(:,5),run2Stats.data(:,5));
    if h == 1
        fprintf('The estimated contrast differences were statistically significant between %s (trial 1) and %s (trial 2) with a calculated probability of %0.3f',run1Stats.condition, run2Stats.condition, p);
    else
        disp('The estimated contrast between %s (trial 1) and %s (trial 2) is not statistically significant', run1Stats.condition, run2Stats.condition)
    end
else 
    error('Error with condition pairings');
end

% trial 2 and 3
if sum(strcmp(run2Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 0 % Neither are perception, do nothing
elseif sum(strcmp(run2Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 1  %One of them is perception, one is working memory
    % ttest and other stats to compare
elseif  sum(strcmp(run2Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 2 % Both are perception
    % run tests between trials
end

% trial 3 and 4
if sum(strcmp(run3Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 0 % Neither are perception, do nothing
elseif sum(strcmp(run3Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 1  %One of them is perception, one is working memory
    % ttest and other stats to compare
elseif  sum(strcmp(run3Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 2 % Both are perception
    % run tests between trials
end

% trial 4 and 1
if sum(strcmp(run4Stats.condition,'Perception') || strcmp(run1Stats.condition,'Perception')) == 0 % Neither are perception, do nothing
elseif sum(strcmp(run4Stats.condition,'Perception') || strcmp(run1Stats.condition,'Perception')) == 1  %One of them is perception, one is working memory
    % ttest and other stats to compare
elseif  sum(strcmp(run4Stats.condition,'Perception') || strcmp(run1Stats.condition,'Perception')) == 2 % Both are perception
    % run tests between trials
end


% trial 2 and 4
if sum(strcmp(run2Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 0 % Neither are perception, do nothing
elseif sum(strcmp(run2Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 1  %One of them is perception, one is working memory
    % ttest and other stats to compare
elseif  sum(strcmp(run2Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 2 % Both are perception
    % run tests between trials
end 

% trial 1 and 3
if sum(strcmp(run1Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 0 % Neither are perception, do nothing
elseif sum(strcmp(run1Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 1  %One of them is perception, one is working memory
    % ttest and other stats to compare
elseif  sum(strcmp(run1Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 2 % Both are perception
    % run tests between trials
end
   
        
% run t tests between conditions of matching and not matching conditions,
% testing differencesin contrast estimations, locations, timing parameters
% Compares WM trials against themselves
% Compares Perception trials against themselves
% Compares WM and Perception trials against each other

%Do some sort of permutation thing where you get all combinations of the 4

[h,p] = ttest(run1Stats.data(:,5),run2Stats.data(:,5));
if h == 1
    fprintf(' The estimated contrast differences were statistically significant between run1 and run2 with a p value of %0.3f',p)
else
    dis('The null hypothesis could not be rejected, and the estimated contrast data is not statistically significant')
end

[h,p] = ttest(run2Stats.data(:,5),run3Stats.data(:,5));
if h == 1
    fprintf(' The estimated contrast differences were statistically significant between run1 and run2 with a p value of %0.3f',p)
else
    disp('The null hypothesis could not be rejected, and the estimated contrast data is not statistically significant')
end



    


