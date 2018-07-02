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
 
 % Trial 3 vs. Trial 4
[results.h34,results.p34] = ttest(run3Stats.data(:,5),run4Stats.data(:,5));
    if results.h34 == 1
        fprintf('\nThe estimated contrast differences were statistically significant between %s (trial 3) and %s (trial 4) with a calculated probability of %0.3f\n',run3Stats.condition, run4Stats.condition, results.p34);
    else
        fprintf('The estimated contrast between %s (trial 3) and %s (trial 4) is not statistically significant\n', run3Stats.condition, run4Stats.condition)
    end
 if sum(strcmp(run3Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 0
     results.trials34 = 0; % =0 is both working memory conditioned trials
 elseif sum(strcmp(run3Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 1
     results.trials34 = 1; % =1 is one working memory trial, and one perception trial
 elseif sum(strcmp(run3Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 2 
     results.trials = 2; % =2 is both perception conditioned trials
 else
     error('Error with condition pairings.');
 end

 % Trial 4 vs. Trial 1
[results.h41,results.p41] = ttest(run4Stats.data(:,5),run1Stats.data(:,5));
    if results.h41 == 1
        fprintf('\nThe estimated contrast differences were statistically significant between %s (trial 4) and %s (trial 1) with a calculated probability of %0.3f\n',run4Stats.condition, run1Stats.condition, results.p41);
    else
        fprintf('The estimated contrast between %s (trial 4) and %s (trial 1) is not statistically significant\n', run4Stats.condition, run1Stats.condition)
    end
 if sum(strcmp(run4Stats.condition,'Perception') || strcmp(run1Stats.condition,'Perception')) == 0
     results.trials41 = 0; % =0 is both working memory conditioned trials
 elseif sum(strcmp(run4Stats.condition,'Perception') || strcmp(run1Stats.condition,'Perception')) == 1
     results.trials41 = 1; % =1 is one working memory trial, and one perception trial
 elseif sum(strcmp(run4Stats.condition,'Perception') || strcmp(run1Stats.condition,'Perception')) == 2 
     results.trials41 = 2; % =2 is both perception conditioned trials
 else
     error('Error with condition pairings.');
 end
 
 % Trial 1 vs. Trial 3
[results.h13,results.p13] = ttest(run1Stats.data(:,5),run3Stats.data(:,5));
    if results.h13 == 1
        fprintf('\nThe estimated contrast differences were statistically significant between %s (trial 1) and %s (trial 3) with a calculated probability of %0.3f\n',run1Stats.condition, run3Stats.condition, results.p23);
    else
        fprintf('The estimated contrast between %s (trial 2) and %s (trial 3) is not statistically significant\n', run1Stats.condition, run3Stats.condition)
    end
 if sum(strcmp(run1Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 0
     results.trials13 = 0; % =0 is both working memory conditioned trials
 elseif sum(strcmp(run1Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 1
     results.trials13 = 1; % =1 is one working memory trial, and one perception trial
 elseif sum(strcmp(run1Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 2 
     results.trials23 = 2; % =2 is both perception conditioned trials
 else
     error('Error with condition pairings.');
 end

% Trial 2 vs. Trial 4
[results.h24,results.p24] = ttest(run2Stats.data(:,5),run4Stats.data(:,5));
    if results.h24 == 1
        fprintf('\nThe estimated contrast differences were statistically significant between %s (trial 2) and %s (trial 4) with a calculated probability of %0.3f\n',run2Stats.condition, run4Stats.condition, results.p24);
    else
        fprintf('The estimated contrast between %s (trial 2) and %s (trial 4) is not statistically significant\n', run2Stats.condition, run4Stats.condition)
    end
 if sum(strcmp(run2Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 0
     results.trials24 = 0; % =0 is both working memory conditioned trials
 elseif sum(strcmp(run2Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 1
     results.trials24 = 1; % =1 is one working memory trial, and one perception trial
 elseif sum(strcmp(run2Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 2 
     results.trials24 = 2; % =2 is both perception conditioned trials
 else
     error('Error with condition pairings.');
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



    


