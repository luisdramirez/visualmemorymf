%% Analysis Designed to Compare statistical results between data sets of different runs of the same subject %%
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
    disp('Not a complete data set.') % Add choice to still run data on not complete runs
end 

possCombs = nchoosek(4,2);
% 6 different combinations of comparison for trials 1-4
% 3 different potential matchups - p & p, wm & wm, wm & p

% Trial 1 vs. Trial 2
[results.trials12.h,results.trials12.p] = ttest(run1Stats.data(:,5),run2Stats.data(:,5));
    if results.h12 == 1
        fprintf('\nThe estimated contrast differences were statistically significant between %s (trial 1) and %s (trial 2) with a calculated probability of %0.3f\n',run1Stats.condition, run2Stats.condition, results.trials12.p);
    else
        fprintf('The estimated contrast between %s (trial 1) and %s (trial 2) is not statistically significant\n', run1Stats.condition, run2Stats.condition)
    end
 if sum(strcmp(run1Stats.condition,'Perception') || strcmp(run2Stats.condition,'Perception')) == 0
     results.trials12.condition = 0; % =0 is both working memory conditioned trials
 elseif sum(strcmp(run1Stats.condition,'Perception') || strcmp(run2Stats.condition,'Perception')) == 1
     results.trials12.condition = 1; % =1 is one working memory trial, and one perception trial
 elseif sum(strcmp(run1Stats.condition,'Perception') || strcmp(run2Stats.condition,'Perception')) == 2 
     results.trials12.condition = 2; % =2 is both perception conditioned trials
 else
     error('Error with condition pairings.');
 end

 % Trial 2 vs. Trial 3
[results.trials23.h,results.trials23.p] = ttest(run2Stats.data(:,5),run3Stats.data(:,5));
    if results.h23 == 1
        fprintf('\nThe estimated contrast differences were statistically significant between %s (trial 2) and %s (trial 3) with a calculated probability of %0.3f\n',run2Stats.condition, run3Stats.condition, results.trials23.p);
    else
        fprintf('The estimated contrast between %s (trial 2) and %s (trial 3) is not statistically significant\n', run2Stats.condition, run3Stats.condition)
    end
 if sum(strcmp(run2Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 0
     results.trials23.condition = 0; % =0 is both working memory conditioned trials
 elseif sum(strcmp(run2Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 1
     results.trials23.condition = 1; % =1 is one working memory trial, and one perception trial
 elseif sum(strcmp(run2Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 2 
     results.trials23.condition = 2; % =2 is both perception conditioned trials
 else
     error('Error with condition pairings.');
 end
 
 % Trial 3 vs. Trial 4
[results.trials34.h,results.trials34.p] = ttest(run3Stats.data(:,5),run4Stats.data(:,5));
if results.trials34.h == 1
    fprintf('\nThe estimated contrast differences were statistically significant between %s (trial 3) and %s (trial 4) with a calculated probability of %0.3f\n',run3Stats.condition, run4Stats.condition, results.trials34.p);
else
    fprintf('The estimated contrast between %s (trial 3) and %s (trial 4) is not statistically significant\n', run3Stats.condition, run4Stats.condition)
end
 if sum(strcmp(run3Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 0
     results.trials34.condition = 0; % =0 is both working memory conditioned trials
 elseif sum(strcmp(run3Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 1
     results.trials34.condition = 1; % =1 is one working memory trial, and one perception trial
 elseif sum(strcmp(run3Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 2 
     results.trials34.condition = 2; % =2 is both perception conditioned trials
 else
     error('Error with condition pairings.');
 end

 % Trial 4 vs. Trial 1
[results.trials41.h,results.trials41.p] = ttest(run4Stats.data(:,5),run1Stats.data(:,5));
    if results.trials41.h == 1
        fprintf('\nThe estimated contrast differences were statistically significant between %s (trial 4) and %s (trial 1) with a calculated probability of %0.3f\n',run4Stats.condition, run1Stats.condition, results.trial41.p);
    else
        fprintf('The estimated contrast between %s (trial 4) and %s (trial 1) is not statistically significant\n', run4Stats.condition, run1Stats.condition)
    end
 if sum(strcmp(run4Stats.condition,'Perception') || strcmp(run1Stats.condition,'Perception')) == 0
     results.trials41.condition = 0; % =0 is both working memory conditioned trials
 elseif sum(strcmp(run4Stats.condition,'Perception') || strcmp(run1Stats.condition,'Perception')) == 1
     results.trials41.condition = 1; % =1 is one working memory trial, and one perception trial
 elseif sum(strcmp(run4Stats.condition,'Perception') || strcmp(run1Stats.condition,'Perception')) == 2 
     results.trials41.condition = 2; %=2 is both perception conditioned trials
 else
     error('Error with condition pairings.');
 end
 
 % Trial 1 vs. Trial 3
[results.trials13.h,results.trials13.p] = ttest(run1Stats.data(:,5),run3Stats.data(:,5));
    if results.trials13.h == 1
        fprintf('\nThe estimated contrast differences were statistically significant between %s (trial 1) and %s (trial 3) with a calculated probability of %0.3f\n',run1Stats.condition, run3Stats.condition, results.trials23.p);
    else
        fprintf('The estimated contrast between %s (trial 2) and %s (trial 3) is not statistically significant\n', run1Stats.condition, run3Stats.condition)
    end
 if sum(strcmp(run1Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 0
     results.trials13.condition = 0; % =0 is both working memory conditioned trials
 elseif sum(strcmp(run1Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 1
     results.trials13.condition= 1; % =1 is one working memory trial, and one perception trial
 elseif sum(strcmp(run1Stats.condition,'Perception') || strcmp(run3Stats.condition,'Perception')) == 2 
     results.trials23.condition = 2; % =2 is both perception conditioned trials
 else
     error('Error with condition pairings.');
 end

% Trial 2 vs. Trial 4
[results.trials24.h,results.trials24.p] = ttest(run2Stats.data(:,5),run4Stats.data(:,5));
    if results.trials24.h == 1
        fprintf('\nThe estimated contrast differences were statistically significant between %s (trial 2) and %s (trial 4) with a calculated probability of %0.3f\n',run2Stats.condition, run4Stats.condition, results.trials24.p);
    else
        fprintf('The estimated contrast between %s (trial 2) and %s (trial 4) is not statistically significant\n', run2Stats.condition, run4Stats.condition)
    end
 if sum(strcmp(run2Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 0
     results.trials24.condition = 0; % =0 is both working memory conditioned trials
 elseif sum(strcmp(run2Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 1
     results.trials24.condition = 1; % =1 is one working memory trial, and one perception trial
 elseif sum(strcmp(run2Stats.condition,'Perception') || strcmp(run4Stats.condition,'Perception')) == 2 
     results.trials24.condition = 2; % =2 is both perception conditioned trials
 else
     error('Error with condition pairings.');
 end

% Organize trials by their conditions 
    % first - both wm
    % 2-5 - mixes of working memoryt and perception
    % 6 - both perception
    % make a counter that matches them equal to 0 1 or 2
    
Afields = fieldnames(results);
Acell = struct2cell(A);
resultSize = size(Acell);
Acell = reshape(Acell,resultSize, []);
Acell = Acell';
Acell = sortrows(Acell,3); %sorts the cell by the condition, 0 - 2
Acell = reshape(Acell',resultSize);
Acell = cell2struct(Acell,Afields,1);
condSet = [0 1 1 1 1 2];
%converted back to struct, sorted by condition 0 - 2 starting with condition.
wm = NaN(1);
wm_p = NaN(1:4);
p = NaN(1);
for i = 1:length(condSet)
    counter = 1; 
    %% trials 1-2
    if condSet(i) == results.trials12.condition
        if ismember(1,counter) == 1
            if isnan(wm) == 1
            wm = results.trials12;
            else 
                error('Index into to a condition vector is full');
            end
        elseif ismember([2,3,4,5],counter) == 1
           if sum(isnan(wm_p)) == 4
                wm_p(1) = results.trials12;
            elseif sum(isnan(wm_p)) == 3
                wm_p(2) = results.trials12;
            elseif sum(isnan(wm_p)) == 2
                wm_p(3) = results.trials12;
           elseif sum(isnan(wm_p)) == 1
               wm_p(4) = results.trials12;
           else
               error('Index into w condition vector is full.')
           end
        elseif ismember(6,counter) == 1
           if isnan(p) == 1
            p = results.trials12;
            else 
                error('Index into to a condition vector is full');
            end
        end 
     %%trials 2-3       
    elseif condSet(i) == results.trials23.condition
         if ismember(1,counter) == 1
            if isnan(wm) == 1
            wm = results.trials23;
            else 
                error('Index into to a condition vector is full');
            end
        elseif ismember([2,3,4,5],counter) == 1
           if sum(isnan(wm_p)) == 4
                wm_p(1) = results.trials23;
            elseif sum(isnan(wm_p)) == 3
                wm_p(2) = results.trials23;
            elseif sum(isnan(wm_p)) == 2
                wm_p(3) = results.trials23;
           elseif sum(isnan(wm_p)) == 1
               wm_p(4) = results.trials23;
           else
               error('Index into w condition vector is full.')
           end
        elseif ismember(6,counter) == 1
           if isnan(p) == 1
            p = results.trials23;
            else 
                error('Index into to a condition vector is full');
            end
         end 
        %%trials 3-4
    elseif condSet(i) == results.trials34.condition
        if ismember(1,counter) == 1
            if isnan(wm) == 1
            wm = results.trials34;
            else 
                error('Index into to a condition vector is full');
            end
        elseif ismember([2,3,4,5],counter) == 1
           if sum(isnan(wm_p)) == 4
                wm_p(1) = results.trials12;
            elseif sum(isnan(wm_p)) == 3
                wm_p(2) = results.trials12;
            elseif sum(isnan(wm_p)) == 2
                wm_p(3) = results.trials12;
           elseif sum(isnan(wm_p)) == 1
               wm_p(4) = results.trials12;
           else
               error('Index into w condition vector is full.')
           end
        elseif ismember(6,counter) == 1
           if isnan(p) == 1
            p = results.trials12;
            else 
                error('Index into to a condition vector is full');
            end
        end
   %% Trials 4-1
    elseif condSet(i) == results.trials41.condition
        if ismember(1,counter) == 1
            if isnan(wm) == 1
            wm = results.trials41;
            else 
                error('Index into to a condition vector is full');
            end
        elseif ismember([2,3,4,5],counter) == 1
           if sum(isnan(wm_p)) == 4
                wm_p(1) = results.trials41;
            elseif sum(isnan(wm_p)) == 3
                wm_p(2) = results.trials41;
            elseif sum(isnan(wm_p)) == 2
                wm_p(3) = results.trials41;
           elseif sum(isnan(wm_p)) == 1
               wm_p(4) = results.trials41;
           else
               error('Index into w condition vector is full.')
           end
        elseif ismember(6,counter) == 1
           if isnan(p) == 1
            p = results.trials41;
            else 
                error('Index into to a condition vector is full');
           end
       end 
        
    %trials 1-3
    elseif condSet(i) == results.trials13.condition
        if ismember(1,counter) == 1
            if isnan(wm) == 1
            wm = results.trials13;
            else 
                error('Index into to a condition vector is full');
            end
        elseif ismember([2,3,4,5],counter) == 1
           if sum(isnan(wm_p)) == 4
                wm_p(1) = results.trials13;
            elseif sum(isnan(wm_p)) == 3
                wm_p(2) = results.trials13;
            elseif sum(isnan(wm_p)) == 2
                wm_p(3) = results.trials13;
           elseif sum(isnan(wm_p)) == 1
               wm_p(4) = results.trials13;
           else
               error('Index into w condition vector is full.')
           end
        elseif ismember(6,counter) == 1
           if isnan(p) == 1
            p = results.trials13;
            else 
                error('Index into to a condition vector is full');
            end
        end 
        
        %trials 2-4
    elseif condSet(i) == results.trials24
        if ismember(1,counter) == 1
            if isnan(wm) == 1
            wm = results.trials12;
            else 
                error('Index into to a condition vector is full');
            end
        elseif ismember([2,3,4,5],counter) == 1
           if sum(isnan(wm_p)) == 4
                wm_p(1) = results.trials12;
            elseif sum(isnan(wm_p)) == 3
                wm_p(2) = results.trials12;
            elseif sum(isnan(wm_p)) == 2
                wm_p(3) = results.trials12;
           elseif sum(isnan(wm_p)) == 1
               wm_p(4) = results.trials12;
           else
               error('Index into w condition vector is full.')
           end
        elseif ismember(6,counter) == 1
           if isnan(p) == 1
            p = results.trials12;
            else 
                error('Index into to a condition vector is full');
            end
        end 
    counter = counter + 1;
    end
end
% Now we have vector arrays organizing what trial consitions were compared
avgP = mean(wm_p(1:4).condition);
fprintf(' The average p value for comparisons of working memory and perception is %.2f',avgP);

%percentages for differences

        
% run t tests between conditions of matching and not matching conditions,
% testing differencesin contrast estimations, locations, timing parameters
% Compares WM trials against themselves
% Compares Perception trials against themselves
% Compares WM and Perception trials against each other




    


