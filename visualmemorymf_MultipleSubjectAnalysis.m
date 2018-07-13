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
files = dir(dataDir);
load('visualmemory_condition_order')
load('visualmemory_subjectsRan')  
possibleFileNames
for i = length(visualmemory_subjectsRan)
    
possibleFileNames = 
for i = 1:length(files)
    if sum(strcmp(files.name(i) == 'data_visualmemorymf_exp_'
       load(['data_visualmemorymf_' p.experiment '_' p.subject '.mat']);