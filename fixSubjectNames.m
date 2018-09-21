% fix subject names

subjects = {'001' '002' '003' '004' '005'};
expDir = pwd;
dataDir = 'data_master';

cd(dataDir)
for nSubj = 1:numel(subjects)
    load(['data_visualmemorymf_exp_' subjects{nSubj} '.mat'])
    for nRun = 1:numel(theData)
        theData(nRun).p.subject = subjects{nSubj};
    end
    save(['data_visualmemorymf_exp_' subjects{nSubj} '.mat'],'theData')
end

cd(expDir)