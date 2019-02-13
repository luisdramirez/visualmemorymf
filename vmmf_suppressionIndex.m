function suppressionIndex = suppression(baselineMat,percMat,wmMat)



for subject = subjects
        subjCount = subjCount + 1; % Subject Counter
        if e == 1
            % Perception Experiment
            baselineMat(subjCount,:) = overallData.baselineForP(subjCount,:);
            variableMat(subjCount,:) = overallData.perceptionmat(subjCount,:);
        elseif e == 2
            % Working Memory Experiment
            baselineMat(subjCount,:) = overallData.baselineForWM(subjCount,:);
            variableMat(subjCount,:) = overallData.workingmemmat(subjCount,:);
        end
        TotalSuppressionIndexVariable(e,1:size(visualmemory_subjectsRan,2),:) = ((variableMat) - (baselineMat))./((variableMat) + (baselineMat)); 
    