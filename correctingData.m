clear all; close all; clc;
p.experiment = 'exp';
p.subject = '006';
expDir = pwd;
dataDir = 'data_master';
cd(dataDir)
load(['data_visualmemorymf_' p.experiment '_' p.subject '.mat'])

for n = 1:numel(theData)
theData_corrected(n).t=theData(n).t;
end

for n = 1:numel(theData)
theData_corrected(n).p=theData(n).p;
end

for n = 1:numel(theData)
theData_corrected(n).data.EstimatedLocation = theData(n).data.EstimatedLocation;
theData_corrected(n).data.DifferenceLocation = theData(n).data.DifferenceLocation;
theData_corrected(n).data.ResponseTime_location = theData(n).data.ResponseTime_location;
theData_corrected(n).data.EstimatedContrast = theData(n).data.EstimatedContrast;
theData_corrected(n).data.DifferenceContrast = theData(n).data.DifferenceContrast;
theData_corrected(n).data.ResponseTime_Contrast = theData(n).data.ResponseTime_Contrast;
theData_corrected(n).data.responseTimes = theData(n).data.responseTime;
end

theData = theData_corrected;
save(['data_visualmemorymf_' p.experiment '_' p.subject '.mat'], 'theData')
cd(expDir)