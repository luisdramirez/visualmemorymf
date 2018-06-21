% run_visualmemorymf.m
% Written by JS and LR June 2018
close all; clear all; clc;
commandwindow;
Screen('Preference', 'SkipSyncTests', 1);
commandwindow;
load('visualmemory_condition_order')
load('visualmemory_subjectsRan')
% visualmemory_condition_order = perms([1 2 1 2]);
% visualmemory_subjectsRan = {};
%% PREPARE
p.repetitions = 15; % data will be saved if > 5
p.numBlocks = p.repetitions;

% Subject Name
p.experiment = 'test_HC'; % 'test' = 5 contrasts ; 'test_HC' = 1 contrast, w/ baseline condition
p.subject = 'LR';

% Set directories
expDir = pwd; % set the experimental directory to the current directory 'pwd'
dataDir = 'data'; %Set the path to a directory called 'data'
t.mySeed = sum(100*clock);
rng(t.mySeed); % start with a random seed
t.theDate = datestr(now, 'yymmdd'); %collect todays date
t.timeStamp = datestr(now,'HHMM'); %collect timestamp

p.numConditions = 2;

cd(dataDir);
if exist(['data_visualmemorymf_' p.experiment '_' p.subject '.mat'],'file') ~= 0
    load(['data_visualmemorymf_' p.experiment '_' p.subject '.mat']);
    p.runNumber = length(theData)+1;
    if strcmp(p.experiment,'test') || strcmp(p.experiment,'test_HC')
        p.testCondition_curr = 1; % fixed to perception condition for hard coded testing
    else
        p.testCondition_curr = theData{1}.p.trialSchedule(p.orderRow,p.runNumber);
    end
else 
    p.runNumber = 1;
    if ~strcmp(p.experiment,'test') && ~strcmp(p.experiment,'test_HC')
        p.orderRow = length(subjectsRan)+1;
        if p.orderRow > length(visualmemory_condition_order)
           p.orderRow = p.orderRow - length(visualmemory_condition_order); 
        end
        subjectsRan{end+1} = p.subject;
        p.trialSchedule = visualmemory_condition_order(p.orderRow,:);
        % Which Test condition, run these test conditions on different days
        p.testCondition_curr = p.trialSchedule(1);
    else
        p.testCondition_curr = 1; % fixed to perception condition
    end
end
cd(expDir);

%% KEYBOARD

deviceNumber = 0;
[keyBoardIndices, ProductNames] = GetKeyboardIndices;
%deviceString = 'Apple Internal Keyboard / Trackpad';
%deviceString = 'Wired USB Keyboard';
%deviceString = 'Apple Keyboard';
%deviceString = 'USB Keyboard';
%deviceString = 'Wired Keyboard 400';
deviceString = 'Lenovo Traditional USB Keyboard';

for nTrial = 1:length(ProductNames)
    if strcmp(ProductNames{nTrial}, deviceString)
        deviceNumber = keyBoardIndices(nTrial);
        break;
    end
end
if deviceNumber == 0
    error('No device by that name was detected');
end

%Check which devicenumber the powermate is assigned to
powermate = PsychPowerMate('Open');
if isempty(powermate)
    error('problem with the powermate');
end
% % Controls the brightness of the powermate color
PsychPowerMate('SetBrightness', powermate, 20);
% 
while 1
    [pmbutton, ~] = PsychPowerMate('Get', powermate);
    if pmbutton == 1
        break;
    end
end

%% SCREEN PARAMTERS
screens=Screen('Screens');
useScreen=max(screens);
p.screenWidthPixels = Screen('Rect', useScreen);
screenWidth = 53; %cm
viewingDistance = 110; %cm
visAngle = (2*atan2(screenWidth/2, viewingDistance))*(180/pi);

p.pixPerDeg = round(p.screenWidthPixels(3)/visAngle);
p.fixation=round(0.085 * p.pixPerDeg);
colors.grey = 128; colors.white = 255; colors.green = [0 255 0]; colors.blue = [0 0 255]; colors.black = [0 0 0]; colors.red = [220 20 60]; 
colors.dimgrey = [105 105 105]; colors.yellow = [255 255 0]; colors.magenta = [255 0 255]; colors.cyan = [0 255 255];

%% GRATING PARAMETERS

% grating contrast for center and surround
p.minContrast = 0.1;
p.maxContrast = 0.75;
if strcmp(p.experiment,'test_HC')
    p.numContrasts = 1;
else
    p.numContrasts = 5;
end
p.centerContrast = [10.^linspace(log10(p.minContrast),log10(p.maxContrast),p.numContrasts)];
% p.centerContrast = rand(1)*(0.55)+0.35;

p.surroundContrast = 1;

% Grating Size 
p.centerSize = round(2*p.pixPerDeg);
p.surroundSize = p.screenWidthPixels(:,3);
p.gapSize = round(0.01*p.pixPerDeg); %space between center and surround annulus

p.ecc = 10;
p.backgroundRadius = p.ecc * p.pixPerDeg; %radius of the background circle
p.eccentricity = round((p.ecc/2) * p.pixPerDeg); %distance from center of screen to center of target
p.innerRadius = p.eccentricity - (p.centerSize/2+p.gapSize); %inner edge of surround annulus
p.outerRadius = p.eccentricity + (p.centerSize/2+p.gapSize); %outer edge of surround annulus

% Fixation dot size
p.innerFixation = round(0.075*p.pixPerDeg);
p.outerFixation = round(0.75*p.pixPerDeg);

% Grating frequency, orientation, phase
p.orientation = 0;
p.stimorientation = 90;
p.freq = 2.75;
p.frequency_center = p.centerSize/p.pixPerDeg * p.freq;
p.frequency_surround = p.surroundSize/p.pixPerDeg * p.freq;
p.numPhases = 1;
p.centerPhase = 360;
p.surroundPhase = p.centerPhase;

%% TRIAL EVENTS
% Create matrix with all unique trial events based on the number of
% repetitions, which will be saved as p.trialEvents. Stimulus
% configurations are the number of possible conditionn, manually inserted
% into BalanceFactors. BalanceFactors will output the possible combinations
% of given parameters, repeated as many times as specified. Grating
% orientation, cue, and which staircase to be fed will be specified here.
% Each of these are included into the p.trialEvents matrix
% Two conditions: perception and working memory,
% Random center grating location per trial, 
% One of five possible center grating contrasts per trial,
% Random probe location and contrast.

%--------------------%
%    Conditions      %
%--------------------%
% 1 - perception: center and surround, mask, blank, mask
% 2 - working memory: center, mask, surround, mask
% 3 - baseline: center, mask, blank, mask

% Baseline number of trials based on the number of center grating contrasts
if strcmp(p.experiment,'test_HC')
    p.numTrialsPerBlock = 5;
    p.numTrialsPerSet = p.numTrialsPerBlock*5;
    p.numBlocksPerSet = p.numTrialsPerSet/p.numTrialsPerBlock;

    col1 = [p.testCondition_curr*ones(p.numTrialsPerBlock,1); 3*ones(p.numTrialsPerBlock,1)];
    col1 = repmat(col1,p.repetitions,1);
%     col1 = repmat([1;3],length(configs)/2,1);
%     col1a = ones(length(configs)/2,1); col1b = 3*ones(length(configs)/2,1);
%     col1 = [col1a;col1b];
    p.numTrials = length(col1);
    p.numSets = p.numTrials/p.numTrialsPerSet;

else
    p.stimConfigurations = 1:length(p.centerContrast);
    [configs] = BalanceFactors(p.numBlocks,0,p.stimConfigurations);
    col1 = repmat(p.testCondition_curr,length(configs),1);
    p.numTrialsPerBlock = p.numContrasts;
    p.numTrialsPerSet = p.numTrialsPerBlock*6;
    p.numBlocksPerSet = p.numTrialsPerSet/p.numTrialsPerBlock;
    p.numTrials = length(col1);
    p.numSets = p.numTrials/p.numTrialsPerSet;
end

%--------------------%
%     Location       %
%--------------------%
%location is random every trial
col2 = randi(360,length(col1),1); % center grating location
col4 =  randi(360,length(col1),1); %probe grating location

%--------------------%
%      Contrast      %
%--------------------%
if strcmp(p.experiment,'test_HC')
    col3 = repmat(p.centerContrast,length(col1),1);
else
    col3 = repmat(p.centerContrast',p.numBlocks,1); % center grating contrast
end
col5 = rand(length(col3),1); % probe grating contrast
col5(col5>p.maxContrast)=p.maxContrast; col5(col5<p.minContrast)=p.minContrast;

% bring all 3 together
p.trialEvents = [col1 col2 col3 col4 col5]; %[condition targetLocation targetContrast probeLocation probeContrast]

if strcmp(p.experiment,'test') || strcmp(p.experiment,'test_HC')
    test = 1;
else
    test = 0;
end
shuffled = 0;
if ~strcmp(p.experiment,'test_HC')
    p.trialEvents = Shuffle(p.trialEvents,2); shuffled = 1;
end
p.trialEvents; % [condition targetLocation targetContrast probeLocation probeContrast configuration]

%% TIMING PARAMETERS
% timing is in seconds
t.stimOn1 = 2; % stimulus 1 duration    
t.retention = 2.2; % memory retention period
t.stimOn2 = 2; % stimulus 2 duration
t.iti = 2;
t.startTime = 2;
t.responseTime = []; t.responseTime_est = 5;
t.flickerTime = 0.4; 
t.flicker = 0.04;

t.trialDur = sum(t.stimOn1 + t.flickerTime + t.stimOn1 + t.flickerTime + t.responseTime_est + t.iti); % (s)
t.runDur = t.trialDur*size(p.trialEvents,1) + t.startTime*p.numBlocks; % (s)
if t.runDur/60 >= 1
    disp(['Total run duration: ' num2str(t.runDur/60) ' min(s).'])
else
    disp(['Total run duration: ' num2str(t.runDur) ' s.'])
end
%% CREATE STIMULI

%%Center
[xc,yc] = meshgrid((-p.centerSize/2):(p.centerSize/2)-1, (-p.centerSize/2):(p.centerSize/2)-1);
eccen = sqrt((xc).^2+(yc).^2); 	% calculate eccentricity of each point in grid relative to center of 2D image
centerGaussian = zeros(p.centerSize); centerGaussian(eccen <= (p.centerSize/2)) = 1;
centerTransparencyMask = zeros(p.centerSize); centerTransparencyMask(eccen <= (p.centerSize/2))=255;

%%Surround
[xs,ys] = meshgrid((-p.surroundSize/2):(p.surroundSize/2)-1, (-p.surroundSize/2):(p.surroundSize/2)-1);
eccen = sqrt((xs).^2+(ys).^2); 	% calculate eccentricity of each point in grid relative to center of 2D image
Annulus = zeros(p.surroundSize); Annulus(eccen <= p.innerRadius) = 1;
Annulus(eccen >= p.outerRadius) = 1;
surroundTransparencyMask = zeros(p.surroundSize); surroundTransparencyMask(eccen <= p.surroundSize/2)=255;

%%Background
bgAnnulus = zeros(p.surroundSize); 
bgAnnulus(eccen <= p.backgroundRadius) = 1;
bgAnnulus(eccen >= p.backgroundRadius) = 0;

% Make unique grating
[Xc,Yc] = meshgrid(0:(p.centerSize-1),0:(p.centerSize-1));
[Xs,Ys] = meshgrid(0:(p.surroundSize-1),0:(p.surroundSize-1));

centerGrating = NaN(p.centerSize,p.centerSize);
surroundGrating = NaN(p.surroundSize,p.surroundSize);

centerPatch = (sin(p.frequency_center*2*pi/p.centerSize*(Xc.*sin(p.orientation*(pi/180))+Yc.*cos(p.orientation*(pi/180)))-p.centerPhase));
centerGrating = (centerPatch .* centerGaussian);

surroundGrating = (sin(p.frequency_surround*2*pi/p.surroundSize*(Xs.*sin(p.orientation*(pi/180))+Ys.*cos(p.orientation*(pi/180)))-p.surroundPhase));
tmpSurround = (surroundGrating .* Annulus);
tmpSurround(bgAnnulus == 0) = -1;
surroundGrating = tmpSurround;
%% MASK

%sf filter
cutoff = [1 3] .*round(p.surroundSize/p.pixPerDeg);
f = freqspace(p.surroundSize);

%bandpass filter
filter = Bandpass2(p.surroundSize, f(cutoff(1)), f(cutoff(2)));

%low pass filter
h = fspecial('disk', 5);
filter = conv2(filter, h, 'same');
filter = filter/max(filter(:));

%noise
maskGrating = NaN(round(t.flickerTime/t.flicker), p.surroundSize, p.surroundSize);
for nTrial = 1:(t.flickerTime/t.flicker)
    noise = -1+2.*rand(p.surroundSize);
    fftNoise = fftshift(fft2(noise));
    filterNoise = fftNoise .* filter;
    newNoise = real(ifft2(fftshift(filterNoise)));
    noiseMask = newNoise./(max(abs(newNoise(:))));
    tempMask = (noiseMask.*bgAnnulus);
    tempMask(bgAnnulus == 0) = -1; %to make black background make == 0
    maskGrating (nTrial,:,:) = tempMask;
end 

%% WINDOW SETUP

if shuffled == 1 || test == 1
    [window,rect] = Screen('OpenWindow', useScreen, colors.black, []); %test screen size [0 0 700 500]
    
    OriginalCLUT = Screen('ReadNormalizedGammaTable', window);
    
    % Enable alpha blending
    Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    % Define coordinates where to draw the fixation
    CenterX = rect(3)/2; CenterY = rect(4)/2;
    % Coordinates for location on left and right side of fixation
    center =  [CenterX CenterY];
    
    Screen('TextStyle', window, 1);
    Screen('TextSize', window, 16);
    t.ifi = Screen('GetFlipInterval',window); % grab screen refresh rate
    
elseif shuffled == 0 && test == 0
    error('Trials not shuffled!')
end
%% EXPERIMENT LOOP
HideCursor;
% Esc to quit
StartKey=zeros(1,256); StartKey(KbName({'ESCAPE'})) = 1;
PsychHID('KbQueueCreate', deviceNumber, StartKey);
PsychHID('KbQueueStart', deviceNumber);

% Generate mask texture
Mask = NaN(round(t.flickerTime/t.flicker),1);
for m = 1:round(t.flickerTime/t.flicker)
    Mask(m,:) = Screen('MakeTexture', window, squeeze(maskGrating(m,:,:))* colors.grey + colors.grey);
end

% Welcome Screen
Screen('TextStyle', window, 1);
Screen('TextSize', window, 16);
welcomeText = ['Hello' '\n' ...
    'Click powermate to start experiment.'];
DrawFormattedText(window, welcomeText, 'center', 'center', 255);
Screen('Flip', window);
while 1
    [pmbutton, ~] = PsychPowerMate('Get', powermate);
    if pmbutton == 1
        break;
    end
end

% Starting Screen
Screen('FillOval', window, colors.grey,  [CenterX-p.backgroundRadius CenterY-p.backgroundRadius CenterX+p.backgroundRadius CenterY+p.backgroundRadius]);
Screen('FillOval', window, colors.black,  [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation]);
Screen('FillOval' , window, colors.green, [CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
Screen('Flip',window);
WaitSecs(t.startTime);

nSet = 1;
for nTrial = 1:size(p.trialEvents,1)
    
    %--------------------%
    %    Trial Settings  %
    %--------------------%  
    p.currentDegreeLocation = round(p.trialEvents(nTrial,2));
    p.centerContrast = p.trialEvents(nTrial,3);
    targetXY = round([CenterX+p.eccentricity*cos(p.currentDegreeLocation*(pi/180))' CenterY-p.eccentricity*sin(p.currentDegreeLocation*(pi/180))']);

    %--------------------%
    %    Stimulus 1      %
    %--------------------%
    %If perception condition, draw surround annulus with center grating
    if p.trialEvents(nTrial,1) == 1 
    surroundTexture = (surroundGrating*p.surroundContrast)*colors.grey + colors.grey;
    surroundTexture(:,:,2) = surroundTransparencyMask;
    surroundStimulus = Screen('MakeTexture', window, surroundTexture);
    Screen('DrawTexture', window, surroundStimulus, [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], CenterX, CenterY), p.stimorientation);
    else
        Screen('FillOval', window, colors.grey,  [CenterX-p.backgroundRadius CenterY-p.backgroundRadius CenterX+p.backgroundRadius CenterY+p.backgroundRadius]);
    end   
    % Draw center grating
    centerTexture = (centerGrating*p.centerContrast)*colors.grey + colors.grey;
    centerTexture(:,:,2) = centerTransparencyMask;
    centerStimulus = Screen('MakeTexture', window, centerTexture);
    Screen('DrawTexture', window, centerStimulus, [], ...
            CenterRectOnPoint([0 0 p.centerSize p.centerSize], targetXY(1), targetXY(2)), p.stimorientation);  
    % Draw fixation
    Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation])
    Screen('FillOval', window,colors.green,[CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
    Screen('Flip',window);
    WaitSecs(t.stimOn1);
    
    %--------------------%
    %       Mask 1       %
    %--------------------%  
    indx = Shuffle(1:round(t.flickerTime/t.flicker));
    StartMask = GetSecs;
    for a = indx
        if GetSecs > StartMask + t.flickerTime - t.flicker
            break
        end        
        Screen('DrawTexture', window, Mask(a,:), [],CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], CenterX, CenterY));
        Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation]);
        Screen('FillOval', window, colors.green, [CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
        Screen('Flip', window);
        WaitSecs(t.flicker);        
    end
    
    %--------------------%
    %    Stimulus 2      %
    %--------------------%
    %%% Retention interval 1
    Screen('FillOval', window, colors.grey,  [CenterX-p.backgroundRadius CenterY-p.backgroundRadius CenterX+p.backgroundRadius CenterY+p.backgroundRadius]);
    % Draw fixation
    Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation])
    Screen('FillOval', window,colors.green,[CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
    Screen('Flip',window);
    WaitSecs(t.retention);
    
    %%% Draw surround if memory condition
    Screen('FillOval', window, colors.grey,  [CenterX-p.backgroundRadius CenterY-p.backgroundRadius CenterX+p.backgroundRadius CenterY+p.backgroundRadius]);
    % If memory condition, display surround annulus alone
    if p.trialEvents(nTrial,1) == 2
        surroundTexture = (surroundGrating*p.surroundContrast)*colors.grey + colors.grey;
        surroundTexture(:,:,2) = surroundTransparencyMask;
        surroundStimulus = Screen('MakeTexture', window, surroundTexture);
        Screen('DrawTexture', window, surroundStimulus, [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], CenterX, CenterY), p.stimorientation);
    end
    % Draw fixation
    Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation])
    Screen('FillOval', window,colors.green,[CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
    Screen('Flip',window);
    WaitSecs(t.stimOn2);

    %%% Retention interval 2
    Screen('FillOval', window, colors.grey,  [CenterX-p.backgroundRadius CenterY-p.backgroundRadius CenterX+p.backgroundRadius CenterY+p.backgroundRadius]);
    % Draw fixation
    Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation])
    Screen('FillOval', window,colors.green,[CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
    Screen('Flip',window);
    WaitSecs(t.retention);
    
    %--------------------%
    %       Mask 2       %
    %--------------------%  
    indx = Shuffle(1:round(t.flickerTime/t.flicker));
    StartMask = GetSecs;
    for a = indx
        if GetSecs > StartMask + t.flickerTime - t.flicker
            break
        end        
        Screen('DrawTexture', window, Mask(a,:), [],CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], CenterX, CenterY));
        Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation]);
        Screen('FillOval', window, colors.green, [CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
        Screen('Flip', window);
        WaitSecs(t.flicker);        
    end

    %--------------------%
    %       Probe        %
    %--------------------%
    
    % Check if ESCAPE has been pressed
    [keyIsDown, keyCode] = PsychHID('KbQueueCheck', deviceNumber); %check response
    key = find(keyCode);
    if key == KbName('ESCAPE') % windows = 'esc', mac = 'ESCAPE' If user presses ESCAPE, exit the program.
        Screen('LoadNormalizedGammaTable', window, OriginalCLUT);
        Screen('CloseAll');
        ListenChar(1); % % Go back to unsuppressed mode
        FlushEvents('keyDown', deviceNumber);
        error('User exited program.');
    end

    %Show probe at random location with random contrast
    Screen('FillOval', window, colors.grey, [CenterX-p.backgroundRadius CenterY-p.backgroundRadius CenterX+p.backgroundRadius CenterY+p.backgroundRadius]);
    
    initialAngle = p.trialEvents(nTrial,4);
    ProbeXY = round([CenterX + p.eccentricity*(cos(initialAngle*(pi/180)))' CenterY - p.eccentricity*(sin(initialAngle*(pi/180)))']);
    intial_contrast = p.trialEvents(nTrial,5); % random start contrast on each trial
    
    
    Probe = Screen('MakeTexture', window, squeeze(centerGrating)* (intial_contrast *colors.grey) + colors.grey);
    Screen('DrawTexture', window, Probe, [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], ProbeXY(1), ProbeXY(2)), p.stimorientation);
    Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation]);
    Screen('FillOval', window, colors.green, [CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
    Screen('Flip', window);
    starttrial = GetSecs; % get the start time of each trial
    
      
%     Allow for dial rotation for location update
    [~, startangle] = PsychPowerMate('Get',powermate);
 
    while 1 %start inf loop
        % Query PowerMate button state and rotation angle in "clicks"
        [pmbutton, angle] = PsychPowerMate('Get', powermate);
        % 1st button is the "or" of the 1st mouse button and the actual PowerMate button
        if startangle ~= angle
             
            % Convert turn of dial first to degrees and then to contrast:
            Screen('FillOval', window, colors.grey, [CenterX-p.backgroundRadius CenterY-p.backgroundRadius CenterX+p.backgroundRadius CenterY+p.backgroundRadius]);
            
            angles = ((startangle-angle)*3.8298);
            changeposition = (angles/(2*pi));
            initialAngle = initialAngle + changeposition; % update the location relative to last dial position
            ProbeXY = round([CenterX + p.eccentricity*(cos(initialAngle*(pi/180)))' CenterY - p.eccentricity*(sin(initialAngle*(pi/180)))']);
            
            Target = Screen('MakeTexture', window, squeeze(centerGrating)* (intial_contrast*colors.grey) + colors.grey);
            Screen('DrawTexture', window, Target, [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], ProbeXY(1), ProbeXY(2)), p.stimorientation)
            Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation]);
            Screen('FillOval', window, colors.green, [CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation])
            Screen('Flip', window);
            
            startangle = angle;
        end
        if pmbutton == 1;
            locationTime = GetSecs;
% %             % make sure angle stays in 0-360 range
            correctedAngle = mod(initialAngle, 360);
            
            data.EstimatedLocation(nTrial) = correctedAngle;
%             
% %             %make sure difference is in the 180 range
            difference = abs(p.trialEvents(nTrial,2) - data.EstimatedLocation(nTrial));
            
            if difference > 180
                difference = abs(difference - 360);
            end
            
            data.DifferenceLocation(nTrial) = difference;
            
            data.ResponseTime_location(nTrial) = (locationTime - starttrial);
            pmbutton = 0;
            break
        end
    end
%     
%     % PowerMate is sampled at 10msec intervals, therefore have a short
%     % break to make sure it doesn't skip the contrast task
    WaitSecs(1);
%     
%     %Button press for contrast
    [~, contrastangle] = PsychPowerMate('Get', powermate);
    while 1 %start inf loop
%         % Query PowerMate button state and rotation angle in "clicks"
        [pmbutton_contrast, angle2] = PsychPowerMate('Get', powermate);
%         % 1st button is the "or" of the 1st mouse button and the actual PowerMate button
        if contrastangle ~= angle2
%             
%             % Convert turn of dial first to degrees and then to contrast:
                        angles = (angle * 3.8298)/360;
            angles = ((contrastangle-angle2)*3.8298);
            changecontrast = angles/360;
            intial_contrast = intial_contrast - changecontrast; % update the contrast relative to last dial position
%             % Make sure we stay in range
%             
            if intial_contrast > 1
                intial_contrast = 1;
            elseif intial_contrast < 0
                intial_contrast = 0.001;
            end
            Screen('FillOval', window, colors.grey, [CenterX-p.backgroundRadius CenterY-p.backgroundRadius CenterX+p.backgroundRadius CenterY+p.backgroundRadius]);
            
            Target = Screen('MakeTexture', window, squeeze(centerGrating)* (intial_contrast*colors.grey) + colors.grey);
            Screen('DrawTexture', window, Target, [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], ProbeXY(1), ProbeXY(2)), p.stimorientation)
            Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation]);
            Screen('FillOval', window, colors.green, [CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation])
            Screen('Flip', window);
            
            contrastangle = angle2;
        end
        if pmbutton_contrast == 1;
            data.EstimatedContrast(nTrial) = intial_contrast;
            data.DifferenceContrast(nTrial) = p.trialEvents(nTrial,3) - data.EstimatedContrast(nTrial);
            data.ResponseTime(nTrial) = (GetSecs - starttrial);
            pmbutton_contrast = 0;
            break
        end     
    end

    %--------------------%
    %       Break        %
    %--------------------%  
    if mod(nTrial,p.numTrialsPerSet) == 0 && nTrial == p.numTrialsPerSet*nSet 
        rest = GetSecs;        
        Screen('FillOval', window, colors.grey, [CenterX-p.backgroundRadius CenterY-p.backgroundRadius CenterX+p.backgroundRadius CenterY+p.backgroundRadius]);
        Screen('TextStyle', window, 1);
        Screen('TextSize', window, 16);
        breakText = ['You make take a short break now.\n Or press the powermate to continue.\n'];
        DrawFormattedText(window, breakText, 'center', 'center', colors.white);
        Screen('Flip', window);
        restText = ['You can take a short break now, ' '' '\n' ...
            'or press the dial to continue' '\n' '\n' ...
            num2str(nSet-1) '/' num2str(p.numSets) ' completed.' ];
        DrawFormattedText(window, restText, 'center', 'center', colors.white);
        Screen('Flip', window);
        pmbuttonbreak = 0;
        WaitSecs(1 );
        while 1
            [pmbuttonbreak, a] = PsychPowerMate('Get', powermate);
            if pmbuttonbreak == 1
                break;
            end
        end
        t.restTime(nSet) = (GetSecs-rest)/60;
        nSet = nSet + 1;
    end

    % ITI
    Screen('FillOval', window, colors.grey,  [CenterX-p.backgroundRadius CenterY-p.backgroundRadius CenterX+p.backgroundRadius CenterY+p.backgroundRadius]);
    Screen('FillOval', window, colors.black,  [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation]);
    Screen('FillOval' , window, colors.green, [CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
    Screen('Flip',window);
    WaitSecs(t.iti);
end

%Draw some more text to the screen outside of the loop:
Screen(window,'TextSize',30);
ByebyeText = ['Great work! You have completed this run.' '\n' '\n' ...
    'Please let the experimenter know you have finished.'];
Screen('FillOval', window, colors.grey, [CenterX-p.backgroundRadius CenterY-p.backgroundRadius CenterX+p.backgroundRadius CenterY+p.backgroundRadius]);
DrawFormattedText(window, ByebyeText, 'center', 'center', colors.white);
Screen('Flip', window);
WaitSecs(3);

Screen('LoadNormalizedGammaTable', window, OriginalCLUT);
Screen('CloseAll')
ShowCursor;
%% SAVE OUT THE DATA FILE
 if p.repetitions > 5
    cd(dataDir);
    theData(p.runNumber).t = t;
    theData(p.runNumber).p = p;
    theData(p.runNumber).data = data;
    save(['data_visualmemorymf_' p.experiment '_' p.subject '.mat'], 'theData')
    cd(expDir);
end