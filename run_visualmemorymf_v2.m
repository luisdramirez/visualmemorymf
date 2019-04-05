%%% run_visualmemorymf_v2

clc; clear all; close all;
commandwindow;
Screen('Preference', 'SkipSyncTests', 1);
test_env = 0;
%% Initial Setup

p.repetitions = 15;
p.subject = '1';

fileName = ['data_vmmf_v2_' p.subject '.mat'];

if strcmp(p.subject,'test')
    test_env = 1;
end

% Set directories
expDir = pwd; % set the experimental directory to the current directory 'pwd'
dataDir = 'data_master'; %Set the path to a directory called 'data'
t.mySeed = sum(100*clock);
rng(t.mySeed); % start with a random seed
t.theDate = datestr(now, 'yymmdd'); %collect todays date
t.timeStamp = datestr(now,'HHMM'); %collect timestamp
cd(dataDir);

if exist(fileName,'file') ~= 0
    load(fileName)
    p.runNumber = length(theData)+1;
else
    p.runNumber = 1;
end
cd(expDir);

if test_env
    disp('!ENTERING TEST ENVIRONMENT!')
end
%% Input Setup

deviceNumber = 0;
[keyBoardIndices, ProductNames] = GetKeyboardIndices;

deviceString = 'Lenovo Traditional USB Keyboard'; %rm208
% deviceString = 'Apple Internal Keyboard / Trackpad';
% deviceString = 'USB-HID Keyboard'; % luis' desk keyboard
% deviceString = 'Wired USB Keyboard';
% deviceString = 'Apple Keyboard';
% deviceString = 'USB Keyboard';
% deviceString = 'Wired Keyboard 400';

for nTrial = 1:length(ProductNames)
    if strcmp(ProductNames{nTrial}, deviceString)
        deviceNumber = keyBoardIndices(nTrial);
        break;
    end
end
if deviceNumber == 0
    error('No device by that name was detected');
end

if ~test_env
    %Check which device number the powermate is assigned to
    powermate = PsychPowerMate('Open');
    if isempty(powermate)
        error('problem with the powermate or test_env has not been set to 0');
    end
    % % Controls the brightness of the powermate color
    PsychPowerMate('SetBrightness', powermate, 20);
    disp('Click powermate to continue.')
    while 1
        [pmbutton, ~] = PsychPowerMate('Get', powermate);
        if pmbutton == 1
            break;
        end
    end
elseif test_env
    disp('Click mouse to continue.')
    GetClicks;
end

%% Screen Setup
screens=Screen('Screens');
useScreen = 0;
p.screenWidthPixels = Screen('Rect', useScreen);
screenWidth = 53; %cm (testing room = 53cm)
viewingDistance = 110; %cm (testing room = 110cm)
visAngle = (2*atan2(screenWidth/2, viewingDistance))*(180/pi);

p.pixPerDeg = round(p.screenWidthPixels(3)/visAngle);
p.fixation=round(0.085 * p.pixPerDeg);
colors.grey = 128; colors.white = 255; colors.green = [0 255 0]; colors.blue = [0 0 255]; colors.black = [0 0 0]; colors.red = [220 20 60];
colors.dimgrey = [105 105 105]; colors.yellow = [255 255 0]; colors.magenta = [255 0 255]; colors.cyan = [0 255 255];

%% Grating Parameters

% Grating contrast for center
p.minContrast = 0.1;
p.maxContrast = 0.75;
p.numContrasts = 5;

% look at old data to choose which contrast to knock out (2 contrasts)
p.centerContrasts = [10.^linspace(log10(p.minContrast),log10(p.maxContrast),p.numContrasts)];
p.centerContrasts(1) = [];

% Grating Size
p.centerSize = round(2*p.pixPerDeg);
p.centerSizeDeg = p.centerSize/p.pixPerDeg;
p.centerRadius = round(p.centerSize/2 + p.pixPerDeg/2);

% Surround Size
p.surroundSize = p.screenWidthPixels(:,3);
p.surroundSizeDeg = p.surroundSize/p.pixPerDeg;
p.surroundRadius = round(p.surroundSize/2 + p.pixPerDeg/2);
if mod(p.surroundRadius,2) ~= 0
    p.surroundRadius = p.surroundRadius-1;
end
p.gapSize = round(0.15*p.pixPerDeg); %space between center and surround annulus (none in this version)

p.ecc = 10;
p.backgroundRadius = round(p.ecc * p.pixPerDeg); %radius of the background circle
p.eccentricity = round((p.ecc/2) * p.pixPerDeg); %distance from center of screen to center of target
p.innerRadius = p.eccentricity - (p.centerSize/2+p.gapSize); %inner edge of surround annulus
p.outerRadius = p.eccentricity + (p.centerSize/2+p.gapSize); %outer edge of surround annulus

% Fixation dot size
p.innerFixation = round(0.075*p.pixPerDeg);
p.outerFixation = round(0.75*p.pixPerDeg);

% Grating frequency, orientation, phase
p.orientation = 0;
p.stimorientation = 90;
p.freq = 1;
p.frequency_center = (2*p.centerRadius/p.pixPerDeg) * p.freq;
p.frequency_surround = (2*p.surroundRadius/p.pixPerDeg) * p.freq;
p.numPhases = 1;
p.centerPhase = 360;
p.surroundPhase = p.centerPhase;
%% Trial Events

% p.trialEvents = [centerContrast, centerLocation, probeContrast,
% probeLocation];
p.numOffsetLoc = 5; % 5 unique offsets to choose from including 0

p.stimConfigurations = 1:length(p.centerContrasts)*p.numOffsetLoc; %incorporate locations
[combs] = BalanceFactors(p.repetitions,0,p.stimConfigurations);

p.numTrials = length(combs);
p.numTrialsPerBlock = length(p.stimConfigurations);
p.numBlocks = p.numTrials/p.numTrialsPerBlock;

col1 = repmat(p.centerContrasts',p.numOffsetLoc,1); % center grating contrast
col1 = repmat(col1,p.repetitions,1);

%--------------------%
%              Location            %
%--------------------%
p.locSpacing = 10; p.locJit = 3;
p.probeLocWidth = 45;

%location is random every trial
col2 = randsample(0:p.locSpacing:359,length(col1),true)'; % center grating location
col2 = col2+randsample(-p.locJit:1:p.locJit,length(col1),true)';
col2 = round(col2); % make sure this is between 0-359

col2(col2>360) = col2(col2>360) - 360;
col2(col2<0) = col2(col2<0) +360;

probeOffset = [0 round(10.^linspace(0,log10(p.probeLocWidth),p.numOffsetLoc-1))]'; % probe grating location
probeOffset = repmat(probeOffset,length(p.centerContrasts),1);
probeOffset = repmat(probeOffset,p.repetitions,1);
probeOffset = probeOffset.*randsample([-1 1],length(probeOffset),true)';
col4 = col2+probeOffset;
col4(col4>360) = col4(col4>360) - 360;
col4(col4<0) = col4(col4<0) + 360;
col4 = round(col4);

%--------------------%
%              Contrast           %
%--------------------%

col3 = randsample(p.minContrast:0.01:p.maxContrast, length(col1),true)'; % probe grating contrast

% Integrating trial events

p.trialEvents = [col1 col2 col3 col4];
shuffled = 0;
p.trialEvents = Shuffle(p.trialEvents,2); shuffled = 1;

if shuffled
    disp('Trials have been shuffled.')
else
    disp('Trials have not been shuffled. Continue?')
    GetClicks;
end


%% TIMING PARAMETERS
% timing is in seconds
t.stimOn1 = 2; % stimulus 1 duration
t.retention = 1; % memory retention period
t.iti = 2;
t.startTime = 2;
t.responseTime = []; t.responseTime_est = 5;
t.flickerTime = 0.4;
t.flicker = 0.04;

t.trialDur = sum(t.stimOn1 + 2*t.flickerTime + t.responseTime_est + t.iti); % (s)
t.runDur = t.trialDur*size(p.trialEvents,1) + t.startTime*p.numBlocks; % (s)
if t.runDur/60 >= 1
    disp(['Total run duration: ' num2str(round(t.runDur/60)) ' min(s).'])
else
    disp(['Total run duration: ' num2str(round(t.runDur)) ' s.'])
end

%% CREATE STIMULI

%%Center
[xc,yc] = meshgrid((-p.centerRadius):(p.centerRadius)-1, (-p.centerRadius):(p.centerRadius)-1);
eccen = sqrt((xc).^2+(yc).^2); 	% calculate eccentricity of each point in grid relative to center of 2D image
centerGaussian = zeros(size(xc)); centerGaussian(eccen <= (p.centerSize/2)) = 1;
centerGaussian = conv2(centerGaussian,fspecial('gaussian',round(p.centerSize/10),p.centerSize/4), 'same');
centerTransparencyMask = zeros(size(xc)); centerTransparencyMask(eccen <= (p.centerSize/2)+p.gapSize/2)=255;

%%Surround
[xs,ys] = meshgrid((-p.surroundRadius):(p.surroundRadius)-1, (-p.surroundRadius):(p.surroundRadius)-1);
eccen = sqrt((xs).^2+(ys).^2); 	% calculate eccentricity of each point in grid relative to center of 2D image
Annulus = zeros(size(xs)); Annulus(eccen <= p.innerRadius) = 1;
Annulus(eccen >= p.outerRadius) = 1;
Annulus = conv2(Annulus, fspecial('gaussian',round(p.centerSize/10),p.centerSize/4),'same');

%%Background
bgAnnulus = zeros(size(Annulus));
bgAnnulus(eccen <= p.backgroundRadius) = 1;
bgAnnulus(eccen >= p.backgroundRadius) = 0;

% Make unique grating
[Xc,Yc] = meshgrid(0:(size(centerGaussian,1)-1),0:(size(centerGaussian,1)-1));
[Xs,Ys] = meshgrid(0:(size(Annulus,1)-1),0:(size(Annulus,1)-1));

centerGrating = NaN(size(centerGaussian));
surroundGrating = NaN(size(Annulus));

centerPatch = (sin(p.frequency_center*2*pi/size(centerGrating,1)*(Xc.*sin(p.orientation*(pi/180))+Yc.*cos(p.orientation*(pi/180)))-p.centerPhase));
centerGrating = (centerPatch .* centerGaussian);

surroundGrating = (sin(p.frequency_surround*2*pi/size(surroundGrating,1)*(Xs.*sin(p.orientation*(pi/180))+Ys.*cos(p.orientation*(pi/180)))-p.surroundPhase));
tmpSurround = (surroundGrating .* Annulus);
tmpSurround(bgAnnulus == 0) = -1;
surroundGrating = tmpSurround;
%% MASK

%sf filter
cutoff = [0.5 1.5] .*round(size(Annulus,1)/p.pixPerDeg);
f = freqspace(size(Annulus,1));

%bandpass filter
filter = Bandpass2(size(Annulus,1), f(cutoff(1)), f(cutoff(2)));

%low pass filter
h = fspecial('disk', 5);
filter = conv2(filter, h, 'same');
filter = filter/max(filter(:));

%noise
maskGrating = NaN(round(t.flickerTime/t.flicker), size(Annulus,1), size(Annulus,1));
for n = 1:(t.flickerTime/t.flicker)
    noise = -1+2.*rand(size(Annulus,1));
    fftNoise = fftshift(fft2(noise));
    filterNoise = fftNoise .* filter;
    newNoise = real(ifft2(fftshift(filterNoise)));
    noiseMask = newNoise./(max(abs(newNoise(:))));
    tempMask = (noiseMask.*bgAnnulus);
    tempMask(bgAnnulus == 0) = -1; %to make black background make == 0
    maskGrating (n,:,:) = tempMask;
end

%% WINDOW SETUP

[window,rect] = Screen('OpenWindow', useScreen, colors.black, []); %test screen size [0 0 700 500]

OriginalCLUT = Screen('ReadNormalizedGammaTable', window);
if ~test_env
    load('linearizedCLUT.mat'); %testing room 208 only
    Screen('LoadNormalizedGammaTable',window, linearizedCLUT); %testing room 208 only
end

% Enable alpha blending
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Define coordinates where to draw the fixation
CenterX = rect(3)/2; CenterY = rect(4)/2;
% Coordinates for location on left and right side of fixation
center = [CenterX CenterY];

Screen('TextStyle', window, 1);
Screen('TextSize', window, 16);
t.ifi = Screen('GetFlipInterval',window); % grab screen refresh rate

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

if ~test_env
    % Welcome Screen
    Screen('TextStyle', window, 1);
    Screen('TextSize', window, 14);
    welcomeText = ['On each trial, you will see a center grating presented at a random location in the periphery.' '\n'...
        'At the end of a trial, a probe grating will appear at a random location.' '\n'...
        'You will then manipulate this probe grating to reconstruct the contrast of the target grating you perceived.' '\n'...
        'Rotating the powermate will manipulate the contrast of this probe.' '\n'... 
        'Clicking the powermate will lock in your contrast estimate and the next trial will begin after your response.' '\n'...
        'Remember to keep your eyes at fixation throughout the entirety of a trial!' '\n'...
        'Additionally, breaks will be provided throughout the session.' '\n'...
        'Click the powermate to start experiment.'];
    DrawFormattedText(window, welcomeText, 'center', 'center', 255);
    Screen('Flip', window);
end


% Check powermate works by forcing a click
while 1
    [pmbutton, ~] = PsychPowerMate('Get', powermate);
    if pmbutton == 1
        break;
    end
end
% Starting Screen
Screen('FillOval', window, colors.grey, [CenterX-p.backgroundRadius CenterY-p.backgroundRadius CenterX+p.backgroundRadius CenterY+p.backgroundRadius]);
Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation]);
Screen('FillOval' , window, colors.green, [CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
Screen('Flip',window);
WaitSecs(t.startTime);


nBlock = 1; % #sets completed tracker
for nTrial = 1:size(p.trialEvents,1)
    %--------------------%
    %         Trial Settings        %
    %--------------------%
    p.centerContrast = p.trialEvents(nTrial,1);
    p.currentDegreeLocation = round(p.trialEvents(nTrial,2));
    targetXY = round([CenterX+p.eccentricity*cos(p.currentDegreeLocation*(pi/180))' CenterY-p.eccentricity*sin(p.currentDegreeLocation*(pi/180))']);
    
    probeContrast = p.trialEvents(nTrial,3);
    probeLocation = p.trialEvents(nTrial,4);
    ProbeXY = round([CenterX + p.eccentricity*(cos(probeLocation*(pi/180)))' CenterY - p.eccentricity*(sin(probeLocation*(pi/180)))']);
    
    %--------------------%
    %           Stimulus 1           %
    %--------------------%
    % Draw center stimulus
    centerTexture = (centerGrating*p.centerContrast)*colors.grey + colors.grey;
    centerTexture(:,:,2) = centerTransparencyMask;
    centerStimulus = Screen('MakeTexture', window, centerTexture);
    Screen('FillOval', window, colors.grey, [CenterX-p.backgroundRadius CenterY-p.backgroundRadius CenterX+p.backgroundRadius CenterY+p.backgroundRadius]);
    Screen('DrawTexture', window, centerStimulus, [], ...
        CenterRectOnPoint([0 0 size(centerTexture,1) size(centerTexture,1)], targetXY(1), targetXY(2)), p.stimorientation);
    % Draw fixation
    Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation])
    Screen('FillOval', window,colors.green,[CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
    Screen('Flip',window);
    WaitSecs(t.stimOn1);
    if test_env
        GetClicks;
    end
    
    %--------------------%
    %               Mask 1             %
    %--------------------%
    indx = Shuffle(1:round(t.flickerTime/t.flicker));
    StartMask = GetSecs;
    for a = indx
        if GetSecs > StartMask + t.flickerTime - t.flicker
            break
        end
        Screen('DrawTexture', window, Mask(a,:), [],CenterRectOnPoint([0 0 size(surroundGrating,1) size(surroundGrating,1)], CenterX, CenterY));
        Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation]);
        Screen('FillOval', window, colors.green, [CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
        Screen('Flip', window);
        WaitSecs(t.flicker);
    end
    
    %--------------------%
    %          Retention             %
    %--------------------%
    %%% Retention interval 1
    Screen('FillOval', window, colors.grey, [CenterX-p.backgroundRadius CenterY-p.backgroundRadius CenterX+p.backgroundRadius CenterY+p.backgroundRadius]);
    % Draw fixation
    Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation])
    Screen('FillOval', window,colors.green,[CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
    Screen('Flip',window);
    WaitSecs(2*t.retention);
    
    %--------------------%
    %               Mask 2             %
    %--------------------%
    indx = Shuffle(1:round(t.flickerTime/t.flicker));
    StartMask = GetSecs;
    for a = indx
        if GetSecs > StartMask + t.flickerTime - t.flicker
            break
        end
        Screen('DrawTexture', window, Mask(a,:), [],CenterRectOnPoint([0 0 size(surroundGrating,1) size(surroundGrating,1)], CenterX, CenterY));
        Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation]);
        Screen('FillOval', window, colors.green, [CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
        Screen('Flip', window);
        WaitSecs(t.flicker);
    end
    
    %--------------------%
    %               Probe               %
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
    
    Probe = Screen('MakeTexture', window, squeeze(centerGrating)* (probeContrast *colors.grey) + colors.grey);
    Screen('DrawTexture', window, Probe, [], CenterRectOnPoint([0 0 size(centerTexture,1) size(centerTexture,1)], ProbeXY(1), ProbeXY(2)), p.stimorientation);
    Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation]);
    Screen('FillOval', window, colors.green, [CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
    Screen('Flip', window);
    
    startResponseTime = GetSecs; % get the start time of response
    if test_env
        GetClicks;
    else
        
        %%% CONTRAST REPORT %%%
        % Allow for dial rotation for contrast update
        [~, contrastangle] = PsychPowerMate('Get', powermate);
        while 1 %start inf loop
            % % Query PowerMate button state and rotation angle in "clicks"
            [pmbutton_contrast, angle2] = PsychPowerMate('Get', powermate);
            % % 1st button is the "or" of the 1st mouse button and the actual PowerMate button
            if contrastangle ~= angle2
                
                % % Convert turn of dial first to degrees and then to contrast:
                angles = (angle2 * 3.8298)/360;
                angles = ((contrastangle-angle2)*3.8298);
                changecontrast = angles/360;
                probeContrast = probeContrast - changecontrast; % update the contrast relative to last dial position
                % % Make sure we stay in range
                
                if probeContrast > 1
                    probeContrast = 1;
                elseif probeContrast < 0
                    probeContrast = 0.001;
                end
                Screen('FillOval', window, colors.grey, [CenterX-p.backgroundRadius CenterY-p.backgroundRadius CenterX+p.backgroundRadius CenterY+p.backgroundRadius]);
                
                Target = Screen('MakeTexture', window, squeeze(centerGrating)* (probeContrast*colors.grey) + colors.grey);
                Screen('DrawTexture', window, Target, [], CenterRectOnPoint([0 0 size(centerTexture,1) size(centerTexture,1)], ProbeXY(1), ProbeXY(2)), p.stimorientation)
                Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation]);
                Screen('FillOval', window, colors.green, [CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation])
                Screen('Flip', window);
                
                contrastangle = angle2;
            end
            if pmbutton_contrast == 1
                contrastTime = GetSecs;
                estimatedContrast(nTrial) = probeContrast;
                differenceContrast(nTrial) = p.trialEvents(nTrial,3) - estimatedContrast(nTrial);
                pmbutton_contrast = 0;
                break
            end
        end
       
    end
    
    responseTime(nTrial) = (GetSecs-startResponseTime);
    
    % ITI
    Screen('FillOval', window, colors.grey, [CenterX-p.backgroundRadius CenterY-p.backgroundRadius CenterX+p.backgroundRadius CenterY+p.backgroundRadius]);
    Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation]);
    Screen('FillOval' , window, colors.green, [CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
    Screen('Flip',window);
    WaitSecs(t.iti);
    
    
    %--------------------%
    %               Break               %
    %--------------------%
    if mod(nTrial,p.numTrialsPerBlock) == 0
        rest = GetSecs;
        Screen('FillOval', window, colors.grey, [CenterX-p.backgroundRadius CenterY-p.backgroundRadius CenterX+p.backgroundRadius CenterY+p.backgroundRadius]);
        Screen('TextStyle', window, 1);
        Screen('TextSize', window, 16);
        breakText = ['You make take a short break now.\n Or press the powermate to continue.\n'];
        DrawFormattedText(window, breakText, 'center', 'center', colors.white);
        Screen('Flip', window);
        restText = ['You can take a short break now, ' '' '\n' ...
            'or press the dial to continue' '\n' '\n' ...
            num2str(nBlock) '/' num2str(p.numBlocks) ' completed.' ];
        DrawFormattedText(window, restText, 'center', 'center', colors.white);
        Screen('Flip', window);
        pmbuttonbreak = 0;
        WaitSecs(1);
        while 1
            [pmbuttonbreak, a] = PsychPowerMate('Get', powermate);
            if pmbuttonbreak == 1
                break;
            end
        end
        t.restTime(nBlock) = (GetSecs-rest)/60;
        nBlock = nBlock + 1;     
    end
    
end

%% SAVE OUT THE DATA FILE

if p.repetitions > 5 && ~test_env
    data.EstimatedContrast = estimatedContrast;
    data.DifferenceContrast = differenceContrast;
    data.responseTime = responseTime;
    
    cd(dataDir);
    theData(p.runNumber).t = t;
    theData(p.runNumber).p = p;
    theData(p.runNumber).data = data;
    save(fileName, 'theData')
    cd(expDir);
end

sca;
ShowCursor;