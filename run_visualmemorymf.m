% Gratings
close all; clear all; clc;
Screen('Preference', 'SkipSyncTests', 1);

%% PREPARE
p.repetitions = 1; % maybe 20?
p.numBlocks = p.repetitions;
%% SCREEN PARAMTERS
screens=Screen('Screens');
useScreen=max(screens);
p.screenWidthPixels = Screen('Rect', useScreen);
screenWidth = 51; %cm
viewingDistance = 50; %cm
visAngle = (2*atan2(screenWidth/2, viewingDistance))*(180/pi);


p.pixPerDeg = round(p.screenWidthPixels(3)/visAngle);
p.fixation=round(0.085 * p.pixPerDeg);
p.Repetitions = 1;
colors.grey = 128; colors.white = 255; colors.green = [0 255 0]; colors.blue = [0 0 255]; colors.black = [0 0 0]; colors.red = [220 20 60]; 
colors.dimgrey = [105 105 105]; colors.yellow = [255 255 0]; colors.magenta = [255 0 255]; colors.cyan = [0 255 255];


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
% Controls the brightness of the powermate color
PsychPowerMate('SetBrightness', powermate, 20);

while 1
    [pmbutton, ~] = PsychPowerMate('Get', powermate);
    if pmbutton == 1;
        break;
    end
end
%% TIMING PARAMETERS
times.stimon = 1;     % in sec
times.retention = 2;    % Different retention intervalst.iti = 0.3;
times.iti = 2;
times.Starttime = 2;
times.ResponseTime = [];
times.flickertime = 0.4; %0.4
times.flicker = 0.04; %0.04

%% GRATING PARAMETERS

% grating contrast for center and surround
p.minContrast = 0.1;
p.maxContrast = 0.75;
p.numContrasts = 5;
%p.centerContrast = [0 10.^linspace(log10(p.minContrast),log10(p.maxContrast),p.numContrasts)];
p.centerContrast = [10.^linspace(log10(p.minContrast),log10(p.maxContrast),p.numContrasts)];
p.surroundContrast = 1;

% Size Parameters
p.centerSize = round(2*p.pixPerDeg);
p.surroundSize = p.screenWidthPixels(:,3);
p.gapSize = round(0.08*p.pixPerDeg);

p.backgroundRadius = p.screenWidthPixels(:,4)/2; %radius of the background circle
p.eccentricity = p.backgroundRadius/4; %distance from center of screen to center of target
p.innerRadius = p.eccentricity - (p.centerSize/2+p.gapSize); %inner edge of surround annulus
p.outerRadius = p.eccentricity + (p.centerSize/2+p.gapSize); %outer edge of surround annulus

p.innerFixation = round(0.075*p.pixPerDeg);
p.outerFixation = round(0.75*p.pixPerDeg);

% Frequency, Orientation, Phase
p.orientation = 0;
p.stimorientation = 90;
p.freq = 2;
p.frequency_center = p.centerSize/p.pixPerDeg * p.freq;
p.frequency_surround = p.surroundSize/p.pixPerDeg * p.freq;
p.numPhases = 1;
p.centerPhase = 360;
p.surroundPhase = p.centerPhase;


%% TRIAL EVENTS %%

%% WINDOW SETUP
[window,rect] = Screen('OpenWindow', useScreen, colors.black, []);
OriginalCLUT = Screen('ReadNormalizedGammaTable', window);

%%% Conditions %%%

%we want 20 repetitions for per condition
p.stimConfigurations = 1:length(p.centerContrast)*2; % 2 because there's perception and memory condition

[col1] = BalanceFactors(p.numBlocks,0,p.stimConfigurations);

% 1 - perception: center and surround, mask, blank, mask
% 2 - working memory: center, mask, surround, mask
col1(1:5) = 1; %half of the events are perception trials
col1(6:10) = 2; %half of the events are working memory trials

%conditions are perception and working memory, hence the first colum of the
%trial events matrix
% p.trialEvents(:,end+1) = col2;

%%% Location %%%
%location is random every time
%incorporate this into whichLocation

col2 = randi(360,1,10)';
col4 =  randi(360,1,10)';
%want to first probably generate a function that selects a random 10
%locations (5 trials per condition), and then have them index into this
%matrix as a location factorn (as the second column)

%%% Contrast %%%
%5 different contrast levels
p.centerContrast = [p.centerContrast p.centerContrast];
p.centerContrast = p.centerContrast';
col3 = p.centerContrast;
col5 = rand(length(col3),1);

% bring all 3 together
p.trialEvents = [col1 col2 col3 col4 col5]; %[condition targetLocation targetContrast probeLocation probeContrast]
% each row, 1-10 represents a different trial event.

% create loop that goes through all different possible trials
[trials, parameters] = size(p.trialEvents);
trials = [1:trials];
config1Shuffle = Shuffle(trials(1:5));
config2Shuffle = Shuffle(trials(6:10));
%if we're going to be testing the configurations on different days find a
%way to make this an input of what trial we want to run or something
p.trialOrder = [config1Shuffle config2Shuffle];
%creates an array thast shuffles the trial numbers, and the code then


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
maskGrating = NaN(round(times.flickertime/times.flicker), p.surroundSize, p.surroundSize);
for nTrial = 1:(times.flickertime/times.flicker)
    noise = -1+2.*rand(p.surroundSize);
    fftNoise = fftshift(fft2(noise));
    filterNoise = fftNoise .* filter;
    newNoise = real(ifft2(fftshift(filterNoise)));
    noiseMask = newNoise./(max(abs(newNoise(:))));
    tempMask = (noiseMask.*bgAnnulus);
    tempMask(bgAnnulus == 0) = -1; %to make black background make == 0
    maskGrating (nTrial,:,:) = tempMask;
end 

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
%% Window
% Esc to quit
StartKey=zeros(1,256); StartKey(KbName({'ESCAPE'})) = 1;
PsychHID('KbQueueCreate', deviceNumber, StartKey);
PsychHID('KbQueueStart', deviceNumber);


% Enable alpha blending
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Define coordinates where to draw the fixation
CenterX = rect(3)/2; CenterY = rect(4)/2;
% Coordinates for location on left and right side of fixation
center =  [CenterX CenterY];

% Center grating coordinates
% circDeg = Shuffle(1:360);
% + or - (p.outerRadius-p.innerRadius)/2
% targetXY = round([CenterX+p.eccentricity*cos(circDeg*(pi/180))' CenterY-p.eccentricity*sin(circDeg*(pi/180))']);

%% EXPERIMENT LOOP
% Starting Screen
Screen('FillOval', window, colors.grey,  [CenterX-p.backgroundRadius CenterY-p.backgroundRadius CenterX+p.backgroundRadius CenterY+p.backgroundRadius]);
Screen('FillOval', window, colors.black,  [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation]);
Screen('FillOval' , window, colors.green, [CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
Screen('Flip',window);
WaitSecs(times.Starttime);

for nTrial = 1:length(trials)
    disp(['Start of trial ' num2str(nTrial)])
    p.currentDegreeLocation = round(p.trialEvents(nTrial,2));
    p.centerContrast = p.trialEvents(nTrial,3);
    targetXY = round([CenterX+p.eccentricity*cos(p.currentDegreeLocation*(pi/180))' CenterY-p.eccentricity*sin(p.currentDegreeLocation*(pi/180))']);

    if p.trialEvents(nTrial,1) == 1 %perception 
    % Surround Grating
    surroundTexture = (surroundGrating*p.surroundContrast)*colors.grey + colors.grey;
    surroundTexture(:,:,2) = surroundTransparencyMask;
    surroundStimulus = Screen('MakeTexture', window, surroundTexture);
    Screen('DrawTexture', window, surroundStimulus, [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], CenterX, CenterY), p.stimorientation);
    else
        Screen('FillOval', window, colors.grey,  [CenterX-p.backgroundRadius CenterY-p.backgroundRadius CenterX+p.backgroundRadius CenterY+p.backgroundRadius]);
    end
    
    % Center Grating
    centerTexture = (centerGrating*p.centerContrast)*colors.grey + colors.grey;
    centerTexture(:,:,2) = centerTransparencyMask;
    centerStimulus = Screen('MakeTexture', window, centerTexture);
    % Screen('DrawTexture', window, centerStimulus, [], ...
    %         CenterRectOnPoint([0 0 p.centerSize p.centerSize], targetXY(p.currentDegreeLocation,1), targetXY(p.currentDegreeLocation,2)), p.stimorientation);

    Screen('DrawTexture', window, centerStimulus, [], ...
            CenterRectOnPoint([0 0 p.centerSize p.centerSize], targetXY(1), targetXY(2)), p.stimorientation);


    % Fixation
    Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation])
    Screen('FillOval', window,colors.green,[CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
    Screen('Flip',window);
    WaitSecs(times.stimon);
    
    % Texture
    Mask = NaN(round(times.flickertime/times.flicker),1);
    for m = 1:round(times.flickertime/times.flicker)
        Mask(m,:) = Screen('MakeTexture', window, squeeze(maskGrating(m,:,:))* colors.grey + colors.grey);
    end


    %Mask 1 appears
    indx = Shuffle(1:round(times.flickertime/times.flicker));
    StartMask = GetSecs;
    for a = indx
        if GetSecs > StartMask + times.flickertime - times.flicker
            break
        end        
        Screen('DrawTexture', window, Mask(a,:), [],CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], CenterX, CenterY));
        Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation]);
        Screen('FillOval', window, colors.green, [CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
        Screen('Flip', window);
        WaitSecs(times.flicker);        
    end
    
    Screen('FillOval', window, colors.grey,  [CenterX-p.backgroundRadius CenterY-p.backgroundRadius CenterX+p.backgroundRadius CenterY+p.backgroundRadius]);
    
    if p.trialEvents(nTrial,1) == 2
        % Surround Grating
        surroundTexture = (surroundGrating*p.surroundContrast)*colors.grey + colors.grey;
        surroundTexture(:,:,2) = surroundTransparencyMask;
        surroundStimulus = Screen('MakeTexture', window, surroundTexture);
        Screen('DrawTexture', window, surroundStimulus, [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], CenterX, CenterY), p.stimorientation);
    end

    % Fixation
    Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation])
    Screen('FillOval', window,colors.green,[CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
    Screen('Flip',window);
    WaitSecs(times.retention);


    %Mask 2 appears
    indx = Shuffle(1:round(times.flickertime/times.flicker));
    StartMask = GetSecs;
    for a = indx
        if GetSecs > StartMask + times.flickertime - times.flicker
            break
        end        
        Screen('DrawTexture', window, Mask(a,:), [],CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], CenterX, CenterY));
        Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation]);
        Screen('FillOval', window, colors.green, [CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
        Screen('Flip', window);
        WaitSecs(times.flicker);        
    end

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
    [~, startangle] = PsychPowerMate('Get', powermate);

    while 1 %start inf loop
%         % Query PowerMate button state and rotation angle in "clicks"
        [pmbutton, angle] = PsychPowerMate('Get', powermate);
%         % 1st button is the "or" of the 1st mouse button and the actual PowerMate button
        if startangle ~= angle
            
%             % Convert turn of dial first to degrees and then to contrast:
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
%             % make sure angle stays in 0-360 range
            correctedAngle = mod(initialAngle, 360);
            
            data.EstimatedLocation(nTrial) = correctedAngle;
            
%             %make sure difference is in the 180 range
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
    
    % PowerMate is sampled at 10msec intervals, therefore have a short
    % break to make sure it doesn't skip the contrast task
    WaitSecs(0.2);
    
    %Button press for contrast
    [~, contrastangle] = PsychPowerMate('Get', powermate);
    while 1 %start inf loop
        % Query PowerMate button state and rotation angle in "clicks"
        [pmbutton_contrast, angle2] = PsychPowerMate('Get', powermate);
        % 1st button is the "or" of the 1st mouse button and the actual PowerMate button
        if contrastangle ~= angle2
            
            % Convert turn of dial first to degrees and then to contrast:
            %             angles = (angle * 3.8298)/360;
            angles = ((contrastangle-angle2)*3.8298);
            changecontrast = angles/360;
            intial_contrast = intial_contrast - changecontrast; % update the contrast relative to last dial position
            % Make sure we stay in range
            
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
    
disp(['End of trial ' num2str(nTrial)])    
%GetClicks;      
WaitSecs(times.iti);
end

%Draw some more text to the screen outside of the loop:
Screen(window,'TextSize',30);
ByebyeText = ['Great work! You have completed this run.' '\n' '\n' ...
    'Please let the experimenter know you have finished.'];
Screen('FillOval', window, colors.grey, [CenterX-p.backgroundRadius CenterY-p.backgroundRadius CenterX+p.backgroundRadius CenterY+p.backgroundRadius]);
DrawFormattedText(window, ByebyeText, 'center', 'center', [255 255 255]);
Screen('Flip', window);
WaitSecs(3);

Screen('LoadNormalizedGammaTable', window, OriginalCLUT);
Screen('CloseAll')