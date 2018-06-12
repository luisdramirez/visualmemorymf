%% PREPARE AND COLLECT INFO
close all
clear all

KbName('UnifyKeyNames');
Screen('Preference', 'SkipSyncTests', 1);

commandwindow;

p.Subject = 'tmp';

p.Repetitions = 1   ;

% Check which devicenumber the powermate is assigned to
% powermate = PsychPowerMate('Open');
% if isempty(powermate)
%     error('problem with the powermate');
% end
% Controls the brightness of the powermate color
% PsychPowerMate('SetBrightness', powermate, 20);

% Check which devicenumber the keyboard is assigned to
deviceNumber = 0;
[keyBoardIndices, ProductNames] = GetKeyboardIndices;
deviceString = 'Apple Internal Keyboard / Trackpad';
% deviceString = 'Wired USB Keyboard';
%deviceString = 'Apple Keyboard';
% deviceString = 'USB Keyboard';
for i = 1:length(ProductNames)
    if strcmp(ProductNames{i}, deviceString)
        deviceNumber = keyBoardIndices(i);
        break;
    end
end
if deviceNumber == 0
    error('No device by that name was detected');
end

%set directory
expdir = pwd; %Set the experimental directory to the current directory 'pwd'
datadir = 'Data'; %Set the path to a directory called 'Data'
t.MySeed = sum(100*clock);
rng(t.MySeed); % make sure we start with a random seed
t.TheDate = datestr(now,'yymmdd'); %Collect todays date
t.TimeStamp = datestr(now,'HHMM'); %Timestamp

%cd(datadir);
if exist(['surrSuppression_spatial_contrast_est_', p.Subject, '.mat'],'file');
    load(['surrSuppression_spatial_contrast_est_', p.Subject, '.mat']);
    runnumber = length(TheData)+1;
else
    runnumber = 1;
end
cd(expdir);


%% SCREEN PARAMETERS
Screens = Screen('Screens'); % look at available screens
p.ScreenWidthPixels = Screen('Rect', Screens(1));
ScreenWidth = 42; % 29 cm macbook air, 40 cm trinitron crt, 60 cm Qnix screen
ViewDistance = 128; % in cm, ideal distance: 1 cm equals 1 visual degree
VisAngle = (2*atan2(ScreenWidth/2, ViewDistance))*(180/pi); % Visual angle of the whole screen
p.ppd = round(p.ScreenWidthPixels(3)/VisAngle); % pixels per degree visual angle
p.Grey = 128;
p.Black = 0;


%% DEFINE MAIN PARAMETERS
%%Main Parameters
p.numContrasts = 5;
p.testContrasts = 10.^linspace(log10(0.1),log10(0.75),p.numContrasts);
p.surroundContrast = 1; % display surround at full contrast
p.stimConfigurations = 2;
p.numTrials = p.numContrasts * p.stimConfigurations *p.Repetitions; % multiple of locations and possible targets

%Size Parameters
p.PatchSize = round(1*p.ppd); %Size of the patch that is drawn on various screen locations to fit gratings into (in pixels):
p.Surround = p.ScreenWidthPixels(:,3); %Length of the screen
p.GapSize = round(0.08 * p.ppd); %size of gap between radii and stimulus
backgroundRadius = p.ScreenWidthPixels(:,4)/2; %radius of the background circle
p.eccentricity = round(backgroundRadius/2); %distance from center of screen to center of target
innerRadius = p.eccentricity - (p.PatchSize/2+p.GapSize); %inner edge of surround annulus
outerRadius = p.eccentricity + (p.PatchSize/2+p.GapSize); %outer edge of surround annulus

p.InnerFixation = round(0.05*p.ppd); %green fixation point
p.OuterFixation = round(0.5*p.ppd); %black patch around fixation point

p.orientation = 0;
p.stimorientation = 90; %orientation of surround and center
freq = 2; %Frequency of the grating in cycles per degree
p.Freq = p.PatchSize/p.ppd * freq;
p.Freq_background = p.Surround/p.ppd * freq;
p.phase = randsample(1:180,p.numTrials*3, true);
p.phase = reshape(p.phase, [p.numTrials 3]);


% TIMING PARAMETERS
t.stimon = 5;     % in sec
t.retention = 2;    % Different retention intervalst.iti = 0.3;
t.iti = 0.5;
t.Starttime = 2;
t.ResponseTime = [];
t.flickertime = 0.4;
t.flicker = 0.04;

%% TRIAL EVENTS

% stimulus can either be collinear (2)or no surround (1)
whichConfigurations = repmat([ones(1,p.numContrasts) 2*ones(1,p.numContrasts)]', [p.numTrials/(p.stimConfigurations*p.numContrasts) 1]);
% Determine different test contrasts
whichContrast = repmat(p.testContrasts', [p.numTrials/p.numContrasts 1]);
% Angle of location
whichLocation = repmat(360, [p.numTrials,1]); %randsample(1:360,p.numTrials, true)';
% whichLocation(1:3,1) = 360;
% whichLocation(4:6,1) = 300;
% whichLocation(7:9,1) = 180;
%Probe contrast
probecontrast = randsample(0.1:0.01:0.9, p.numTrials, true)';
%probe location
probelocation = randsample(1:360,p.numTrials, true)';
% probelocation(1:3,1) = 360;
% probelocation(4:6,1) = 300;
% probelocation(7:9,1) = 180;
% Orientation
% combine all possible trials together in p.TrialEvents
% 1: stimulus configuration; 2:test contrast;  3: probe contrast(random);
% 4: location angle 5:probe location (randon)
tempTrialEvents = [whichConfigurations whichContrast(:) probecontrast whichLocation probelocation];
% shuffle all the rows, but keep the location in this order, so that conditions are balanced within within location from fixation
p.TrialEvents = flipud(tempTrialEvents); %Shuffle(tempTrialEvents,2);

%% CREATE STIMULI

%%Patch
[x,y] = meshgrid((-p.PatchSize/2):(p.PatchSize/2)-1, (-p.PatchSize/2):(p.PatchSize/2)-1);
eccen = sqrt((x).^2+(y).^2); 	% calculate eccentricity of each point in grid relative to center of 2D image
centerGaussian = zeros(p.PatchSize); centerGaussian(eccen <= (p.PatchSize/2)) = 1;

%%Surround
[x,y] = meshgrid((-p.Surround/2):(p.Surround/2)-1, (-p.Surround/2):(p.Surround/2)-1);
eccen = sqrt((x).^2+(y).^2); 	% calculate eccentricity of each point in grid relative to center of 2D image
Annulus = zeros(p.Surround); Annulus(eccen <= innerRadius) = 1;
Annulus(eccen >= outerRadius) = 1;

%%Background
bgAnnulus = zeros(p.Surround); 
bgAnnulus(eccen <= backgroundRadius) = 1;
bgAnnulus(eccen >= backgroundRadius) = 0;

%% Mask
% Create SF filter
Cutoff = [1 3] .*round(p.Surround/p.ppd);
f = freqspace(p.Surround);
% bandpass filter
Filter = Bandpass2(p.Surround, f(Cutoff(1)), f(Cutoff(2)));
% low pass filter
% Filter = ones(p.Surround); Filter(eccen > Cutoff(2)) = 0;
h = fspecial('disk', 5);
Filter = conv2(Filter, h, 'same');
Filter = Filter/max(Filter(:));
% Make noise
mask_grating = NaN(round(t.flickertime/t.flicker),p.Surround,p.Surround);
for ii = 1:(t.flickertime/t.flicker)
    noise = -1 + 2.*rand(p.Surround);
    fft_noise = fftshift(fft2(noise));
    filt_noise = fft_noise .* Filter;
    new_noise = real(ifft2(fftshift(filt_noise)));
    noise_Mask = new_noise./(max(abs(new_noise(:))));
    tmp_mask = (noise_Mask.*bgAnnulus);
    tmp_mask(bgAnnulus == 0) = -1; %to make black background make == 0
    mask_grating(ii,:,:) = tmp_mask;
end



% Make unique grating for every trial
[Xc,Yc] = meshgrid(0:(p.PatchSize-1),0:(p.PatchSize-1));
[Xs,Ys] = meshgrid(0:(p.Surround-1),0:(p.Surround-1));

patch_grating = NaN(p.numTrials*2,p.PatchSize,p.PatchSize);
background_grating = NaN(p.numTrials*2,p.Surround,p.Surround);
patch_target = NaN(p.numTrials*2,p.PatchSize,p.PatchSize);

for s = 1:p.numTrials
    targetpatch = (sin(p.Freq*2*pi/p.PatchSize*(Xc.*sin(p.orientation*(pi/180))+Yc.*cos(p.orientation*(pi/180)))-p.phase(s,1)));
    patch_grating(s,:,:) = (targetpatch .* centerGaussian);
    
    targetsurround = (sin(p.Freq_background*2*pi/p.Surround*(Xs.*sin(p.orientation*(pi/180))+Ys.*cos(p.orientation*(pi/180)))-p.phase(s,2)));
    tmp_background = (targetsurround .* Annulus);
    tmp_background(bgAnnulus == 0) = -1;
    background_grating(s,:,:) = tmp_background;
    
    target = (sin(p.Freq*2*pi/p.PatchSize*(Xc.*sin(p.orientation*(pi/180))+Yc.*cos(p.orientation*(pi/180)))-p.phase(s,3)));
    patch_target(s,:,:) = (target .* centerGaussian);
    
end

%% Window Setup
[window,rect] = Screen('OpenWindow', Screens(1), p.Black, [0 0 600 400]);
OriginalCLUT = Screen('ReadNormalizedGammaTable', window);
% load('linearizedCLUT.mat');
% Screen('LoadNormalizedGammaTable', window, linearizedCLUT);
%HideCursor;
white = 255; green = [0 255 0];

% Enable alpha blending
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


% Define coordinates where to draw the fixation
CenterX = rect(3)/2; CenterY = rect(4)/2;
% Coordinates for location on left and right side of fixation
Patch =  [CenterX CenterY];

%Make textures early
Mask = NaN(round(t.flickertime/t.flicker),1);
for m = 1:round(t.flickertime/t.flicker)  
    Mask(m,:) = Screen('MakeTexture', window, squeeze(mask_grating(m,:,:))* p.Grey +p.Grey );
end

%% START THE EXPERIMENT
% Draw some text to the screen first outside of the experimental loop:

% experiment setup
space = zeros(1,256); space(KbName('Space')) = 1;
PsychHID('KbQueueCreate',deviceNumber, space);
PsychHID('KbQueueStart', deviceNumber);

Screen('FillOval', window, p.Grey, [CenterX-backgroundRadius CenterY-backgroundRadius CenterX+backgroundRadius CenterY+backgroundRadius]);

Screen('TextStyle', window, 1);
Screen('TextSize', window, 16);
WelcomeText = ['Welcome!' '\n' '\n'...
   '' '\n' '\n'...
    'On every trial a stimulus is presented at a random location, and varies in its intensity.' '\n' '\n' ...
    'Your task is to maintain a precise representation of the location and intensity of the stimulus. ' '\n' '\n' ...
    'After a 2 second delay period you are asked to dial the probe to match first, the location of the stimulus, ' '\n' '\n' ...
    'and then the intensity of the stimulus to the one you have in memory as closely as possible.' '\n' '\n' ...
    'At the start and end of the delay period a flickering mask will appear, ' '\n' '\n' ...
    'cueing you to dial the probe to make your response. ' '\n' '\n' ...
    '' '\n' '\n'...
    'Be sure to always maintain steady fixation on the green dot! ' '\n' '\n' '\n' ...
    'Click the dial to continue.' '\n' '\n' ];
DrawFormattedText(window, WelcomeText, 'center', 'center', 255);
Screen('Flip', window);

% while 1
%     [pmbutton, ~] = PsychPowerMate('Get', powermate);
%     if pmbutton == 1;
%         break;
%     end
% end


%Starting Screen
Screen('FillOval', window, p.Grey, [CenterX-backgroundRadius CenterY-backgroundRadius CenterX+backgroundRadius CenterY+backgroundRadius]);
Screen('FillOval', window, p.Black, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation]);
Screen('FillOval', window, green, [CenterX-p.InnerFixation CenterY-p.InnerFixation CenterX+p.InnerFixation CenterY+p.InnerFixation]);
Screen('Flip', window);

StartTime = GetSecs;
WaitSecs(t.Starttime);

% Make sure we can press esc to quit the experiment
StartKey=zeros(1,256); StartKey(KbName({'ESCAPE'})) = 1;
PsychHID('KbQueueCreate', deviceNumber, StartKey);
PsychHID('KbQueueStart', deviceNumber);

% Preallocate some variables
data.EstimatedContrast = NaN(p.numTrials, 1);
data.DifferenceContrast = NaN(p.numTrials, 1);
data.EstimatedLocation = NaN(p.numTrials, 1);
data.DifferenceLocation = NaN(p.numTrials, 1);
data.ResponseTime = NaN(p.numTrials, 1);

%%Trial Loop%%
for n = 1:p.numTrials
    
    %Define coordinates for drawing the stimuli every trial
    TargetXY = round([CenterX + p.eccentricity*cos(p.TrialEvents(n,4)*(pi/180))' CenterY - p.eccentricity*sin(p.TrialEvents(n,4)*(pi/180))']);
    
  
    %If loop for background grating
    if p.TrialEvents(n,1) == 2
        
        surroundtext = squeeze(background_grating(n,:,:))* (p.surroundContrast*p.Grey) + p.Grey;
        
        SurroundGrating = Screen('MakeTexture', window, surroundtext);
        Screen('DrawTexture', window, SurroundGrating, [], ...
            CenterRectOnPoint([0 0 p.Surround p.Surround], CenterX, CenterY), p.stimorientation);
    end
    
     %Show Grating at random location
    if p.TrialEvents(n,1) == 1
        Screen('FillOval', window, p.Grey, [CenterX-backgroundRadius CenterY-backgroundRadius CenterX+backgroundRadius CenterY+backgroundRadius]);
    end
    
    temp_grating = ones(p.PatchSize,p.PatchSize,2);
    temp_grating(:,:,1) = squeeze(patch_grating(n,:,:))*(p.TrialEvents(n,2) * p.Grey) + p.Grey;
    temp_grating(:,:,2) = centerGaussian*255;
    
    CenterStimulus = Screen('MakeTexture', window, temp_grating);
    Screen('DrawTexture', window, CenterStimulus, [], ...
        CenterRectOnPoint([0 0 p.PatchSize p.PatchSize], TargetXY(1), TargetXY(2)), p.stimorientation);
    
    Screen('FillOval', window, p.Black, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation]);
    Screen('FillOval', window, green, [CenterX-p.InnerFixation CenterY-p.InnerFixation CenterX+p.InnerFixation CenterY+p.InnerFixation]);
    Screen('Flip', window);
    WaitSecs(t.stimon);
    
    
        
    %Mask appears
    indx = Shuffle(1:round(t.flickertime/t?.flicker));
    StartMask = GetSecs;
    for a = indx
        if GetSecs > StartMask + t.flickertime - t.flicker
            break
        end        
        Screen('DrawTexture', window, Mask(a,:), [],CenterRectOnPoint([0 0 p.Surround p.Surround], CenterX, CenterY));
        Screen('FillOval', window, p.Black, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation]);
        Screen('FillOval', window, green, [CenterX-p.InnerFixation CenterY-p.InnerFixation CenterX+p.InnerFixation CenterY+p.InnerFixation]);
        Screen('Flip', window);
        WaitSecs(t.flicker);        
    end
    
    
    %Retention Interval
    Screen('FillOval', window, p.Grey, [CenterX-backgroundRadius CenterY-backgroundRadius CenterX+backgroundRadius CenterY+backgroundRadius]);
    Screen('FillOval', window, p.Black, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation]);
    Screen('FillOval', window, green, [CenterX-p.InnerFixation CenterY-p.InnerFixation CenterX+p.InnerFixation CenterY+p.InnerFixation]);
    Screen('Flip', window);
    WaitSecs(t.retention-t.flickertime);
    
    
    %Second Mask
    indx = Shuffle(1:round(t.flickertime/t.flicker));
    StartMask = GetSecs; 
    for a = indx        
        if GetSecs > StartMask + t.flickertime - t.flicker
            break
        end       
        Screen('DrawTexture', window, Mask(a,:), [],CenterRectOnPoint([0 0 p.Surround p.Surround], CenterX, CenterY));
        Screen('FillOval', window, p.Black, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation]);
        Screen('FillOval', window, green, [CenterX-p.InnerFixation CenterY-p.InnerFixation CenterX+p.InnerFixation CenterY+p.InnerFixation])
        Screen('Flip', window);
        WaitSecs(t.flicker);       
    end
    
    %Show probe at random location with random contrast
    Screen('FillOval', window, p.Grey, [CenterX-backgroundRadius CenterY-backgroundRadius CenterX+backgroundRadius CenterY+backgroundRadius]);
    
    initialAngle = p.TrialEvents(n,5);
    ProbeXY = round([CenterX + p.eccentricity*(cos(initialAngle*(pi/180)))' CenterY - p.eccentricity*(sin(initialAngle*(pi/180)))']);
    intial_contrast = p.TrialEvents(n,3); % random start contrast on each trial
    
    
    Probe = Screen('MakeTexture', window, squeeze(patch_target(n,:,:))* (intial_contrast *p.Grey) + p.Grey);
    Screen('DrawTexture', window, Probe, [], CenterRectOnPoint([0 0 p.PatchSize p.PatchSize], ProbeXY(1), ProbeXY(2)), p.stimorientation);
    Screen('FillOval', window, p.Black, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation]);
    Screen('FillOval', window, green, [CenterX-p.InnerFixation CenterY-p.InnerFixation CenterX+p.InnerFixation CenterY+p.InnerFixation]);
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
            Screen('FillOval', window, p.Grey, [CenterX-backgroundRadius CenterY-backgroundRadius CenterX+backgroundRadius CenterY+backgroundRadius]);
%             
            angles = ((startangle-angle)*3.8298);
            changeposition = (angles/(2*pi));
            initialAngle = initialAngle + changeposition; % update the location relative to last dial position
            ProbeXY = round([CenterX + p.eccentricity*(cos(initialAngle*(pi/180)))' CenterY - p.eccentricity*(sin(initialAngle*(pi/180)))']);
            
            Target = Screen('MakeTexture', window, squeeze(patch_target(n,:,:))* (intial_contrast*p.Grey) + p.Grey);
            Screen('DrawTexture', window, Target, [], CenterRectOnPoint([0 0 p.PatchSize p.PatchSize], ProbeXY(1), ProbeXY(2)), p.stimorientation)
            Screen('FillOval', window, p.Black, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation]);
            Screen('FillOval', window, green, [CenterX-p.InnerFixation CenterY-p.InnerFixation CenterX+p.InnerFixation CenterY+p.InnerFixation])
            Screen('Flip', window);
            
            startangle = angle;
        end
        if pmbutton == 1;
            locationTime = GetSecs;
%             % make sure angle stays in 0-360 range
            correctedAngle = mod(initialAngle, 360);
            
            data.EstimatedLocation(n) = correctedAngle;
%             
%             %make sure difference is in the 180 range
            difference = abs(p.TrialEvents(n,4) - data.EstimatedLocation(n));
%             
            if difference > 180
                difference = abs(difference - 360);
            end
            
            data.DifferenceLocation(n) = difference;
            
            data.ResponseTime_location(n) = (locationTime - starttrial);
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
%         % Query PowerMate button state and rotation angle in "clicks"
        [pmbutton_contrast, angle2] = PsychPowerMate('Get', powermate);
%         % 1st button is the "or" of the 1st mouse button and the actual PowerMate button
        if contrastangle ~= angle2
%             
%             % Convert turn of dial first to degrees and then to contrast:
%             %             angles = (angle * 3.8298)/360;
            angles = ((contrastangle-angle2)*3.8298);
            changecontrast = angles/360;
            intial_contrast = intial_contrast - changecontrast; % update the contrast relative to last dial position
            % Make sure we stay in range
            
            if intial_contrast > 1
                intial_contrast = 1;
            elseif intial_contrast < 0
                intial_contrast = 0.001;
            end
            Screen('FillOval', window, p.Grey, [CenterX-backgroundRadius CenterY-backgroundRadius CenterX+backgroundRadius CenterY+backgroundRadius]);
            
            Target = Screen('MakeTexture', window, squeeze(patch_target(n,:,:))* (intial_contrast*p.Grey) + p.Grey);
            Screen('DrawTexture', window, Target, [], CenterRectOnPoint([0 0 p.PatchSize p.PatchSize], ProbeXY(1), ProbeXY(2)), p.stimorientation)
            Screen('FillOval', window, p.Black, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation]);
            Screen('FillOval', window, green, [CenterX-p.InnerFixation CenterY-p.InnerFixation CenterX+p.InnerFixation CenterY+p.InnerFixation])
            Screen('Flip', window);
            
            contrastangle = angle2;
        end
        if pmbutton_contrast == 1;
            data.EstimatedContrast(n) = intial_contrast;
            data.DifferenceContrast(n) = p.TrialEvents(n,2) - data.EstimatedContrast(n);
            data.ResponseTime(n) = (GetSecs - starttrial);
            pmbutton_contrast = 0;
            break
        end
        
    end
    
GetClicks;

    % Check if esc button has been pressed
    [keyIsDown, keyCode] = PsychHID('KbQueueCheck', deviceNumber); %check response
    key = find(keyCode);
    if key == KbName('ESCAPE') % windows = 'esc', mac = 'ESCAPE' If user presses ESCAPE, exit the program.
        Screen('LoadNormalizedGammaTable', window, OriginalCLUT);
        Screen('CloseAll');
        ListenChar(1); % % Go back to unsuppressed mode
        FlushEvents('keyDown', deviceNumber);
        error('User exited program.');
    end
    
    
    %Center fixation
    Screen('FillOval', window, p.Grey, [CenterX-backgroundRadius CenterY-backgroundRadius CenterX+backgroundRadius CenterY+backgroundRadius]);
    Screen('FillOval', window, p.Black, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation]);
    Screen('FillOval', window, green, [CenterX-p.InnerFixation CenterY-p.InnerFixation CenterX+p.InnerFixation CenterY+p.InnerFixation]);
    Screen('Flip', window);
    WaitSecs(t.iti);
    
    %Break
    if n== ceil(p.numTrials/2)
        rest = GetSecs;
        
        Screen('FillOval', window, p.Grey, [CenterX-backgroundRadius CenterY-backgroundRadius CenterX+backgroundRadius CenterY+backgroundRadius]);

        RestText = ['Halfway there! You can take a short break now, ' '' '\n' ...
            'or press the dial to continue' '\n' '\n' ];
        DrawFormattedText(window, RestText, 'center', 'center', white);
        Screen('Flip', window);
        pmbuttonbreak = 0;
        
        while 1
            [pmbuttonbreak, a] = PsychPowerMate('Get', powermate);
            if pmbuttonbreak == 1;
                break;
            end
        end
        Screen('FillOval', window, p.Grey, [CenterX-backgroundRadius CenterY-backgroundRadius CenterX+backgroundRadius CenterY+backgroundRadius]);
        Screen('FillOval', window, p.Black, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation]);
        Screen('FillOval', window, green, [CenterX-p.InnerFixation CenterY-p.InnerFixation CenterX+p.InnerFixation CenterY+p.InnerFixation]);            Screen('Flip', window);
        WaitSecs(t.iti);
        t.Resttime = (GetSecs-rest)/60;
    end
end


t.EndTime = (GetSecs-StartTime)/60; %Get endtime of the experiment in seconds

%Draw some more text to the screen outside of the loop:
Screen(window,'TextSize',30);
ByebyeText = ['Great work! You have completed this run.' '\n' '\n' ...
    'Please let the experimenter know you have finished.'];
Screen('FillOval', window, p.Grey, [CenterX-backgroundRadius CenterY-backgroundRadius CenterX+backgroundRadius CenterY+backgroundRadius]);
DrawFormattedText(window, ByebyeText, 'center', 'center', [255 255 255]);
Screen('Flip', window);
WaitSecs(3);

Screen('LoadNormalizedGammaTable', window, OriginalCLUT);
Screen('CloseAll')
%
%% SAVE OUT THE DATA FILE
cd(datadir);
TheData(runnumber).t = t;
TheData(runnumber).p = p;
TheData(runnumber).data = data;
eval(['save surrSuppression_spatial_contrast_est_', p.Subject, '.mat TheData'])

cd(expdir);
