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

p.centerContrast = 0.75;
p.surroundContrast = 1;

p.centerSize = round(2*p.pixPerDeg);
p.surroundSize = p.screenWidthPixels(:,3);
p.gapSize = round(0.01*p.pixPerDeg); %space between center and surround annulus

p.backgroundRadius = 10 * p.pixPerDeg; %radius of the background circle
p.eccentricity = round(6 * p.pixPerDeg); %distance from center of screen to center of target
p.innerRadius = p.eccentricity - (p.centerSize/2+p.gapSize); %inner edge of surround annulus
p.outerRadius = p.eccentricity + (p.centerSize/2+p.gapSize); %outer edge of surround annulus

% Fixation dot size
p.innerFixation = round(0.075*p.pixPerDeg);
p.outerFixation = round(0.75*p.pixPerDeg);

% Grating frequency, orientation, phase
p.orientation = 90;
p.stimorientation = 90;
p.freq = 2.75;
p.frequency_center = p.centerSize/p.pixPerDeg * p.freq;
p.frequency_surround = p.surroundSize/p.pixPerDeg * p.freq;
p.numPhases = 1;
p.centerPhase = 360;
p.surroundPhase = p.centerPhase;


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
bgAnnulus(:,p.surroundSize/2:end) = 0;

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
foo = (eccen(:,p.surroundSize/2:end) <= p.backgroundRadius) -1;
tmpSurround(:,p.surroundSize/2:end) = foo;
surroundGrating = tmpSurround;


[window,rect] = Screen('OpenWindow', useScreen, colors.black, []); %test screen size [0 0 700 500]
    
OriginalCLUT = Screen('ReadNormalizedGammaTable', window);

% Enable alpha blending
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Define coordinates where to draw the fixation
CenterX = rect(3)/2; CenterY = rect(4)/2;
% Coordinates for location on left and right side of fixation
center =  [CenterX CenterY];
targetXY = round([CenterX+p.eccentricity*cos(90*(pi/180))' CenterY-p.eccentricity*sin(90*(pi/180))']);

Screen('TextStyle', window, 1);
Screen('TextSize', window, 16);
t.ifi = Screen('GetFlipInterval',window); % grab screen refresh rate


surroundTexture = (surroundGrating*p.surroundContrast)*colors.grey + colors.grey;
surroundTexture(:,:,2) = surroundTransparencyMask;
surroundStimulus = Screen('MakeTexture', window, surroundTexture);
Screen('DrawTexture', window, surroundStimulus, [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], CenterX, CenterY));

% Draw center grating
centerTexture = (centerGrating*p.centerContrast)*colors.grey + colors.grey;
centerTexture(:,:,2) = centerTransparencyMask;
centerStimulus = Screen('MakeTexture', window, centerTexture);
Screen('DrawTexture', window, centerStimulus, [], ...
        CenterRectOnPoint([0 0 p.centerSize p.centerSize], CenterX-p.eccentricity, CenterY));  
Screen('DrawTexture', window, centerStimulus, [], ...
        CenterRectOnPoint([0 0 p.centerSize p.centerSize], CenterX+p.eccentricity, CenterY));  

    % Draw fixation
Screen('FillOval', window, colors.black, [CenterX-p.outerFixation CenterY-p.outerFixation CenterX+p.outerFixation CenterY+p.outerFixation])
Screen('FillOval', window,colors.green,[CenterX-p.innerFixation CenterY-p.innerFixation CenterX+p.innerFixation CenterY+p.innerFixation]);
Screen('Flip',window);

GetClicks;
Screen('CloseAll')


