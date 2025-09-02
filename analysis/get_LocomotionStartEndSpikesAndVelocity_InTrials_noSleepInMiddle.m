clear, clc, clf
cd /Users/phb_lopes/Documents/MATLAB/xmltree-main/@xmltree/private/

mex -O xml_findstr.c

cd /Users/phb_lopes/Library/CloudStorage/OneDrive-Pessoal/projetoEletrofisiologiaMacbook/
SetCurrentSession


%%
clc
timeStimulus = num2cell(GetEvents);
stims = GetEvents('output','descriptions');
aux = cell2mat(timeStimulus(:,1))/60;
minutes = fix(aux);
seconds = round(((round(((aux-fix(aux))*100)))*60/100)/100,2);
timeInMinutes = num2cell(minutes + seconds);
logStimulus = [stims, timeStimulus, timeInMinutes];

%% Processing name rat session
stims = GetEvents('output','descriptions');
info = split(stims(1:1),'-');
nameRat = split(info{1});
nameRat = nameRat(end);
nameRat = [nameRat, info(2)];
nameRat = join(nameRat);   

%% Get general data
clc
[allChoices, allTimeCueOn, allTimeEndTrial] = getTrialsChoice(logStimulus);
allEvents = getTrialsEvents(logStimulus);
allRatChoice = getTrialsRatChoices(allEvents);
allCues = getTrialsCues(allEvents);
[allTaskIndex, allRewardByTrial] = getTrialsReward(nameRat, allEvents, allRatChoice, allCues);
[mapPointsRight, mapPointsLeft, mapCentralArm] = getPosition(nameRat);

%%
clc, clear velocityAndBinPosition
dtPosition =  0.0256000000000001;
    
mainArm = 100;          % in cm. 
dist_firstToLastBin = abs(mapCentralArm(1)-mapCentralArm(end));   % Distance between first and last spatial bins
factor = mainArm/dist_firstToLastBin; 

pLab = GetPositions;
xSmoth = smooth(pLab(:,2));
ySmoth = smooth(pLab(:,3));
auxInd = 1;
% Calculate the animal velocity and position
velocityAndBinPosition_WholeSession = (pLab(:,1));

for i=2:length(pLab(:,1))
    dist = sqrt(((xSmoth(i-1) - xSmoth(i))^2) + ((ySmoth(i-1) - ySmoth(i))^2));
    auxVelocity = (dist*factor)/dtPosition;     % Converteing to cm/s
    velocityAndBinPosition_WholeSession(i,2) =  auxVelocity;
    
    if auxVelocity < 0.1
        velocityAndBinPosition_WholeSession(i,2) = 0;
    end
    if auxVelocity > 80
        velocityAndBinPosition_WholeSession(i,2) =  50;
    end
    positionMause = sqrt((pLab(i-1,2) - mapPointsRight(1,1)) ^2 + (pLab(i-1,3) - mapPointsRight(1,2)) ^2);
    
    % Get the map bin position in each time moment during the task
    for j=1:length(mapPointsRight)
        auxPositionR = sqrt((pLab(i-1,2) - mapPointsRight(j,1)) ^2 + (pLab(i-1,3) - mapPointsRight(j,2)) ^2);
        auxPositionL = sqrt((pLab(i-1,2) - mapPointsLeft(j,1)) ^2 + (pLab(i-1,3) - mapPointsLeft(j,2)) ^2);
        if auxPositionR<=auxPositionL
            auxPosition = auxPositionR;
        else
            auxPosition = auxPositionL;
        end
        if auxPosition <= positionMause
            positionMause = auxPosition;
            auxInd = j;
        end
    end
%     if pLab(:,1) > allTimeCueOn(k)
    velocityAndBinPosition_WholeSession(i,3) = auxInd;
end
velocityAndBinPosition = velocityAndBinPosition_WholeSession;

%% Removing sleep time from Session
clc, clear taskInfo 
velocityAndBinPosition = [];

for k=2:length(allTimeCueOn+1)
    if abs(allTimeCueOn(k-1) - allTimeCueOn(k)) > 3000 % This number is in seconds and is used here to eliminate sleep period during the trials session
        taskInfo(k-1,:) = [allTaskIndex(k-1), allTimeCueOn(k-1), allTimeEndTrial(k-1), allRewardByTrial(k-1), k-1];
    else
        taskInfo(k-1,:) = [allTaskIndex(k-1), allTimeCueOn(k-1), allTimeCueOn(k), allRewardByTrial(k-1), k-1];
    end
    if k == length(allTimeCueOn)
        taskInfo(k,:) = [allTaskIndex(k), allTimeCueOn(k), allTimeEndTrial(k), allRewardByTrial(k-1),k];
    end
end

taskTotalTime = sum(taskInfo(:,3)-taskInfo(:,2));

% Get velocity only in Trials
for i = 1:length(allTimeCueOn+1)
    I = find(velocityAndBinPosition_WholeSession(:,1) > taskInfo(i,2) & velocityAndBinPosition_WholeSession(:,1) < taskInfo(i,3));
    velocityAndBinPosition = [velocityAndBinPosition; velocityAndBinPosition_WholeSession(I,:)];
end

%% 
clc
acertos_no_intervalo = NaN(size(velocityAndBinPosition(:,1), 1), 1);
tipo_de_trial = NaN(size(velocityAndBinPosition(:,1), 1), 1);
lado_do_labirinto = NaN(size(velocityAndBinPosition(:,1), 1), 1);
num_trial = NaN(size(velocityAndBinPosition(:,1), 1), 1);

for k=2:length(allTimeCueOn)
    intervalo_inicio = allTimeCueOn(k-1);
    intervalo_fim = allTimeCueOn(k);
    I = velocityAndBinPosition(:,1) >= intervalo_inicio & velocityAndBinPosition(:,1) <= intervalo_fim;
    acertos_no_intervalo(I) = allRewardByTrial(k-1);
    tipo_de_trial(I) = allTaskIndex(k-1);
    lado_do_labirinto(I) = allRatChoice(k-1);
    num_trial(I) = k-1;
    
    if k == length(allTimeCueOn)
        intervalo_inicio = allTimeCueOn(k);
        intervalo_fim = allTimeEndTrial(k);
        I = velocityAndBinPosition(:,1) >= intervalo_inicio & velocityAndBinPosition(:,1) <= intervalo_fim;
        acertos_no_intervalo(I) = allRewardByTrial(k);
        tipo_de_trial(I) = allTaskIndex(k);
        lado_do_labirinto(I) = allRatChoice(k);
        num_trial(I) = k;
    end
end

velocityAndBinPosition = [velocityAndBinPosition, acertos_no_intervalo, tipo_de_trial, lado_do_labirinto, num_trial];

%%
timeInfoSession = 1;
lag = 3;
overlap = 0.9;

folder = strcat('/Users/phb_lopes/Library/CloudStorage/OneDrive-Pessoal/projetoEletrofisiologiaMacbook/',num2str(timeInfoSession),'sInterval/', num2str(lag),'sLag/',num2str(overlap),'_overlap/');
dataFsiAndLocomotionCells_all = readData(folder);

cells = {};
for i=1:length(dataFsiAndLocomotionCells_all)
    cells = [cells; dataFsiAndLocomotionCells_all{1,i}.data];
end
strCells = cells(cell2mat(cells(:,6)) >=21,:);

load('./ccgLagFrSpeedCellsInfoPosAndNeg_1sInterval_0.9overlap_3sLagNoSCwith2.5Lag.mat')
locomotioncellIndex = ccgLagFrLocomotionCellsData{4};

I = strcmp(strCells(:,1), nameRat{1});
neuronInfo = strCells(I,[1,2,3,6]);
%%
spikes = GetSpikeTimes('output','full');
spikes_OnlyInTrial = [];
for i = 1:length(taskInfo)
    I = find(spikes(:,1) > taskInfo(i,2) & spikes(:,1) < taskInfo(i,3));
    spikes_OnlyInTrial = [spikes_OnlyInTrial; spikes(I,:)];
end
spikes = spikes_OnlyInTrial;
col2 = cell2mat(neuronInfo(:, 2));
col3 = cell2mat(neuronInfo(:, 3));
col4 = cell2mat(neuronInfo(:, 4));
tetrodoAndNeuro = [col2, col3, col4];
allNeuronSpikesPerdBinInfo = {};

%% VELOCITY ABOVE 5cm/s IN FIRST SECONDS AND LATER INCREASES

time = velocityAndBinPosition(:, 1);
velocity = velocityAndBinPosition(:, 2);
threshold_peak = 5; % cm/s
minimum_portion = 0.5; % Minimum time portion where velocity is above 5 cm/s

% Find the indices where the velocity is below the threshold
indicesAbaixoDoLimiar = find(velocity < threshold_peak);
windows = []; 
lastWindowEnd = -Inf; 

% Check each index below the threshold
for i = 1:length(indicesAbaixoDoLimiar)
    startIndex = indicesAbaixoDoLimiar(i); % Start index of the window
    windowStart = time(startIndex); % Start of the window
    windowEnd = windowStart + 2; % End of the window (2 seconds later)

    % Ensure the new window does not overlap with the last one
    if windowStart < lastWindowEnd
        continue; % Skip if there is an overlap
    end

    % Find all indices within the window
    j = startIndex;
    while j <= length(time) && time(j) <= windowEnd
        j = j + 1;
    end

    % Check if all velocities are below the threshold in the window
    if all(velocity(startIndex:j-1) < threshold_peak)
        % Check if the velocity stays above 5 cm/s for most of the next 3 seconds
        nextWindowStart = windowEnd; % Start of the next window
        nextWindowEnd = nextWindowStart + 6; % End of the next window
        
        % Find indices for the next window
        indicesProximaJanela = find(time >= nextWindowStart & time < nextWindowEnd);
        
        % Calculate the proportion of time where velocity is above 5 cm/s
        if ~isempty(indicesProximaJanela)
            tempoAcima = sum(velocity(indicesProximaJanela) > threshold_peak);
            proporcao = tempoAcima / length(indicesProximaJanela);
            
            % Store the 6-second window if the proportion is greater than the minimum
            if proporcao >= minimum_portion
                windows = [windows; windowStart, nextWindowEnd]; 
                lastWindowEnd = nextWindowEnd; 
            end
        end
    end
end

% Display the found windows
disp('3-second windows where all velocities were below 5 cm/s or were above 5 cm/s for less than 1 second, with no overlap:');
windowsStart = windows;
length(windowsStart)
clear windows
%% VELOCITY ABOVE 5cm/s IN THE FIRST SECONDS AND THEN DECREASES

time = velocityAndBinPosition(:, 1);
velocity = velocityAndBinPosition(:, 2);
threshold_peak = 5; % cm/s
limite_baixo = 2; % cm/s
minimum_portion = 0.5; % Minimum proportion of time when velocity should be below 2 cm/s

% Find the indices where the velocity is above the threshold
indicesAcimaDoLimiar = find(velocity > threshold_peak);
windows = []; % To store the windows
lastWindowEnd = -Inf; % To track the end of the last window

% Check each index above the threshold
for i = 1:length(indicesAcimaDoLimiar)
    startIndex = indicesAcimaDoLimiar(i); % Start index of the window
    windowStart = time(startIndex); % Start of the window
    windowEnd = windowStart + 2; % End of the window (3 seconds later)

    % Ensure the new window does not overlap with the last one
    if windowStart < lastWindowEnd
        continue; % Skip if there is an overlap
    end

    % Find all indices within the window
    j = startIndex;
    while j <= length(time) && time(j) <= windowEnd
        j = j + 1;
    end

    % Check if all velocities are above the threshold in the window
    if all(velocity(startIndex:j-1) > threshold_peak)
        % Check if the velocity decreases below 2 cm/s for most of the next 3 seconds
        nextWindowStart = windowEnd; % Start of the next window
        nextWindowEnd = nextWindowStart + 6; % End of the next window
        
        % Find indices for the next window
        indicesProximaJanela = find(time >= nextWindowStart & time < nextWindowEnd);
        
        % Calculate the proportion of time where the velocity is below 2 cm/s
        if ~isempty(indicesProximaJanela)
            tempoAbaixo = sum(velocity(indicesProximaJanela) < limite_baixo);
            proporcao = tempoAbaixo / length(indicesProximaJanela);
            
            % Store the 6-second window if the proportion is greater than the minimum
            if proporcao >= minimum_portion
                windows = [windows; windowStart, nextWindowEnd]; 
                lastWindowEnd = nextWindowEnd; 
            end
        end
    end
end
windowsEnd = windows;
length(windowsEnd)
clear windows

%%
clc, clf
subplot(2,1,1)
hold on;

I = find(time > windowsStart(1,1) & time < windowsStart(1,2));
tam = length(I); % 273 for 7s 
soma = zeros(tam, 1);
count = 0;
for i = 1:length(windowsStart)
    I = find(time > windowsStart(i,1) & time < windowsStart(i,2));
    if length(I) ~= tam
        count = count + 1;
        continue
    else
        soma = soma + velocity(I);
        plot(velocity(I))
    end
end
count
meanStartMovment = soma / length(windowsStart);

plot(meanStartMovment,'r','lineW',5)
ylim([0 50])
xlabel('Time (s)');
ylabel('Mean velocity (cm/s)');
plot([78 78], [0 50], 'k--','linew',2)


subplot(2,1,2)
hold on;

I = find(time > windowsEnd(1,1) & time < windowsEnd(1,2));
tam = length(I); % 273 for 7s 
soma = zeros(tam, 1);
count = 0;

for i = 1:length(windowsEnd)
    I = find(time > windowsEnd(i,1) & time < windowsEnd(i,2));
    if length(I) ~= tam
        count = count + 1;
        continue
    else
        soma = soma + velocity(I);
        plot(velocity(I))
    end
end
count
meeanEndMovment = soma / length(windowsEnd);

plot(meeanEndMovment,'r','lineW',5)

plot([78 78], [0 50], 'k--','linew',2)
ylim([0 50])
xlabel('Time (s)');
ylabel('Mean velocity (cm/s)');
%%
clc
spkSAux = {};
spkRAux = {};
spkEAux = {};
spks = {};
secondsAnalysed = 3;
for neuron=1:size(tetrodoAndNeuro, 1)
    I = find(spikes(:, 2) == tetrodoAndNeuro(neuron, 1) & spikes(:, 3) == tetrodoAndNeuro(neuron, 2));
    auxAllSpk = spikes(I, 1);
    for i=1:length(windowsStart)
        spkSAux = [spkSAux; auxAllSpk(auxAllSpk(:,1)>windowsStart(i,1) & auxAllSpk(:,1)<windowsStart(i,2))-windowsStart(i,1)];
    end
    for i=1:length(windowsEnd)
        spkEAux = [spkEAux; auxAllSpk(auxAllSpk(:,1)>windowsEnd(i,1) & auxAllSpk(:,1)<windowsEnd(i,2))-windowsEnd(i,1)];
    end
    spks = [spks; {cell2mat(nameRat), tetrodoAndNeuro(neuron, 1), tetrodoAndNeuro(neuron, 2), tetrodoAndNeuro(neuron, 3), spkSAux, spkEAux}];
    spkSAux = {};
    spkEAux = {};
end

%%
clc
if ~isempty(spks)
    saveName = strcat(cell2mat(nameRat),'_spks_OC.mat');
    saveName = strcat('./startEndSpksAndVelocity/8sAnalise/spks/',saveName);
    save(saveName, 'spks', '-v7.3');
    disp('saved')
end
if ~isempty(meanStartMovment)
    saveName = strcat(cell2mat(nameRat),'_meanVelocityStart_OC.mat');
    saveName = strcat('./startEndSpksAndVelocity/8sAnalise/meanVelocity/start/',saveName);
    save(saveName, 'meanStartMovment', '-v7.3');
    disp('saved')
end
if ~isempty(meeanEndMovment)
    saveName = strcat(cell2mat(nameRat),'_meanVelocityEnd_OC.mat');
    saveName = strcat('./startEndSpksAndVelocity/8sAnalise/meanVelocity/end/',saveName);
    save(saveName, 'meeanEndMovment', '-v7.3');
    disp('saved')
end
cell2mat(nameRat)
%% ---------- Helper functions ----------
function [data] = readData(directory)
    filePattern = fullfile(directory, '*.mat'); % Change to whatever pattern you need.
    theFiles = dir(filePattern);

    for LengthEachCell = 1 : length(theFiles)
        baseFileName = theFiles(LengthEachCell).name;
        fullFileName = fullfile(theFiles(LengthEachCell).folder, baseFileName);
        fprintf(1, 'Reading %s\n', fullFileName);
        data{LengthEachCell} = load(fullFileName);
    end
end


function [allChoices, timeCueOn, timeEndTrial] = getTrialsChoice(logStimulus)
    newStr = logStimulus;

    timeEndTrial = []; 
    timeCueOn = []; 

    for n=2:length(newStr) 
        if (strcmp(newStr{n,1}, '6'))         %% get time of Beginning of  (SoundCue)
            timeCueOn = [timeCueOn; newStr{n,2}]; %#ok
        end

        if (strcmp(newStr{n-1,1}, '-56'))      %% get time of End Trial (EndOfCue - 1 position)
            timeEndTrial = [timeEndTrial; newStr{n,2}]; %#ok
        end

        if cell2mat(newStr(n)) == "1"
            choiceStr{n} = 'choiceRight'; %#ok
        end

        if cell2mat(newStr(n)) == "0"
            choiceStr{n} = 'choiceLeft'; %#ok
        end

    end

    allChoices = string(choiceStr(~cellfun('isempty',choiceStr)));

end

% Processing Data
function [newStr] = getTrialsEvents(logStimulus)
    newStr = logStimulus(:,1);

    sleepBegEvents = regexp(logStimulus(:,1),'beginning of .*Sleep.*', 'match');
    sleepBegEvents = cellfun(@(x)str2double(x), sleepBegEvents, 'UniformOutput', false);

    sleepEndEvents = regexp(logStimulus(:,1),'end of .*Sleep.*', 'match');
    sleepEndEvents = cellfun(@(x)str2double(x), sleepEndEvents, 'UniformOutput', false);

    for n=1:length(newStr)
        if ~isempty(sleepBegEvents{n})
            newStr{n} = 'begSleep';
        end
        if ~isempty(sleepEndEvents{n})
            newStr{n} = 'endSleep';
        end
        if cell2mat(newStr(n)) == "100"
            newStr{n} = 'EndOfTrial';
        end
         if cell2mat(newStr(n)) == "-56"
            newStr{n} = 'EndOfCue';
        end

        if cell2mat(newStr(n)) == "11"
            newStr{n} = 'rewardRight';
        end

        if cell2mat(newStr(n)) == "10"
            newStr{n} = 'rewardLeft';
        end

        if cell2mat(newStr(n)) == "1"
            newStr{n} = 'choiceRight';
        end

        if cell2mat(newStr(n)) == "0"
            newStr{n} = 'choiceLeft';
        end

        if cell2mat(newStr(n)) == "6"
            newStr{n} = 'SoundCue';
        end

        if cell2mat(newStr(n)) == "51"
            newStr{n} = 'taskA';
        end

        if cell2mat(newStr(n)) == "53" || cell2mat(newStr(n)) == "54"
            newStr{n} = 'taskB';
        end
    end
end

% Processing Rat choice Events data
function [ratChoice] = getTrialsRatChoices(events)
    stringChoice = regexp(events,'choiceLeft|choiceRight', 'match');
    ratChoice = string(stringChoice(~cellfun('isempty',stringChoice)));
    ratChoice = strrep(ratChoice,'choiceLeft','0');
    ratChoice = strrep(ratChoice,'choiceRight','1');
    ratChoice = double(ratChoice);
end


% Processing Cue Events data
function [cues] = getTrialsCues(events)
    stringCues = regexp(events,'Cue left|Cue right', 'match');
    cues = string(stringCues(~cellfun('isempty',stringCues)));
    cues = strrep(cues,'Cue left','0');
    cues = strrep(cues,'Cue right','1');
    cues = double(cues);
end

% Processing Reward data
function [taskIndex, rewardByTrial] = getTrialsReward(nameRat, events, ratChoice, cues)
    
    cuesAux = 2*cues-1;
    choiceAux = 2*ratChoice-1;
    choiceTv = cuesAux.*choiceAux*-1;    
    cont=0;
    rewardByTrial = zeros();

    for n=1:length(events)
        if (ismember (events(n), 'SoundCue')) && ((ismember (events(n+3), 'rewardLeft')) || (ismember (events(n+3), 'rewardRight')))
            cont=cont+1;
            rewardByTrial(cont) = 1;
        end
        if (ismember (events(n), 'SoundCue')) && (ismember (events(n+3), 'EndOfCue'))
            cont=cont+1;
            rewardByTrial(cont) = 0;
        end
    end

    contTaskSequenceRight = 0;
    contTaskSession = 0;
    taskIndex = [];
    contAuxSession = 0;
    contHit = 0;

    contPercente = 0;
    contVDSession = 0;
    contSDSession = 0;

    auxTask = 1.2;

    contVDReward = 0;
    contSDReward = 0;
    auxReward = 0;
    
    
    for n=1:length(rewardByTrial)
        if strcmp(nameRat,'R205 20120904')
            if contTaskSequenceRight >= 8 || contTaskSession >= 12 && contPercente >= 0.78

                if auxTask == 1.2   % Check if it is in Visual Task
                    auxTask = 1.6;  % Change to Spatial Taks

                    contVDSession = contAuxSession + contVDSession;
                    contAuxSession = contSDSession; %#ok
                    contVDReward = contVDReward+auxReward;

                else                % Change to Visual Task
                    auxTask = 1.2;

                    contSDSession = contAuxSession + contSDSession;
                    contAuxSession = contVDSession; %#ok
                    contSDReward = contSDReward+auxReward;
                end

                contTaskSession = 0;
                contTaskSequenceRight=0;  
                contAuxSession = 0;
                contPercente = 0; %#ok
                contHit = 0;
                auxReward = 0;
            end
        else 
            if contTaskSequenceRight >= 8

                if auxTask == 1.2   % Check if it is in Visual Task
                    auxTask = 1.6;  % Change to Spatial Taks

                    contVDSession = contAuxSession + contVDSession;
                    contAuxSession = contSDSession; %#ok
                    contVDReward = contVDReward+auxReward;

                else                % Change to Visual Task
                    auxTask = 1.2;

                    contSDSession = contAuxSession + contSDSession;
                    contAuxSession = contVDSession; %#ok
                    contSDReward = contSDReward+auxReward;
                end

                contTaskSession = 0;
                contTaskSequenceRight=0;  
                contAuxSession = 0;
                contPercente = 0; %#ok
                contHit = 0;
                auxReward = 0;
            end
        end
        
        
        if rewardByTrial(n) == 1
            contTaskSequenceRight = contTaskSequenceRight+1;
            contHit = contHit+1;
            auxReward = auxReward+1;
        else
            contTaskSequenceRight = 0;
        end
        contTaskSession = contTaskSession+1;
        contAuxSession = contAuxSession+1;
        contPercente = contHit/contAuxSession;

        taskIndex(n)= auxTask; %#ok
    end

    for n=1:length(taskIndex)
        if taskIndex(n) == 1.6
            taskIndex(n) = 1; %#ok
        else
            taskIndex(n) = 0; %#ok
        end
    end
    
    
    %%%%%%%%%%% R209 - 20121113 %%%%%%%%%%%%
    if strcmp(nameRat,'R209 20121113')
        taskIndex(1:end) = 0;
        taskIndex(16:29) = 1;
        taskIndex(30:43) = 0;
        taskIndex(44:56) = 1;
        taskIndex(57:62) = 0;
        taskIndex(63:75) = 1;
    end
    
    %%%%%%%%%%%% R209 - 20121114 %%%%%%%%%%%%
    if strcmp(nameRat,'R209 20121114')
        taskIndex(1:end) = 0;
        taskIndex(17:26) = 1;
        taskIndex(27:39) = 0;
        taskIndex(40:46) = 1;
        taskIndex(47:54) = 0;
        taskIndex(55:59) = 1;
        taskIndex(60:82) = 0;
        taskIndex(83:111) = 1;
    end
end


function [mapPointsRightWay, mapPointsLeftWay, xMapCentralArm] = getPosition(nameRat)
    
    %%% [R205 - 20120902, 20120903 & 20120904] %%%
    r205 = {'R205 20120902', 'R205 20120903', 'R205 20120904'};
    if contains(nameRat, r205)
        xMapCentralArm = (0.2012:0.042433:0.7104);
        yMapCentralArm(1,1:13) = 0.4625;

        %%%%%%%%%%%%%%%%%% [Right Arm] %%%%%%%%%%%%%%%%%%%%%%%%%%%
        xMapRightArm = (0.6829:0.00535:0.7043);
        yMapRightArm = (0.5167:0.0375:0.6667);

        xMapRightLateralArm = (0.2774:0.0342:0.6546);
        yMapRightLateralArm = (0.5458:0.01113:0.6683);

        xMapRightLateralArm = xMapRightLateralArm(end:-1:1);
        yMapRightLateralArm = yMapRightLateralArm(end:-1:1);
        xMapRightArm = xMapRightArm(end:-1:1);


        %%%%%%%%%%%%%%%%%% [Left Arm] %%%%%%%%%%%%%%%%%%%%%%%%%%%
        xMapLeftArm = (0.689:0.004575:0.7104);
        yMapLeftArm = (0.2167:0.0479:0.4083);

        xMapLeftLateralArm = (0.2774:0.0342:0.6546);
        yMapLeftLateralArm = (0.2208:0.0125:0.3583);

        xMapLeftLateralArm = xMapLeftLateralArm(end:-1:1);
        xMapLeftArm = xMapLeftArm(end:-1:1);
        yMapLeftArm = yMapLeftArm(end:-1:1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        mapPointsLeftWay(:,1) = [xMapCentralArm, xMapLeftArm, xMapLeftLateralArm];
        mapPointsLeftWay(:,2) = [yMapCentralArm, yMapLeftArm, yMapLeftLateralArm];

        mapPointsRightWay(:,1) = [xMapCentralArm, xMapRightArm, xMapRightLateralArm];
        mapPointsRightWay(:,2) = [yMapCentralArm, yMapRightArm, yMapRightLateralArm];
    end
    
    %%% [R209 - 20121113, 20121114] %%% 
    r209_13_14 = {'R209 20121113', 'R209 20121114'};
    if contains(nameRat, r209_13_14)

        xMapCentralArm = (0.1646:0.045225:0.7073);
        yMapCentralArm = (0.2875:0.00041:0.2925);
        yMapCentralArm = yMapCentralArm(end:-1:1);

        %%%%%%%%%%%%%%%%%% [Right Arm] %%%%%%%%%%%%%%%%%%%%%%%%%%%
        xMapRightArm = (0.689:0.003825:0.7043);
        yMapRightArm = (0.3292:0.0427:0.5);

        xMapRightLateralArm = (0.2195:0.0382:0.6402);
        yMapRightLateralArm = (0.5458:0.01113:0.6683)-.185;

        xMapRightLateralArm = xMapRightLateralArm(end:-1:1);
        yMapRightLateralArm = yMapRightLateralArm(end:-1:1);
        xMapRightArm = xMapRightArm(end:-1:1);

        %%%%%%%%%%%%%%%%%% [Left Arm] %%%%%%%%%%%%%%%%%%%%%%%%%%%
        xMapLeftArm = (0.689:0.004575:0.7073);
        yMapLeftArm = (0.0467:0.0518:0.2542);

        xMapLeftLateralArm = (0.2195:0.0382:0.6402);
        yMapLeftLateralArm = (0.05417:0.0147:0.2167);

        xMapLeftLateralArm = xMapLeftLateralArm(end:-1:1);
        xMapLeftArm = xMapLeftArm(end:-1:1);
        yMapLeftArm = yMapLeftArm(end:-1:1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        mapPointsLeftWay(:,1) = [xMapCentralArm, xMapLeftArm, xMapLeftLateralArm];
        mapPointsLeftWay(:,2) = [yMapCentralArm, yMapLeftArm, yMapLeftLateralArm];

        mapPointsRightWay(:,1) = [xMapCentralArm, xMapRightArm, xMapRightLateralArm];
        mapPointsRightWay(:,2) = [yMapCentralArm, yMapRightArm, yMapRightLateralArm];

    end
    
    %%% [R209 - 20121115] %%% 
    r209_15 = {'R209 20121115'};
    if contains(nameRat, r209_15)
        xMapCentralArm = (0.1646: 0.04624:0.7195);
        yMapCentralArm = (0.4625:0.00041:0.4675);
        yMapCentralArm = yMapCentralArm(end:-1:1);

        %%%%%%%%%%%%%%%%%% [Right Arm] %%%%%%%%%%%%%%%%%%%%%%%%%%%
        xMapRightArm = (0.689:0.003825:0.7043)+.005;
        yMapRightArm = (0.3292:0.0427:0.5)+.18;

        xMapRightLateralArm = (0.2195:0.0382:0.6402);
        yMapRightLateralArm = (0.5458:0.01113:0.6683);

        xMapRightLateralArm = xMapRightLateralArm(end:-1:1);
        yMapRightLateralArm = yMapRightLateralArm(end:-1:1);
        xMapRightArm = xMapRightArm(end:-1:1);
        
        %%%%%%%%%%%%%%%%%% [Left Arm] %%%%%%%%%%%%%%%%%%%%%%%%%%%
        xMapLeftArm = (0.689:0.004575:0.7073)+.005;
        yMapLeftArm = (0.2267:0.04437:0.4042);

        xMapLeftLateralArm = (0.2195:0.0382:0.6402);
        yMapLeftLateralArm = (0.05417:0.0147:0.2167)+.18;

        xMapLeftLateralArm = xMapLeftLateralArm(end:-1:1);
        xMapLeftArm = xMapLeftArm(end:-1:1);
        yMapLeftArm = yMapLeftArm(end:-1:1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        mapPointsLeftWay(:,1) = [xMapCentralArm, xMapLeftArm, xMapLeftLateralArm];
        mapPointsLeftWay(:,2) = [yMapCentralArm, yMapLeftArm, yMapLeftLateralArm];

        mapPointsRightWay(:,1) = [xMapCentralArm, xMapRightArm, xMapRightLateralArm];
        mapPointsRightWay(:,2) = [yMapCentralArm, yMapRightArm, yMapRightLateralArm];
    end

    %%% [R268 - 20140930] %%%
    R287_30 = {'R268 20140930'};
    if contains(nameRat, R287_30)  
        xMapCentralArm = (0.1646: 0.04624:0.7195)-.01;
        yMapCentralArm = (0.4625:0.00041:0.4675)+.01;
        yMapCentralArm = yMapCentralArm(end:-1:1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%% [Right Arm] %%%%%%%%%%%%%%%%%%%%%%%%%%%
        xMapRightArm = (0.689:0.003825:0.7043)+.002;
        yMapRightArm = (0.3292:0.0427:0.5)+.19;

        xMapRightLateralArm = (0.2195:0.0382:0.6402);
        yMapRightLateralArm = (0.5458:0.01113:0.6683);

        xMapRightLateralArm = xMapRightLateralArm(end:-1:1);
        yMapRightLateralArm = yMapRightLateralArm(end:-1:1);
        xMapRightArm = xMapRightArm(end:-1:1);


        %%%%%%%%%%%%%%%%%%%%%%%%%%% [Left Arm] %%%%%%%%%%%%%%%%%%%%%%%%%%%
        xMapLeftArm = (0.689:0.004575:0.7073)+.002;
        yMapLeftArm = (0.2267:0.04437:0.4042)+.01;

        xMapLeftLateralArm = (0.2195:0.0382:0.6402);
        yMapLeftLateralArm = (0.05417:0.0147:0.2167)+.19;

        xMapLeftLateralArm = xMapLeftLateralArm(end:-1:1);
        xMapLeftArm = xMapLeftArm(end:-1:1);
        yMapLeftArm = yMapLeftArm(end:-1:1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        mapPointsLeftWay(:,1) = [xMapCentralArm, xMapLeftArm, xMapLeftLateralArm];
        mapPointsLeftWay(:,2) = [yMapCentralArm, yMapLeftArm, yMapLeftLateralArm];

        mapPointsRightWay(:,1) = [xMapCentralArm, xMapRightArm, xMapRightLateralArm];
        mapPointsRightWay(:,2) = [yMapCentralArm, yMapRightArm, yMapRightLateralArm];
    end
    
    %%% [R268 - 20141016] %%%
    R287_16 = {'R268 20141016'};
    if contains(nameRat, R287_16)  
        xMapCentralArm = (0.1646: 0.04624:0.7195)-.02;
        yMapCentralArm = (0.4625:0.00041:0.4675)-.17;
        yMapCentralArm = yMapCentralArm(end:-1:1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%% [Right Arm] %%%%%%%%%%%%%%%%%%%%%%%%%%%
        xMapRightArm = (0.6768:0.00687:0.7043)+.005;
        yMapRightArm = (0.3292:0.0427:0.5)+.01;

        xMapRightLateralArm = (0.2226:0.0376:0.6372);
        yMapRightLateralArm = (0.3625:0.01288:0.5042);

        xMapRightLateralArm = xMapRightLateralArm(end:-1:1);
        yMapRightLateralArm = yMapRightLateralArm(end:-1:1);
        xMapRightArm = xMapRightArm(end:-1:1);


        %%%%%%%%%%%%%%%%%%%%%%%%%%% [Left Arm] %%%%%%%%%%%%%%%%%%%%%%%%%%%
        xMapLeftArm = (0.689:0.004575:0.7073);
        yMapLeftArm = (0.2267:0.04437:0.4042)-.16;

        xMapLeftLateralArm = (0.2195:0.0382:0.6402);
        yMapLeftLateralArm = (0.07917:0.0140:0.2333);

        xMapLeftLateralArm = xMapLeftLateralArm(end:-1:1);
        xMapLeftArm = xMapLeftArm(end:-1:1);
        yMapLeftArm = yMapLeftArm(end:-1:1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        mapPointsLeftWay(:,1) = [xMapCentralArm, xMapLeftArm, xMapLeftLateralArm];
        mapPointsLeftWay(:,2) = [yMapCentralArm, yMapLeftArm, yMapLeftLateralArm];

        mapPointsRightWay(:,1) = [xMapCentralArm, xMapRightArm, xMapRightLateralArm];
        mapPointsRightWay(:,2) = [yMapCentralArm, yMapRightArm, yMapRightLateralArm];
    end
    
    %%% [R287 - 20150105, 20150107, 20150113, 20150114, 20150115] %%%
    R287 = {'R287 20150105', 'R287 20150107', 'R287 20150113', 'R287 20150114', 'R287 20150115'};
    if contains(nameRat, R287) 
        xMapCentralArm = (0.186: 0.0498:0.7845);
        yMapCentralArm = (0.4625:0.00041:0.4675)+.13;
        yMapCentralArm = yMapCentralArm(end:-1:1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%% [Right Arm] %%%%%%%%%%%%%%%%%%%%%%%%%%%
        xMapRightArm = (0.759:0.00612:0.7835);
        yMapRightArm = (0.6417:0.052:0.85);

        xMapRightLateralArm = (0.2317:0.0443:0.7195);
        yMapRightLateralArm = (0.6625:0.01688:0.8482);

        xMapRightLateralArm = xMapRightLateralArm(end:-1:1);
        yMapRightLateralArm = yMapRightLateralArm(end:-1:1);
        xMapRightArm = xMapRightArm(end:-1:1);


        %%%%%%%%%%%%%%%%%%%%%%%%%%% [Left Arm] %%%%%%%%%%%%%%%%%%%%%%%%%%%
        xMapLeftArm = (0.7622:0.005325:0.7835);
        yMapLeftArm = (0.2958:0.05835:0.5292);

        xMapLeftLateralArm = (0.2317:0.0443:0.7195);
        yMapLeftLateralArm = (0.3042:0.0185:0.5083);

        xMapLeftLateralArm = xMapLeftLateralArm(end:-1:1);
        xMapLeftArm = xMapLeftArm(end:-1:1);
        yMapLeftArm = yMapLeftArm(end:-1:1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        mapPointsLeftWay(:,1) = [xMapCentralArm, xMapLeftArm, xMapLeftLateralArm];
        mapPointsLeftWay(:,2) = [yMapCentralArm, yMapLeftArm, yMapLeftLateralArm];

        mapPointsRightWay(:,1) = [xMapCentralArm, xMapRightArm, xMapRightLateralArm];
        mapPointsRightWay(:,2) = [yMapCentralArm, yMapRightArm, yMapRightLateralArm];
    end


    %%%[R294 - 20150313, 20150317, 20150401] %%%
    R287 = {'R294 20150313', 'R294 20150317', 'R294 20150401'};
    if contains(nameRat, R287) 

        xMapCentralArm = (0.186:0.05055:0.7927);
        yMapCentralArm = (0.5792:0.00151:0.5974);
        yMapCentralArm = yMapCentralArm(end:-1:1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%% [Right Arm] %%%%%%%%%%%%%%%%%%%%%%%%%%%
        xMapRightArm = (0.7561:0.006849:0.7835);
        yMapRightArm = (0.625:0.059375:0.8625);

        xMapRightLateralArm = (0.2317:0.0443:0.7195);
        yMapRightLateralArm = (0.6833:0.01499:0.8482);

        xMapRightLateralArm = xMapRightLateralArm(end:-1:1);
        yMapRightLateralArm = yMapRightLateralArm(end:-1:1);
        xMapRightArm = xMapRightArm(end:-1:1);


        %%%%%%%%%%%%%%%%%%%%%%%%%%% [Left Arm] %%%%%%%%%%%%%%%%%%%%%%%%%%%
        xMapLeftArm = (0.7622:0.005325:0.7835);
        yMapLeftArm = (0.2958:0.05835:0.5292);

        xMapLeftLateralArm = (0.2317:0.0443:0.7195);
        yMapLeftLateralArm = (0.3042:0.0185:0.5083)+.02;

        xMapLeftLateralArm = xMapLeftLateralArm(end:-1:1);
        xMapLeftArm = xMapLeftArm(end:-1:1);
        yMapLeftArm = yMapLeftArm(end:-1:1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        mapPointsLeftWay(:,1) = [xMapCentralArm, xMapLeftArm, xMapLeftLateralArm];
        mapPointsLeftWay(:,2) = [yMapCentralArm, yMapLeftArm, yMapLeftLateralArm];

        mapPointsRightWay(:,1) = [xMapCentralArm, xMapRightArm, xMapRightLateralArm];
        mapPointsRightWay(:,2) = [yMapCentralArm, yMapRightArm, yMapRightLateralArm];
    end
    
    clf
    pLab = GetPositions;
    
    f = gcf;
    set(f,'position',[0, 0, 1400, 1050]);
    plot(pLab(:,2),pLab(:,3), 'k.','markersize', 16)


    hold on

    for i=1:length(xMapCentralArm)
        text(xMapCentralArm(i),yMapCentralArm(i)-.01,num2str(i),'fontweight', 'bold', 'FontSize', 18,'Color', 'r')
    end

    for i=1:length(xMapRightArm)
        text(xMapRightArm(i),yMapRightArm(i)-.01,num2str(i+length(xMapCentralArm)),'fontweight', 'bold', 'FontSize', 18,'Color', 'r')
    end

    for i=1:length(xMapLeftArm)
        text(xMapLeftArm(i),yMapLeftArm(i)-.01,num2str(i+length(xMapCentralArm)),'fontweight', 'bold', 'FontSize', 18,'Color', 'r')
    end

    for i=1:length(xMapLeftLateralArm)
        text(xMapLeftLateralArm(i),yMapLeftLateralArm(i)-.01,num2str(i+length(xMapCentralArm)+length(xMapLeftArm)),'fontweight', 'bold', 'FontSize', 18,'Color', 'r')
    end

    for i=1:length(xMapRightLateralArm)
        text(xMapRightLateralArm(i),yMapRightLateralArm(i)-.01,num2str(i+length(xMapCentralArm)+length(xMapRightArm)),'fontweight', 'bold', 'FontSize', 18,'Color', 'r')
    end

    plot(xMapCentralArm, yMapCentralArm,'o','MarkerEdgeColor','g','MarkerFaceColor','g','markersize', 7)
    plot(xMapRightArm, yMapRightArm,'o','MarkerEdgeColor','b','MarkerFaceColor','b','markersize', 7)
    plot(xMapLeftArm, yMapLeftArm,'o','MarkerEdgeColor','y','MarkerFaceColor','y','markersize', 7)
    plot(xMapLeftLateralArm, yMapLeftLateralArm,'o','MarkerEdgeColor','#FFA500','MarkerFaceColor','#FFA500','markersize', 7)
    plot(xMapRightLateralArm, yMapRightLateralArm,'o','MarkerEdgeColor','#8A2BE2','MarkerFaceColor','#8A2BE2','markersize', 7)

    %%%%%%%%%%%%%%%%%%%%%%%%%%% [R205 - 20120902, 20120903 & 20120904] %%%%%%%%%%%%%%%%%%%%%%%%%%%
    text(.7,.77,'Right','fontweight', 'bold', 'FontSize', 18)
    text(.7,.15,'Left','fontweight', 'bold', 'FontSize', 18)

    title(nameRat,'fontweight', 'bold', 'FontSize', 19)
    xticks([])
    yticks([])

end