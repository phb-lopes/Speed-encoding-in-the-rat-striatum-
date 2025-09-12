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
    auxVelocity = (dist*factor)/dtPosition;     % Convertendo pra cm/s
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
    velocityAndBinPosition_WholeSession(i,3) = auxInd;
end
velocityAndBinPosition_beforeSleep = velocityAndBinPosition_WholeSession;

%% Removing sleep time from Session
clc, clear taskInfo 
velocityAndBinPosition_beforeSleep = [];
velocityAndBinPosition_afterSleep = [];

for k=2:length(allTimeCueOn+1)
    if abs(allTimeCueOn(k-1) - allTimeCueOn(k)) > 150 % This number (150) is 150 seconds and is used here to eliminate sleep period during the trials session
        taskInfo(k-1,:) = [allTaskIndex(k-1), allTimeCueOn(k-1), allTimeEndTrial(k-1), allRewardByTrial(k-1), k-1];
        lastTrialBeforeSleep = k-1;
    else
        taskInfo(k-1,:) = [allTaskIndex(k-1), allTimeCueOn(k-1), allTimeCueOn(k), allRewardByTrial(k-1), k-1];
    end
    if k == length(allTimeCueOn)
        taskInfo(k,:) = [allTaskIndex(k), allTimeCueOn(k), allTimeEndTrial(k), allRewardByTrial(k-1),k];
    end
end

% Get velocity only in Trials
for i = 1:lastTrialBeforeSleep
    I = find(velocityAndBinPosition_WholeSession(:,1) > taskInfo(i,2) & velocityAndBinPosition_WholeSession(:,1) < taskInfo(i,3));
    velocityAndBinPosition_beforeSleep = [velocityAndBinPosition_beforeSleep; velocityAndBinPosition_WholeSession(I,:)];
end

% Get velocity only in Trials
for i = lastTrialBeforeSleep+1:length(taskInfo)
    I = find(velocityAndBinPosition_WholeSession(:,1) > taskInfo(i,2) & velocityAndBinPosition_WholeSession(:,1) < taskInfo(i,3));
    velocityAndBinPosition_afterSleep = [velocityAndBinPosition_afterSleep; velocityAndBinPosition_WholeSession(I,:)];
end

%% BEFORE SLEEP
clc
acertos_no_intervalo_BF = NaN(size(velocityAndBinPosition_beforeSleep(:,1), 1), 1);
tipo_de_trial_BF = NaN(size(velocityAndBinPosition_beforeSleep(:,1), 1), 1);
lado_do_labirinto_BF = NaN(size(velocityAndBinPosition_beforeSleep(:,1), 1), 1);
num_trial_BF = NaN(size(velocityAndBinPosition_beforeSleep(:,1), 1), 1);

acertos_no_intervalo_AF = NaN(size(velocityAndBinPosition_afterSleep(:,1), 1), 1);
tipo_de_trial_AF = NaN(size(velocityAndBinPosition_afterSleep(:,1), 1), 1);
lado_do_labirinto_AF = NaN(size(velocityAndBinPosition_afterSleep(:,1), 1), 1);
num_trial_AF = NaN(size(velocityAndBinPosition_afterSleep(:,1), 1), 1);


for k=2:length(allTimeCueOn)
    intervalo_inicio = allTimeCueOn(k-1);
    intervalo_fim = allTimeCueOn(k);
    I = velocityAndBinPosition_beforeSleep(:,1) >= intervalo_inicio & velocityAndBinPosition_beforeSleep(:,1) <= intervalo_fim;
    acertos_no_intervalo_BF(I) = allRewardByTrial(k-1);
    tipo_de_trial_BF(I) = allTaskIndex(k-1);
    lado_do_labirinto_BF(I) = allRatChoice(k-1);
    num_trial_BF(I) = k-1;

    I = velocityAndBinPosition_afterSleep(:,1) >= intervalo_inicio & velocityAndBinPosition_afterSleep(:,1) <= intervalo_fim;
    acertos_no_intervalo_AF(I) = allRewardByTrial(k-1);
    tipo_de_trial_AF(I) = allTaskIndex(k-1);
    lado_do_labirinto_AF(I) = allRatChoice(k-1);
    num_trial_AF(I) = k-1;

    if k == length(allTimeCueOn)
        intervalo_inicio = allTimeCueOn(k);
        intervalo_fim = allTimeEndTrial(k);
        I = velocityAndBinPosition_beforeSleep(:,1) >= intervalo_inicio & velocityAndBinPosition_beforeSleep(:,1) <= intervalo_fim;
        acertos_no_intervalo_BF(I) = allRewardByTrial(k);
        tipo_de_trial_BF(I) = allTaskIndex(k-1);
        lado_do_labirinto_BF(I) = allRatChoice(k);
        num_trial_BF(I) = k;
        
        I = velocityAndBinPosition_afterSleep(:,1) >= intervalo_inicio & velocityAndBinPosition_afterSleep(:,1) <= intervalo_fim;
        acertos_no_intervalo_AF(I) = allRewardByTrial(k);
        tipo_de_trial_AF(I) = allTaskIndex(k-1);
        lado_do_labirinto_AF(I) = allRatChoice(k);
        num_trial_AF(I) = k;       
    end
end

velocityAndBinPosition_beforeSleep = [velocityAndBinPosition_beforeSleep, acertos_no_intervalo_BF, tipo_de_trial_BF, lado_do_labirinto_BF, num_trial_BF];
velocityAndBinPosition_afterSleep = [velocityAndBinPosition_afterSleep, acertos_no_intervalo_AF, tipo_de_trial_AF, lado_do_labirinto_AF, num_trial_AF];

timeInTrial = (velocityAndBinPosition_beforeSleep(end,1) - velocityAndBinPosition_beforeSleep(1,1)) + (velocityAndBinPosition_afterSleep(end,1) - velocityAndBinPosition_afterSleep(1,1));

%%
load('./IntPrincCellList.mat')

typeNeuronStr = arrayfun(@(n, d) sprintf('R%03d %d', n, d), TypeNeuron(:,1), TypeNeuron(:,2), 'UniformOutput', false);
TypeNeuronCell = num2cell(TypeNeuron);
TypeNeuronCell(:,1) = typeNeuronStr;
TypeNeuronNew = [TypeNeuronCell(:,1), TypeNeuronCell(:,3:5)];
%%
I = strcmp(TypeNeuronNew(:,1), nameRat);
neuronInfo = TypeNeuronNew(I,:);
%%
spikes = GetSpikeTimes('output','full');
spikes_OnlyInTrial = [];
for i = 1:lastTrialBeforeSleep
    I = find(spikes(:,1) > taskInfo(i,2) & spikes(:,1) < taskInfo(i,3));
    spikes_OnlyInTrial = [spikes_OnlyInTrial; spikes(I,:)];
end
spikes_beforeSleep = spikes_OnlyInTrial;

spikes_OnlyInTrial = [];
for i = lastTrialBeforeSleep+1:length(taskInfo)
    I = find(spikes(:,1) > taskInfo(i,2) & spikes(:,1) < taskInfo(i,3));
    spikes_OnlyInTrial = [spikes_OnlyInTrial; spikes(I,:)];
end
spikes_afterSleep = spikes_OnlyInTrial;

col2 = cell2mat(neuronInfo(:, 2));
col3 = cell2mat(neuronInfo(:, 3));
col4 = cell2mat(neuronInfo(:, 4));
tetrodoAndNeuro = [col2, col3, col4];
allNeuronSpikesPerdBinInfo_beforeSleep = {};

tic
count = 1;
neuronIndex = {};
for i = 1:size(tetrodoAndNeuro, 1)
    
    I = find(spikes_beforeSleep(:, 2) == tetrodoAndNeuro(i, 1) & spikes_beforeSleep(:, 3) == tetrodoAndNeuro(i, 2));

    auxAllSpk = spikes_beforeSleep(I, 1);  
    if length(auxAllSpk)/timeInTrial > .5
        spikePerSpatialBin = [];
        for j=1:length(auxAllSpk)
            [m,idx] = min(abs(auxAllSpk(j)-velocityAndBinPosition_beforeSleep(:,1)));
            spikePerSpatialBin(j,1) = auxAllSpk(j);                                 % Spike time
            spikePerSpatialBin(j,2) = velocityAndBinPosition_beforeSleep(idx,2);    % Velocity
            spikePerSpatialBin(j,3) = velocityAndBinPosition_beforeSleep(idx,3);    % Spatial Bin
            spikePerSpatialBin(j,4) = velocityAndBinPosition_beforeSleep(idx,4);    % Correct Choice (0 - Wrong | 1 - Right)
            spikePerSpatialBin(j,5) = velocityAndBinPosition_beforeSleep(idx,5);    % Cue (0 - Visual | 1 - Spatial)
            spikePerSpatialBin(j,6) = velocityAndBinPosition_beforeSleep(idx,6);    % Arm Choice (0 - Left | 1 - Right)
            spikePerSpatialBin(j,7) = velocityAndBinPosition_beforeSleep(idx,7);    % Trial number
        end
        allNeuronSpikesPerdBinInfo_beforeSleep{count} = spikePerSpatialBin;
        neuronIndex{1, count} = [tetrodoAndNeuro(i,1), tetrodoAndNeuro(i,2), tetrodoAndNeuro(i,3), length(auxAllSpk)/timeInTrial];
    else
        continue
    end
    I = find(spikes_afterSleep(:, 2) == tetrodoAndNeuro(i, 1) & spikes_afterSleep(:, 3) == tetrodoAndNeuro(i, 2));
    auxAllSpk = spikes_afterSleep(I, 1);  
    if length(auxAllSpk)/timeInTrial > .5
        spikePerSpatialBin = [];
        for j=1:length(auxAllSpk)
            [m,idx] = min(abs(auxAllSpk(j)-velocityAndBinPosition_afterSleep(:,1)));
            spikePerSpatialBin(j,1) = auxAllSpk(j);                                 % Spike time
            spikePerSpatialBin(j,2) = velocityAndBinPosition_afterSleep(idx,2);    % Velocity
            spikePerSpatialBin(j,3) = velocityAndBinPosition_afterSleep(idx,3);    % Spatial Bin
            spikePerSpatialBin(j,4) = velocityAndBinPosition_afterSleep(idx,4);    % Correct Choice (0 - Wrong | 1 - Right)
            spikePerSpatialBin(j,5) = velocityAndBinPosition_afterSleep(idx,5);    % Cue (0 - Visual | 1 - Spatial)
            spikePerSpatialBin(j,6) = velocityAndBinPosition_afterSleep(idx,6);    % Arm Choice (0 - Left | 1 - Right)
            spikePerSpatialBin(j,7) = velocityAndBinPosition_afterSleep(idx,7);    % Trial number
        end
        allNeuronSpikesPerdBinInfo_afterSleep{count} = spikePerSpatialBin;
        neuronIndex{2, count} = [tetrodoAndNeuro(i,1), tetrodoAndNeuro(i,2), tetrodoAndNeuro(i,3), length(auxAllSpk)/timeInTrial];
        count = count+1;
    else
        continue
    end
end
toc
%%
interval_time = 1;
overlap = 0.9;
step = interval_time - overlap;
if overlap == 0
    step = interval_time;
end

time_firing_rate_beforeSleep = min(velocityAndBinPosition_beforeSleep(:,1)):step:(max(velocityAndBinPosition_beforeSleep(:,1)) - interval_time);
time_firing_rate_afterSleep  = min(velocityAndBinPosition_afterSleep(:,1)):step:(max(velocityAndBinPosition_afterSleep(:,1)) - interval_time);

allNeuronsFiringRate_beforeSleep = {};

tic
for neuron=1:length(allNeuronSpikesPerdBinInfo_beforeSleep)
    
    spike_times = allNeuronSpikesPerdBinInfo_beforeSleep{neuron};
    auxFiringRate = [];
    
    for i = 1:length(time_firing_rate_beforeSleep)
        begin_interval = time_firing_rate_beforeSleep(i);
        end_interval = begin_interval + interval_time;

        spikes_no_intervalo = spike_times(:,1) >= begin_interval & spike_times(:,1) < end_interval;
        if any(spikes_no_intervalo)
            auxFiringRate(i,1) = sum(spikes_no_intervalo)/(end_interval-begin_interval);     % Spike FR/s
            auxFiringRate(i,2) = processVector(spike_times(spikes_no_intervalo,3));             % Spatial bin
            auxFiringRate(i,3) = processVector(spike_times(spikes_no_intervalo,4));             % Correct Choice (0 - Wrong | 1 - Right)
            auxFiringRate(i,4) = processVector(spike_times(spikes_no_intervalo,5));             % Cue (0 - Visual | 1 - Spatial)
            auxFiringRate(i,5) = processVector(spike_times(spikes_no_intervalo,6));             % Arm Choice (0 - Left | 1 - Right)
            auxFiringRate(i,6) = processVector(spike_times(spikes_no_intervalo,7));             % Trial number
        else            
            I = (velocityAndBinPosition_beforeSleep(:,1) >= begin_interval) & (velocityAndBinPosition_beforeSleep(:,1) <= end_interval);
            auxFiringRate(i,1) = sum(spikes_no_intervalo)/(end_interval-begin_interval);     % Spike FR/s
            auxFiringRate(i,2) = processVector(velocityAndBinPosition_beforeSleep(I,3));        % Spatial bin
            auxFiringRate(i,3) = processVector(velocityAndBinPosition_beforeSleep(I,4));        % Correct Choice (0 - Wrong | 1 - Right)
            auxFiringRate(i,4) = processVector(velocityAndBinPosition_beforeSleep(I,5));        % Cue (0 - Visual | 1 - Spatial)
            auxFiringRate(i,5) = processVector(velocityAndBinPosition_beforeSleep(I,6));        % Arm Choice (0 - Left | 1 - Right)
            auxFiringRate(i,6) = processVector(velocityAndBinPosition_beforeSleep(I,7));        % Trial number
        end
    end
    allNeuronsFiringRate_beforeSleep{neuron} = auxFiringRate;
    
    spike_times = allNeuronSpikesPerdBinInfo_afterSleep{neuron};
    auxFiringRate = [];
    for i = 1:length(time_firing_rate_afterSleep)
        begin_interval = time_firing_rate_afterSleep(i);
        end_interval = begin_interval + interval_time;

        spikes_no_intervalo = spike_times(:,1) >= begin_interval & spike_times(:,1) < end_interval;
        if any(spikes_no_intervalo)
            auxFiringRate(i,1) = sum(spikes_no_intervalo)/(end_interval-begin_interval);     % Spike FR/s
            auxFiringRate(i,2) = processVector(spike_times(spikes_no_intervalo,3));             % Spatial bin
            auxFiringRate(i,3) = processVector(spike_times(spikes_no_intervalo,4));             % Correct Choice (0 - Wrong | 1 - Right)
            auxFiringRate(i,4) = processVector(spike_times(spikes_no_intervalo,5));             % Cue (0 - Visual | 1 - Spatial)
            auxFiringRate(i,5) = processVector(spike_times(spikes_no_intervalo,6));             % Arm Choice (0 - Left | 1 - Right)
            auxFiringRate(i,6) = processVector(spike_times(spikes_no_intervalo,7));             % Trial number
        else 
            I = (velocityAndBinPosition_afterSleep(:,1) >= begin_interval) & (velocityAndBinPosition_afterSleep(:,1) <= end_interval);
            auxFiringRate(i,1) = sum(spikes_no_intervalo)/(end_interval-begin_interval);     % Spike FR/s
            auxFiringRate(i,2) = processVector(velocityAndBinPosition_afterSleep(I,3));         % Spatial bin
            auxFiringRate(i,3) = processVector(velocityAndBinPosition_afterSleep(I,4));         % Correct Choice (0 - Wrong | 1 - Right)
            auxFiringRate(i,4) = processVector(velocityAndBinPosition_afterSleep(I,5));        % Cue (0 - Visual | 1 - Spatial)
            auxFiringRate(i,5) = processVector(velocityAndBinPosition_afterSleep(I,6));         % Arm Choice (0 - Left | 1 - Right)
            auxFiringRate(i,6) = processVector(velocityAndBinPosition_afterSleep(I,7));         % Trial number
        end
    end
    allNeuronsFiringRate_afterSleep{neuron} = auxFiringRate;
    
    neuron 
end


data_velocity = velocityAndBinPosition_beforeSleep;
for i = 1:length(time_firing_rate_beforeSleep)
    begin_interval = time_firing_rate_beforeSleep(i);
    end_interval = begin_interval + interval_time;

    index_interval_beforeSleep = data_velocity(:, 1) >= begin_interval & data_velocity(:, 1) < end_interval;

    if any(index_interval_beforeSleep)
        velocity_bins_beforeSleep(i,1) = mean(data_velocity(index_interval_beforeSleep, 2));
    else
        velocity_bins_beforeSleep(i,1) = NaN; 
    end
end


data_velocity = velocityAndBinPosition_afterSleep;
for i = 1:length(time_firing_rate_afterSleep)
    begin_interval = time_firing_rate_afterSleep(i);
    end_interval = begin_interval + interval_time;

    index_interval_afterSleep = data_velocity(:, 1) >= begin_interval & data_velocity(:, 1) < end_interval;

    if any(index_interval_afterSleep)
        velocity_bins_afterSleep(i,1) = mean(data_velocity(index_interval_afterSleep, 2));
    else
        velocity_bins_afterSleep(i,1) = NaN;
    end
end
toc
%% REMOVING NAN ACCELERATION INDEX
indexToAvoidNaN_beforeSleep = find(~isnan(velocity_bins_beforeSleep));
if ~isempty(indexToAvoidNaN_beforeSleep)
    indices_selecionadas_beforeSleep = indexToAvoidNaN_beforeSleep;
    velocity_bins_beforeSleep = velocity_bins_beforeSleep(indices_selecionadas_beforeSleep);
end

indexToAvoidNaN_afterSleep = find(~isnan(velocity_bins_afterSleep));
if ~isempty(indexToAvoidNaN_afterSleep)
    indices_selecionadas_afterSleep = indexToAvoidNaN_afterSleep;
    velocity_bins_afterSleep = velocity_bins_afterSleep(indices_selecionadas_afterSleep);
end

% REMOVING VELOCITY == 0
indexToAvoidNozero_beforeSleep = find(velocity_bins_beforeSleep ~= 0);
if ~isempty(indexToAvoidNozero_beforeSleep)
    indices_nao_zero_beforeSleep = indexToAvoidNozero_beforeSleep;
    velocity_bins_beforeSleep = velocity_bins_beforeSleep(indices_nao_zero_beforeSleep);
end

indexToAvoidNozero_afterSleep = find(velocity_bins_afterSleep ~= 0);
if ~isempty(indexToAvoidNozero_afterSleep)
    indices_nao_zero_afterSleep = indexToAvoidNozero_afterSleep;
    velocity_bins_afterSleep = velocity_bins_afterSleep(indices_nao_zero_afterSleep);
end

% REMOVING VELOCITY > 60
inidexBellow60_beforeSleep = find(velocity_bins_beforeSleep < 60);
if ~isempty(inidexBellow60_beforeSleep)
    indicesBellow60_beforeSleep = inidexBellow60_beforeSleep;
    velocity_bins_beforeSleep = velocity_bins_beforeSleep(indicesBellow60_beforeSleep);
end

inidexBellow60_afterSleep = find(velocity_bins_afterSleep < 60);
if ~isempty(inidexBellow60_afterSleep)
    indicesBellow60_afterSleep = inidexBellow60_afterSleep;
    velocity_bins_afterSleep = velocity_bins_afterSleep(indicesBellow60_afterSleep);
end

velocity_bins = vertcat(velocity_bins_beforeSleep,velocity_bins_afterSleep);
 
%%
clc, clf


brainRegion = 'striatum';
directory = 'onlyInTrials';
colorScatter = '#7575a3';
colorCirclScatter = '#7575a3';

lag = 3;
bin_step = step;  
surNumber = 10000;
maxLag_seconds = lag;
maxLag_bins = round(maxLag_seconds/bin_step); 

data = {};

for i = 1:length(allNeuronsFiringRate_beforeSleep)
    if ~isempty(indexToAvoidNaN_beforeSleep)
        firing_rate_beforeSleep = allNeuronsFiringRate_beforeSleep{i}(indices_selecionadas_beforeSleep,:);
    else
        firing_rate_beforeSleep = allNeuronsFiringRate_beforeSleep{i};
    end
        
    if ~isempty(indexToAvoidNaN_afterSleep)
        firing_rate_afterSleep = allNeuronsFiringRate_afterSleep{i}(indices_selecionadas_afterSleep,:);
    else
        firing_rate_afterSleep = allNeuronsFiringRate_afterSleep{i};
    end

    if ~isempty(indexToAvoidNozero_beforeSleep)
        firing_rate_beforeSleep = firing_rate_beforeSleep(indices_nao_zero_beforeSleep,:);
    else
        firing_rate_beforeSleep = allNeuronsFiringRate_beforeSleep{i};
    end
    
    if ~isempty(indexToAvoidNozero_afterSleep)
        firing_rate_afterSleep = firing_rate_afterSleep(indices_nao_zero_afterSleep,:);
    else
        firing_rate_afterSleep = allNeuronsFiringRate_afterSleep{i};
    end
    
     if ~isempty(inidexBellow60_beforeSleep)
        firing_rate_beforeSleep = firing_rate_beforeSleep(indicesBellow60_beforeSleep,:);
    else
        firing_rate_beforeSleep = allNeuronsFiringRate_beforeSleep{i};
    end
    
    if ~isempty(inidexBellow60_afterSleep)
        firing_rate_afterSleep = firing_rate_afterSleep(indicesBellow60_afterSleep,:);
    else
        firing_rate_afterSleep = allNeuronsFiringRate_afterSleep{i};
    end
    
    firing_rate = vertcat(firing_rate_beforeSleep,firing_rate_afterSleep);
    
    
    time_firing_rate = 1:1:length(firing_rate);
    time_firing_rate = (1:1:length(firing_rate))*step;

    ini = 350;
    fim = 550;
    fig1 = gcf;
    set(gcf,'color','w')
    set(fig1,'position',[0, 0, 1400, 1050]);
    pos1 = [0.13,0.489517819706499,0.775,0.435482180293501];
    

    dataFr = firing_rate(:,1);
    dataVel = velocity_bins;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(5,3,[1:3])
    sgtitle(strcat(nameRat, ', FR=',num2str(sum(firing_rate(:,1))/length(firing_rate(:,1))),strcat(' Tetrode=',num2str(cell2mat(neuronInfo(i,2))),' Neuron=',num2str(cell2mat(neuronInfo(i,3))))),'fontweight', 'bold', 'FontSize', 19);
    hold on
    plot(time_firing_rate, dataVel, 'color','#e6ac00','lineW',3);
    xlim([ini fim])
    ylabel({'Velocidade','Média'});
    set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 10);
    set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 10);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(5,3,[4:6])
    plot(time_firing_rate, dataFr, 'color','#8300b4','lineW',3);
    xlim([ini fim])
    ylabel({'Firing', 'Rate (Hz)'});
    set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 10);
    set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 10);
    box off

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(5,3,[7:9])
    hold on
    plot(time_firing_rate, zscore(dataFr),  'color','#8300b4','lineW',3);
    plot(time_firing_rate, zscore(dataVel), 'color','#e6ac00','lineW',3);

    xlim([ini fim])
    ylabel('Z-score');
    xlabel('Tempo (s)');
    set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 10);
    set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 10);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(5,3,[[10:11],[13:14]])
    scatter(dataVel, dataFr, 'MarkerEdgeColor',colorScatter,'MarkerFaceColor',colorCirclScatter)
    hold on
    [y_fit, R2_1] = getFit_polyfit(dataVel, dataFr', 1);
    [r, p] = get_pearson(dataVel, dataFr);
    
    text(50,65,strcat(['r = ' num2str(round(r,3))]), 'color','k','FontWeight', 'bold','FontSize', 16)
    plot(dataVel, y_fit, 'k-', 'LineWidth', 3);
    ylabel('Firing Rate (Hz)');
    xlabel('Velocidade Média');
    set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 10);
    set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 10);
    yy = ylim;
    yTicks = yticks;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(5,3,[[12],[15]])
    
    for surIdx=1:surNumber
        shiftAmount = randi([20,50],1)/100;
        circularShift = round(length(dataVel)*shiftAmount,0);
        [ccg, ~] = xcorr(zscore(circshift(dataFr',circularShift)), zscore(dataVel), maxLag_bins, 'normalized');
        allSurCCG(:,surIdx) = ccg;
    end
    hold on;
    
    [ccg, lags] = xcorr(zscore(dataFr), zscore(dataVel), maxLag_bins, 'coeff');
    
    [m, ] = max(ccg);
    Threshold = max(mean(allSurCCG, 2)+2*std(allSurCCG,1, 2));
    m > Threshold;
    valoresLag_ms = lags*bin_step;  
    maxPointLength = 2*maxLag_bins+1;
    halfPointLength = maxLag_bins+1;
    
    plot(valoresLag_ms, ccg,'k','linew',3)
    plot(valoresLag_ms, mean(allSurCCG, 2),'r-','linew',3)
    
    A = mean(allSurCCG, 2)+2*std(allSurCCG,1, 2);
    B = mean(allSurCCG, 2)-2*std(allSurCCG,1, 2);   
    
    fill([valoresLag_ms, fliplr(valoresLag_ms)],[A' , fliplr(B')],'r','facealpha',0.2,'EdgeColor','none', 'HandleVisibility', 'off')
    
    plot([0 0], [-.2 .6], 'color','#9F9E9E','lineS','--','linew',2, 'HandleVisibility', 'off')
    xlabel('Lag (ms)')
    ylabel('CCG (r)')
    set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
    set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
    
    xticks(-maxLag_seconds:.5:maxLag_seconds)
    xlim([-maxLag_seconds, maxLag_seconds])
    yticks([-.2:.1:.6])
    ylim([-.2 .6])
    legend('Real','Sur (± 2 SD)', 'FontWeight', 'bold', 'FontSize', 11)
    set(legend,'Position',[0.81 0.38 0.1 0.01])
    
    box off
    data = [data; {nameRat{1}, neuronIndex{1,i}(1), neuronIndex{1,i}(2), firing_rate, velocity_bins, neuronIndex{1,i}(3), sum(firing_rate(:,1))/length(firing_rate(:,1))}];

    % pause();
    clf
    i
end

disp('Done!')
disp(nameRat{1})
%%
if ~isempty(data)
    saveName = strcat(replace(cell2mat(nameRat),' ','_'),'_veloAndFR_',num2str(interval_time),'sInterval_',num2str(maxLag_seconds),'sLag',num2str(overlap),'msOverlap.mat');
    saveName = strcat('./',strcat(num2str(interval_time),'sInterval/'),num2str(maxLag_seconds),'sLag/',num2str(overlap),'_overlap/',saveName);
    save(saveName, 'data', '-v7.3');
    disp('saved')
end
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

function result = processVector(inputVec)
    if isempty(inputVec)
        result = NaN;
        return;
    end
    
    numNaNs = sum(isnan(inputVec));
    validVec = inputVec(~isnan(inputVec));
    numValid = numel(validVec);
    
    if numNaNs > numValid
        result = NaN;
    else
        result = mode(validVec);
    end
end

function [y_fit, R2] = getFit_polyfit(data_vel, data_fr, degree)
    
    if isrow(data_vel)
        data_vel = data_vel'; 
    end
    if isrow(data_fr)
        data_fr = data_fr'; 
    end
    
    p = polyfit(data_vel, data_fr, degree);
    y_fit = polyval(p, data_vel);
    
    SSR = sum((data_fr - y_fit).^2);
    SST = sum((data_fr - mean(data_fr)).^2);
    
    R2 = 1 - (SSR / SST);
end

function [r, p] = get_pearson(data_vel, data_fr)
    if length(data_vel) ~= length(data_fr)
        error('As duas variáveis de entrada devem ter o mesmo tamanho');
    end
    [r, p] = corr(data_vel(:), data_fr(:), 'Type', 'Pearson');
end