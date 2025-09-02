%% Clear environment and set working directory
clear, clc, clf
cd /Users/phb_lopes/Library/CloudStorage/OneDrive-Pessoal/projetoEletrofisiologiaMacbook

%% Load firing rate data by trial (with choice info)
folder = strcat('./firingRateByTrial_withChoiceAdded/');
dataFRByTrial = readData(folder);

%% Load session-specific striatal data
clc
timeInfoSession = 1;
overlap = 0.9;
lag = 3;

folder = strcat('/Users/phb_lopes/Library/CloudStorage/OneDrive-Pessoal/projetoEletrofisiologiaMacbook/', ...
    num2str(timeInfoSession),'sInterval/', num2str(lag),'sLag/', num2str(overlap),'_overlap/');
dataFsiAndLocomotionCells_all = readData(folder);

% Collect all cells from all sessions into one array
allCells = {};
for i=1:length(dataFsiAndLocomotionCells_all)
    allCells = [allCells; dataFsiAndLocomotionCells_all{1,i}.data];
end

% Keep only striatal cells (cell type >= 21)
strCells = allCells(cell2mat(allCells(:,6)) >=21,:);

%% Load cross-correlation results (lag firing rate vs speed correlation)
load('./ccgLagFrSpeedCellsInfoPosAndNeg_1sInterval_0.9overlap_3sLagNoSCwith2.5Lag.mat')
strCells(:, end+1) = num2cell(ccgLagFrLocomotionCellsData{1,1}); % lag info
strCells(:, end+1) = num2cell(ccgLagFrLocomotionCellsData{1,4}); % speed correlation info

%% Combine firing rate by trial data
clc
allData = [];
for i = 1:length(dataFRByTrial)
     aux = dataFRByTrial{i};
    
    % Get field name of struct
    nome_campo = fieldnames(aux);
    
    % Access trial matrix
    celula = aux.(nome_campo{1});
    
    % Select identifying columns (1–3) and trial condition column (22)
    selecte_columns = celula(:,[1:3,22]);

    % Append to combined dataset
    allData = [allData; selecte_columns];
end

%% Merge trial condition info into striatal cell data
lastColoum = length(strCells(1,:));
for i = 1:size(strCells, 1)
    for j = 1:size(allData, 1)
        % Compare first 3 identifying columns (animal, session, cellID)
        if isequal(strCells{i,1}, allData{j,1}) && ...
           isequal(strCells{i,2}, allData{j,2}) && ...
           isequal(strCells{i,3}, allData{j,3})
       
            % Add trial condition column to strCells
            strCells{i,lastColoum+1} = allData{j,4};
            break; % Found match, stop searching
        end
    end
end

%% Define task condition markers
clc
visualMarker = 0;
spatialMarker = 1;
rewardMarker = 1;
noRewardMarker = 0;
leftChoice = 0;
rightChoice = 1;

surNumber = 10000; % number of surrogate permutations

% Containers for neurons modulated by different conditions
neuronModulate_v_s = [];   % visual vs spatial
neuronModulate_r_nr = [];  % reward vs no-reward
neuronModulate_lc_rc = []; % left vs right

%% Test modulation of individual neurons
for i = 1:length(strCells)
    data = cell2mat(strCells(i,end));
    
    % Separate firing rates by trial condition
    frAux_v  = data(data(:,3)==visualMarker);
    frAux_s  = data(data(:,3)==spatialMarker);
    frAux_r  = data(data(:,4)==rewardMarker);
    frAux_nr = data(data(:,4)==noRewardMarker);
    frAux_cl = data(data(:,5)==leftChoice);
    frAux_cr = data(data(:,5)==rightChoice);

    % Initialize surrogate deltas
    allSurDelta_v_s = [];
    allSurDelta_r_nr = [];
    allSurDelta_lc_rc = [];
    
    % Shuffle labels to create surrogate distributions
    for surIdx=1:surNumber
        surData = data;
        
        % Visual vs spatial surrogate
        I = randperm(length(data(:,3)));
        frAux_v_sur = surData(surData(I,3)==visualMarker);
        frAux_s_sur = surData(surData(I,3)==spatialMarker);
        allSurDelta_v_s(:,surIdx) = mean(frAux_v_sur)-mean(frAux_s_sur);
        
        % Reward vs no-reward surrogate
        I = randperm(length(data(:,4)));
        frAux_r_sur = surData(surData(I,4)==rewardMarker);
        frAux_nr_sur = surData(surData(I,4)==noRewardMarker);
        allSurDelta_r_nr(:,surIdx) = mean(frAux_r_sur)-mean(frAux_nr_sur);

        % Left vs right choice surrogate
        I = randperm(length(data(:,5)));
        frAux_lc_sur = surData(surData(I,5)==leftChoice);
        frAux_rc_sur = surData(surData(I,5)==rightChoice);
        allSurDelta_lc_rc(:,surIdx) = mean(frAux_lc_sur)-mean(frAux_rc_sur);
    end

    % --- Compare real data against surrogate thresholds ---
    % Visual vs Spatial
    deltaMean = mean(frAux_v)-mean(frAux_s);
    threshold_up   = mean(allSurDelta_v_s, 2)+2*std(allSurDelta_v_s,1, 2);
    threshold_down = mean(allSurDelta_v_s, 2)-2*std(allSurDelta_v_s,1, 2);
    if deltaMean > threshold_up || deltaMean < threshold_down
        neuronModulate_v_s = [neuronModulate_v_s; strCells(i,[1:3,6,8,9])];
    end    

    % Reward vs No Reward
    deltaMean = mean(frAux_r)-mean(frAux_nr);
    threshold_up   = mean(allSurDelta_r_nr, 2)+2*std(allSurDelta_r_nr,1, 2);
    threshold_down = mean(allSurDelta_r_nr, 2)-2*std(allSurDelta_r_nr,1, 2);
    if deltaMean > threshold_up || deltaMean < threshold_down
        neuronModulate_r_nr = [neuronModulate_r_nr; strCells(i,[1:3,6,8,9])];
    end   

    % Left vs Right
    deltaMean = mean(frAux_cl)-mean(frAux_cr);
    threshold_up   = mean(allSurDelta_lc_rc, 2)+2*std(allSurDelta_lc_rc,1, 2);
    threshold_down = mean(allSurDelta_lc_rc, 2)-2*std(allSurDelta_lc_rc,1, 2);
    if deltaMean > threshold_up || deltaMean < threshold_down
        neuronModulate_lc_rc  = [neuronModulate_lc_rc ; strCells(i,[1:3,6,8,9])];
    end
end

% Count modulated neurons
length(neuronModulate_v_s(:,1))
length(neuronModulate_r_nr(:,1))
length(neuronModulate_lc_rc(:,1))
%% FSI positive
clc
fsiToken = 22;
speedCellToken = 1;
[A,B,C,AB,AC,BC,ABC] = neuron_venn(neuronModulate_v_s, neuronModulate_r_nr, neuronModulate_lc_rc, fsiToken, +1, speedCellToken)

%% FSI negative
clc
[A,B,C,AB,AC,BC,ABC] = neuron_venn(neuronModulate_v_s, neuronModulate_r_nr, neuronModulate_lc_rc, fsiToken, -1, speedCellToken)

%% MSN positive
clc
msnToken = 21;
[A,B,C,AB,AC,BC,ABC] = neuron_venn(neuronModulate_v_s, neuronModulate_r_nr, neuronModulate_lc_rc, msnToken, +1, speedCellToken)

%% MSN negative
clc
[A,B,C,AB,AC,BC,ABC] = neuron_venn(neuronModulate_v_s, neuronModulate_r_nr, neuronModulate_lc_rc, msnToken, -1, speedCellToken)

%% ---------- Helper functions ----------
function [data] = readData(directory)
    filePattern = fullfile(directory, '*.mat'); % Change to whatever pattern you need
    theFiles = dir(filePattern);

    for LengthEachCell = 1 : length(theFiles)
        baseFileName = theFiles(LengthEachCell).name;
        fullFileName = fullfile(theFiles(LengthEachCell).folder, baseFileName);
        fprintf(1, 'Reading %s\n', fullFileName);
        data{LengthEachCell} = load(fullFileName);
    end
end


function [A,B,C,AB,AC,BC,ABC] = neuron_venn(neuronModulate_v_s, neuronModulate_r_nr, neuronModulate_lc_rc, token, signFlag, speedCellToken)
% neuron_venn: compute Venn counts for neuron modulation analysis
%
% Inputs:
%   neuronModulate_v_s   - cell array/matrix (v_s condition)
%   neuronModulate_r_nr  - cell array/matrix (r_nr condition)
%   neuronModulate_lc_rc - cell array/matrix (lc_rc condition)
%   token                - neuron type token (22 = FSI, 21 = MSN, etc.)
%   signFlag             - +1 for positive, -1 for negative correlation
%   speedCellToken       - token for speed cell (usually = 1)
%
% Outputs:
%   A   - count only in v_s
%   B   - count only in r_nr
%   C   - count only in lc_rc
%   AB  - count in v_s & r_nr only
%   AC  - count in v_s & lc_rc only
%   BC  - count in r_nr & lc_rc only
%   ABC - count in all three

    % 1. Define filter condition
    if signFlag > 0
        signCond = @(x) x > 0;
    else
        signCond = @(x) x < 0;
    end

    % 2. Filter data
    I = find(cell2mat(neuronModulate_v_s(:,4)) == token & ...
             signCond(cell2mat(neuronModulate_v_s(:,5))) & ...
             cell2mat(neuronModulate_v_s(:,6)) == speedCellToken);
    v_s = neuronModulate_v_s(I,:);

    I = find(cell2mat(neuronModulate_r_nr(:,4)) == token & ...
             signCond(cell2mat(neuronModulate_r_nr(:,5))) & ...
             cell2mat(neuronModulate_r_nr(:,6)) == speedCellToken);
    r_nr = neuronModulate_r_nr(I,:);

    I = find(cell2mat(neuronModulate_lc_rc(:,4)) == token & ...
             signCond(cell2mat(neuronModulate_lc_rc(:,5))) & ...
             cell2mat(neuronModulate_lc_rc(:,6)) == speedCellToken);
    lc_rc = neuronModulate_lc_rc(I,:);

    % 3. Create keys from first 3 columns
    format_key = @(row) sprintf('%s_%d_%d', row{1}, row{2}, row{3});
    keys_v_s   = unique(cellfun(format_key, num2cell(v_s(:,1:3), 2), 'UniformOutput', false));
    keys_r_nr  = unique(cellfun(format_key, num2cell(r_nr(:,1:3), 2), 'UniformOutput', false));
    keys_lc_rc = unique(cellfun(format_key, num2cell(lc_rc(:,1:3), 2), 'UniformOutput', false));

    % 4. Set operations
    only_v_s   = setdiff(keys_v_s, union(keys_r_nr, keys_lc_rc));
    only_r_nr  = setdiff(keys_r_nr, union(keys_v_s, keys_lc_rc));
    only_lc_rc = setdiff(keys_lc_rc, union(keys_v_s, keys_r_nr));

    v_s_and_r_nr   = intersect(keys_v_s, keys_r_nr);
    v_s_and_lc_rc  = intersect(keys_v_s, keys_lc_rc);
    r_nr_and_lc_rc = intersect(keys_r_nr, keys_lc_rc);

    % 5. Intersection of all three
    all_three = intersect(v_s_and_r_nr, keys_lc_rc);

    % 6. Subtract all_three from double intersections
    v_s_and_r_nr_only   = setdiff(v_s_and_r_nr, all_three);
    v_s_and_lc_rc_only  = setdiff(v_s_and_lc_rc, all_three);
    r_nr_and_lc_rc_only = setdiff(r_nr_and_lc_rc, all_three);

    % 7. Counts
    A   = length(only_v_s);
    B   = length(only_r_nr);
    C   = length(only_lc_rc);
    AB  = length(v_s_and_r_nr_only);
    AC  = length(v_s_and_lc_rc_only);
    BC  = length(r_nr_and_lc_rc_only);
    ABC = length(all_three);
end