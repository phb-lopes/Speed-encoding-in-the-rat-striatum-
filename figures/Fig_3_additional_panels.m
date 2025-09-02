%%% Initialization
clc, clf, clear
timeInfoSession = 1;   % Size of time window (s)
lag = 3;               % Maximum lag for CCG analysis (s)
overlap = .9;          % Fraction of overlap between windows
step = timeInfoSession - overlap;   % Step size for sliding window
if overlap == 0
    step = timeInfoSession;
end
bin_step = timeInfoSession-overlap;
maxLag_seconds = lag;

% Define folder path based on parameters
folder = strcat('/Users/phb_lopes/Library/CloudStorage/OneDrive-Pessoal/projetoEletrofisiologiaMacbook/', ...
    num2str(timeInfoSession),'sInterval/', num2str(lag),'sLag/',num2str(overlap),'_overlap/');

% Load all data files
dataFsiAndLocomotionCells_all = readData(folder);

% Move to base project directory
cd /Users/phb_lopes/Library/CloudStorage/OneDrive-Pessoal/projetoEletrofisiologiaMacbook


%% Collect all neurons into a single structure
clc
allCells = {};
for i = 1:length(dataFsiAndLocomotionCells_all)
    allCells = [allCells; dataFsiAndLocomotionCells_all{1,i}.data];
end

% Select cells with classification >= 21
strCells = allCells(cell2mat(allCells(:,6)) >=21,:);

% Load pre-computed CCG results (change filename if needed)
load('./ccgLagFrSpeedCellsInfoPosAndNeg_1sInterval_0.9overlap_3sLagNoSCwith2.5Lag.mat') % Change to whatever pattern you need
valorCcg = ccgLagFrLocomotionCellsData{1};
valoresLag = ccgLagFrLocomotionCellsData{2};
valoresFr = ccgLagFrLocomotionCellsData{3};
speedCell = ccgLagFrLocomotionCellsData{4};
lag_ms = valoresLag*bin_step;
bin_step = timeInfoSession-overlap;           % Bin size in ms
maxLag_seconds = lag;
maxLag_bins = round(maxLag_seconds/bin_step);  
valoresLag_ms = valoresLag*bin_step;  % convert bin position to ms
binHisto = -maxLag_seconds:bin_step:maxLag_seconds;
%% Get CCG from positive and negative cell
maxLag_seconds = lag;
maxLag_bins = round(maxLag_seconds/bin_step);
positive = strCells(28,:);
dataFR = positive{4}(:,1);
dataVel = positive{5};
[ccg_pos, lags] = xcorr(zscore(dataFR), zscore(dataVel), maxLag_bins, 'coeff');

negative = strCells(100,:);
dataFR = negative{4}(:,1);
dataVel = negative{5};
[ccg_neg, lags] = xcorr(zscore(dataFR), zscore(dataVel), maxLag_bins, 'coeff');
%%
clf, clc
fig1 = gcf;
set(gcf,'color','w')
set(fig1,'position',[0, 0, 1300, 1050]);
pos1 = [0.13,0.489517819706499,0.775,0.435482180293501];
colorCirclScatterSignificant = '#37a259';
colorScatterSignificant = '#37a259';
lags_ms = lags*bin_step;

% Positive speed cell
subplot(2,2,1)
[~, idx] = max(abs(ccg_pos));
hold on;
plot(lags_ms, ccg_pos, 'k','linew',3)
xticks(-2.5:.5:2.5)
scatter(lags_ms(idx), ccg_pos(idx),1000,'MarkerEdgeColor',colorScatterSignificant,'MarkerFaceColor',colorCirclScatterSignificant);
yticks(0:.2:.4)
xlim([-2.5 2.5])
xticks([])
ylim([0 .4])
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 35);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 35);
box off

% Negative speed cell
subplot(2,2,3)
[~, idx] = max(abs(ccg_neg));
hold on;
plot(lags*bin_step, ccg_neg, 'k','linew',3)
scatter(lags_ms(idx), ccg_neg(idx),1000,'MarkerEdgeColor',colorScatterSignificant,'MarkerFaceColor',colorCirclScatterSignificant);
xticks(-2.5:2.5:2.5)
xlim([-2.5 2.5])
yticks(-.4:.2:0)
ylim([-.4 0])
xlabel('Lag (s)')
ax = gca;
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 35);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 35);
xticklabel = ax.XTickLabel;  % Get the tick labels
xtickangle(0);
box off
f = gcf;
%% ---------- Helper functions ----------
function [data] = readData(directory)
    % Reads all .mat files in a folder into a cell array
    filePattern = fullfile(directory, '*.mat');
    theFiles = dir(filePattern);

    for LengthEachCell = 1 : length(theFiles)
        baseFileName = theFiles(LengthEachCell).name;
        fullFileName = fullfile(theFiles(LengthEachCell).folder, baseFileName);
        fprintf(1, 'Reading %s\n', fullFileName);
        data{LengthEachCell} = load(fullFileName);
    end
end