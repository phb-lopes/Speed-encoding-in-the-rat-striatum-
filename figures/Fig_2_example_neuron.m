%% Initialization
clc, clf, clear
timeInfoSession = 1;   % Size of time window (s)
lag = 3;               % Maximum lag for CCG analysis (s)
overlap = .9;          % Fraction of overlap between windows
step = timeInfoSession - overlap;   % Step size for sliding window
if overlap == 0
    step = timeInfoSession;
end

% Define folder path based on parameters
folder = strcat('/Users/phb_lopes/Library/CloudStorage/OneDrive-Pessoal/projetoEletrofisiologiaMacbook/', ...
    num2str(timeInfoSession),'sInterval/', num2str(lag),'sLag/',num2str(overlap),'_overlap/');

% Load all data files
dataFsiAndLocomotionCells_all = readData(folder);

% Move to base project directory
cd /Users/phb_lopes/Library/CloudStorage/OneDrive-Pessoal/projetoEletrofisiologiaMacbook


%% Collect all neurons into a single structure
clc
cells = {};
for i = 1:length(dataFsiAndLocomotionCells_all)
    cells = [cells; dataFsiAndLocomotionCells_all{1,i}.data];
end


%% Example neuron (Figure 2)
neuron = 114;   % Example neuron index
data = cells(neuron,:);
dataFR = data{4}(:,1);  % Firing rate (Hz)
dataVel = data{5};      % Locomotion speed (cm/s)


%% Figure setup
clc, clf, clear allSurCCG
fig1 = gcf;
set(gcf,'color','w')
set(fig1,'position',[0, 0, 1400, 1050]);   % Figure size

% Subplot grid
rows = 5;
cols = 3;

% Plotting settings
std_gaussian = .5;
colorScatter = '#7575a3';
colorCircleScatter = '#7575a3';
time_firing_rate = (1:length(dataFR)) * step;   % Time axis

% Time window for display
ini = 350;
fim = 550;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(rows, cols, [1:3])
plot(time_firing_rate, dataVel, 'color','#e6ac00','lineW',3);
xticks([])
xlim([ini fim])
ylabel({'Locomotion','speed (cm/s)'});
set(gca, 'FontWeight', 'bold','FontSize', 12);
box off
text(-.07, 1.2, 'A', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(rows, cols, [4:6])
plot(time_firing_rate, dataFR, 'color','#8300b4','lineW',3);
xticks([])
xlim([ini fim])
ylabel({'Firing rate (Hz)'});
set(gca, 'FontWeight', 'bold','FontSize', 12);
box off
text(-.07, 1.2, 'B', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(rows, cols, [7:9])
hold on
plot(time_firing_rate, zscore(dataFR), 'color','#8300b4','lineW',3);
plot(time_firing_rate, zscore(dataVel), 'color','#e6ac00','lineW',3);
xlim([ini fim])
ylabel('Z-score');
xlabel('Time (s)');
xticks([ini:25:fim])
xticklabels([0:25:(fim-ini)])
set(gca, 'FontWeight', 'bold','FontSize', 12);
text(-.07, 1.2, 'C', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(rows, cols, [[10 11],[13 14]])
hold on
n = length(dataVel);
sample_size = round(0.15 * n);   % Use only 15% of samples for clarity
sample_idx = randperm(n, sample_size);
vel_sampled = dataVel(sample_idx);
fr_sampled = dataFR(sample_idx);

scatter(vel_sampled, fr_sampled, 'MarkerEdgeColor',colorScatter,'MarkerFaceColor',colorCircleScatter)

% Fit (linear regression)
[y_fit, R2_1] = getFit_polyfit(dataVel, dataFR', 1);

% Pearson correlation
[r, p] = coeficiente_pearson(dataVel, dataFR);
text(55,55,strcat(['r = ' num2str(round(r,3))]), 'color','k','FontWeight', 'bold','FontSize', 16)

plot(dataVel, y_fit, 'k-', 'LineWidth', 3);
ylim([0 60])
ylabel('Firing rate (Hz)');
xlabel({'Locomotion speed (cm/s)'});
set(gca, 'FontWeight', 'bold','FontSize', 12);
text(-.11, 1.05, 'D', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(rows, cols, [[12],[15]])
surNumber = 10000;       % Number of surrogate shifts
bin_step = step;         % Bin size (s)

maxLag_seconds = lag;
maxLag_bins = round(maxLag_seconds/bin_step);  % Convert lag to bins

% Generate surrogate CCGs by circular shifting
for surIdx = 1:surNumber
    shiftAmount = randi([20,50],1)/100;
    circularShift = round(length(dataVel)*shiftAmount,0);
    [ccg, ~] = xcorr(zscore(circshift(dataFR',circularShift)), zscore(dataVel), maxLag_bins, 'normalized');
    allSurCCG(:,surIdx) = ccg;
end

% Real CCG
hold on;
[ccg, lags] = xcorr(zscore(dataFR), zscore(dataVel), maxLag_bins, 'coeff');

% Threshold for significance
Threshold = max(mean(allSurCCG, 2) + 2*std(allSurCCG,1, 2));

valoresLag_ms = lags*bin_step;

plot(valoresLag_ms, ccg,'k','linew',3)
plot(valoresLag_ms, mean(allSurCCG, 2),'r-','linew',3)

% Shaded area = ± 2 SD of surrogate distribution
A = mean(allSurCCG, 2)+2*std(allSurCCG,1, 2);
B = mean(allSurCCG, 2)-2*std(allSurCCG,1, 2);   
fill([valoresLag_ms, fliplr(valoresLag_ms)], [A' , fliplr(B')], 'r', ...
    'facealpha',0.2,'EdgeColor','none', 'HandleVisibility', 'off')

plot([0 0], [-.2 .6], 'color','#9F9E9E','lineS','--','linew',2, 'HandleVisibility', 'off')
xlabel('Lag (s)')
ylabel('CCG (r)')
set(gca, 'FontWeight', 'bold','FontSize', 13);
xticks(-2.5:.5:2.5)
xlim([-2.5 2.5])
yticks([-.2:.1:.6])
ylim([-.2 .6])
legend('Real','Sur (± 2 SD)', 'FontWeight', 'bold', 'FontSize', 11)
set(legend,'Position',[0.81 0.38 0.1 0.01])
box off
text(-.2, 1.05, 'E', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

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

function [y_fit, R2] = getFit_polyfit(data_vel, data_fr, degree)
    % Polynomial fit between velocity and firing rate
    if isrow(data_vel), data_vel = data_vel'; end
    if isrow(data_fr), data_fr = data_fr'; end
    
    % Polynomial fit
    p = polyfit(data_vel, data_fr, degree);
    y_fit = polyval(p, data_vel);
    
    % R² computation
    SSR = sum((data_fr - y_fit).^2);
    SST = sum((data_fr - mean(data_fr)).^2);
    R2 = 1 - (SSR / SST);
end

function [r, p] = coeficiente_pearson(data_vel, data_fr)
    % Pearson correlation coefficient
    if length(data_vel) ~= length(data_fr)
        error('Both variables must have the same length');
    end
    [r, p] = corr(data_vel(:), data_fr(:), 'Type', 'Pearson');
end
