%%% Initialization
clc, clf, clear
timeInfoSession = 1;   % Size of time window (s)
lag = 3;               % Maximum lag for CCG analysis (s)
overlap = .9;          % Fraction of overlap between windows
step = timeInfoSession - overlap;   % Step size for sliding window
if overlap == 0
    step = timeInfoSession;
end
bin_step = timeInfoSession-overlap; % Bin size in seconds
maxLag_seconds = lag;               % Maximum lag in seconds

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
load('./ccgLagFrSpeedCellsInfoPosAndNeg_1sInterval_0.9overlap_3sLagNoSCwith2.5Lag.mat') % Change if filename differs
valorCcg = ccgLagFrLocomotionCellsData{1};
valoresLag = ccgLagFrLocomotionCellsData{2};
valoresFr = ccgLagFrLocomotionCellsData{3};
speedCell = ccgLagFrLocomotionCellsData{4};

%%
clc
maxLag = lag;
% std_gaussiana = .5;
tic
numPortion = 5;

valorCcg_portion = {};
valoresLag_portion = {};

valorCcg_taskTrials = {};
valoresLag_taskTrials = {};

valorCcg_mazeSide = {};
valoresLag_mazeSide = {};

valorCcg_Rew = {};
valoresLag_Rew = {};

for neuron = 1:length(strCells)
    data = strCells(neuron,:);
    maxTrial = max(data{4}(:,6));

    part = floor(maxTrial / numPortion);
    rest = mod(maxTrial, numPortion);
    matrix_vectors = cell(numPortion, 1);
    
    beg = 1;
    for i = 1:numPortion
        if i <= rest
            end_ = beg + part;
        else
            end_ = beg + part - 1;
        end
        matrix_vectors{i} = beg:end_;
        beg = end_ + 1;
    end
    
    for portion = 1:numPortion
        dataFR = data{4}(ismember(data{4}(:,6), matrix_vectors{portion}),1);
        dataVel = data{5}(ismember(data{4}(:,6), matrix_vectors{portion}),1);
        time_firing_rate = 1:1:length(dataFR);

        [ccg, lags] = xcorr(zscore(dataFR), zscore(dataVel), maxLag, 'coeff');
        [m i] = max(abs(ccg));
        valorCcg_portion{neuron,portion} = ccg(i);
    end
    
    % Trial outcome
    I = find(data{4}(:,3) == 0); % Wrong
    dataFR = data{4}(I);
    dataVel = data{5}(I);
    time_firing_rate = 1:1:length(dataFR);

    dataFR = dataFR;
    dataVel = dataVel;

    [ccg, lags] = xcorr(zscore(dataFR), zscore(dataVel), maxLag, 'coeff');
    [m i] = max(abs(ccg));
    valorCcg_Rew{neuron,1} = ccg(i);

    
    I = find(data{4}(:,3) == 1); % Right
    dataFR = data{4}(I);
    dataVel = data{5}(I);
    time_firing_rate = 1:1:length(dataFR);
    dataFR = dataFR;
    dataVel = dataVel;

    [ccg, lags] = xcorr(zscore(dataFR), zscore(dataVel), maxLag, 'coeff');
    [m i] = max(abs(ccg));
    valorCcg_Rew{neuron,2} = ccg(i);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cue
    I = find(data{4}(:,4) == 0);    % Visual
    dataFR = data{4}(I);
    dataVel = data{5}(I);
    time_firing_rate = 1:1:length(dataFR);
    dataFR = dataFR;
    dataVel = dataVel;

    [ccg, lags] = xcorr(zscore(dataFR), zscore(dataVel), maxLag, 'coeff');
    [m i] = max(abs(ccg));
    valorCcg_taskTrials{neuron,1} = ccg(i);
    
    I = find(data{4}(:,4) == 1);    % Spatial
    dataFR = data{4}(I);
    dataVel = data{5}(I);
    time_firing_rate = 1:1:length(dataFR);
    dataFR = dataFR;
    dataVel = dataVel;

    [ccg, lags] = xcorr(zscore(dataFR), zscore(dataVel), maxLag, 'coeff');
    [m i] = max(abs(ccg));
    valorCcg_taskTrials{neuron,2} = ccg(i);
    
    
    % Spatial choice
    I = find(data{4}(:,5) == 0); % left
    dataFR = data{4}(I);
    dataVel = data{5}(I);
    time_firing_rate = 1:1:length(dataFR);
    dataFR = dataFR;
    dataVel = dataVel;

    [ccg, lags] = xcorr(zscore(dataFR), zscore(dataVel), maxLag, 'coeff');
    [m i] = max(abs(ccg));
    valorCcg_mazeSide{neuron,1} = ccg(i);
    
    I = find(data{4}(:,5) == 1); % right
    dataFR = data{4}(I);
    dataVel = data{5}(I);
    time_firing_rate = 1:1:length(dataFR);

    dataFR = dataFR;
    dataVel = dataVel;

    [ccg, lags] = xcorr(zscore(dataFR), zscore(dataVel), maxLag, 'coeff');
    [m i] = max(abs(ccg));
    valorCcg_mazeSide{neuron,2} = ccg(i);
    
end
%% 
clc, clf
fig1 = gcf;
set(gcf,'color','w')
set(fig1,'position',[0, 0, 1200, 1050]);
pos1 = [0.13,0.489517819706499,0.775,0.435482180293501];
lines = 6;
columns = 2;
pValueThreshold = .01;
colorCirclScatterFSI = '#FC4366';
colorScatterFSI = '#FC4366';

colorCirclScatterMSN = '#AEB2FF';
colorScatterMSN = '#AEB2FF';

FSI_pos = (cell2mat(strCells(:,6)) == 22)' & (speedCell==1) & ccgLagFrLocomotionCellsData{1}>0;
MSN_pos = (cell2mat(strCells(:,6)) == 21)' & (speedCell==1) & ccgLagFrLocomotionCellsData{1}>0;

FSI_neg = (cell2mat(strCells(:,6)) == 22)' & (speedCell==1) & ccgLagFrLocomotionCellsData{1}<0;
MSN_neg = (cell2mat(strCells(:,6)) == 21)' & (speedCell==1) & ccgLagFrLocomotionCellsData{1}<0;

subplot(lines, columns, [1, 3, 5])
hold on;
% Positivos
dataFSI_pos = cell2mat(valorCcg_portion(FSI_pos,:));
scatter(dataFSI_pos, [1:5]+.2, 100, [.988 .263 .3], 's', 'filled', 'MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
plot([mean(ccgLagFrLocomotionCellsData{1}(:,FSI_pos)) mean(ccgLagFrLocomotionCellsData{1}(:,FSI_pos))], [0 6.5],'color',[.988 .263 .4],'lineS','--', 'lineW',2, 'HandleVisibility', 'off')
scatter(mean(cell2mat(valorCcg_portion(FSI_pos,:))',2), [1:length(mean(cell2mat(valorCcg_portion(FSI_pos,:))))]+.1,300, 's', 'MarkerEdgeColor', colorScatterFSI, 'MarkerFaceColor', colorCirclScatterFSI)
% Negativos
dataFSI_neg = cell2mat(valorCcg_portion(FSI_neg,:));
scatter(dataFSI_neg, [1:5]-.2, 100, [.988 .263 .3], 'o', 'filled', 'MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
plot([mean(ccgLagFrLocomotionCellsData{1}(:,FSI_neg)) mean(ccgLagFrLocomotionCellsData{1}(:,FSI_neg))], [0 6.5],'color',[.988 .263 .4],'lineS','--', 'lineW',2, 'HandleVisibility', 'off')
scatter(mean(cell2mat(valorCcg_portion(FSI_neg,:))',2), [1:length(mean(cell2mat(valorCcg_portion(FSI_neg,:))))]-.1, 300, 'o', 'MarkerEdgeColor', colorScatterFSI, 'MarkerFaceColor', colorCirclScatterFSI)

plot([0 0], [0 6.5], 'color',[0 0 0  0.2],'lineS','-', 'lineW',2, 'HandleVisibility', 'off')

[h, p, ci, stats] = ttest(dataFSI_pos, mean(ccgLagFrLocomotionCellsData{1}(:,FSI_pos)), 'alpha', pValueThreshold/5);
p < pValueThreshold/5;
errorbar(mean(dataFSI_pos), [1:length(mean(cell2mat(valorCcg_portion(FSI_pos,:))))]+.1, ...
    (ci(2,1) - ci(1,:))/2, ...
    'horizontal', 'k', 'LineStyle', 'none','lineW',1.5, 'HandleVisibility', 'off');

[h, p, ci, stats] = ttest(dataFSI_neg, mean(ccgLagFrLocomotionCellsData{1}(:,FSI_neg)), 'alpha', pValueThreshold/5);
p < pValueThreshold/5;
errorbar(mean(dataFSI_neg), [1:length(mean(cell2mat(valorCcg_portion(FSI_neg,:))))]-.1, ...
    (ci(2,1) - ci(1,:))/2, ...
    'horizontal', 'k', 'LineStyle', 'none','lineW',1.5, 'HandleVisibility', 'off');
ylim([.8 5.2])
yticks(1:5)
yticklabels({'1st','2nd','3rd','4th','5th'})
xlabel([])
ylabel('Trial block')
set(gca, 'FontWeight', 'bold', 'FontSize', 13);
xlim([-.7, .7])
xticks([])
legend({'Positive FSI', 'Negative FSI'}, 'FontWeight', 'bold', 'FontSize', 15, 'NumColumns', 1);
set(legend,'Position', [0.177 0.917 0.03 0.02], 'FontWeight','bold');
text(-.1, 1.06, 'A', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);


subplot(lines, columns, 2)
hold on;
% Positivos
dataFSI_pos = cell2mat(valorCcg_taskTrials(FSI_pos,:));
scatter(dataFSI_pos, [1:2]+.2, 100, [.988 .263 .3], 's', 'filled','MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
plot([mean(ccgLagFrLocomotionCellsData{1}(:,FSI_pos)) mean(ccgLagFrLocomotionCellsData{1}(:,FSI_pos))],[0 6.5], 'color', [.988 .263 .4], 'lineS', '--', 'lineW', 2)
scatter(mean(cell2mat(valorCcg_taskTrials(FSI_pos,:))',2), [1:length(mean(cell2mat(valorCcg_taskTrials(FSI_pos,:))))]+.1, 300, 's', 'MarkerEdgeColor', colorScatterFSI, 'MarkerFaceColor', colorCirclScatterFSI)
% Negativos
dataFSI_neg = cell2mat(valorCcg_taskTrials(FSI_neg,:));
scatter(dataFSI_neg, [1:2]-.2, 100, [.988 .263 .3], 'o', 'filled','MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
plot([mean(ccgLagFrLocomotionCellsData{1}(:,FSI_neg)) mean(ccgLagFrLocomotionCellsData{1}(:,FSI_neg))],[0 6.5], 'color', [.988 .263 .4], 'lineS', '--', 'lineW', 2)
scatter(mean(cell2mat(valorCcg_taskTrials(FSI_neg,:))',2), [1:length(mean(cell2mat(valorCcg_taskTrials(FSI_neg,:))))]-.1, 300, 'o', 'MarkerEdgeColor', colorScatterFSI, 'MarkerFaceColor', colorCirclScatterFSI)

plot([0 0], [0 6.5], 'color',[0 0 0  0.2],'lineS','-', 'lineW',2, 'HandleVisibility', 'off')

[h, p, ci, stats] = ttest(dataFSI_pos, mean(ccgLagFrLocomotionCellsData{1}(:,FSI_pos)), ...
    'alpha', pValueThreshold/2);
p < pValueThreshold/2;
errorbar(mean(dataFSI_pos), ...
    [1:length(mean(cell2mat(valorCcg_taskTrials(FSI_pos,:))))]+.1, ...
    (ci(2,1) - ci(1,:))/2, ...
    'horizontal', 'k', 'LineStyle', 'none', 'lineW', 1.5, 'HandleVisibility', 'off');

[h, p, ci, stats] = ttest(dataFSI_neg, mean(ccgLagFrLocomotionCellsData{1}(:,FSI_neg)), ...
    'alpha', pValueThreshold/2);
p < pValueThreshold/2;
errorbar(mean(dataFSI_neg), ...
    [1:length(mean(cell2mat(valorCcg_taskTrials(FSI_neg,:))))]-.1, ...
    (ci(2,1) - ci(1,:))/2, ...
    'horizontal', 'k', 'LineStyle', 'none', 'lineW', 1.5, 'HandleVisibility', 'off');

ylim([.8 2.2])
yticks(1:2)
yticklabels({'Visual','Spatial'})
xlabel([])
ylabel('Cue')
set(gca, 'FontWeight', 'bold', 'FontSize', 13);
xlim([-.6, .6])
xticks([])
text(-.25, 1.2, 'B', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);


subplot(lines, columns,4)
hold on;
% Positivos
dataFSI_pos = cell2mat(valorCcg_mazeSide(FSI_pos,:));
scatter(dataFSI_pos, [1:2]+.2, 100, [.988 .263 .3], 's', 'filled','MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
plot([mean(ccgLagFrLocomotionCellsData{1}(:,FSI_pos)) mean(ccgLagFrLocomotionCellsData{1}(:,FSI_pos))],[0 6.5], 'color', [.988 .263 .4], 'lineS', '--', 'lineW', 2)
scatter(mean(cell2mat(valorCcg_mazeSide(FSI_pos,:))',2), [1:length(mean(cell2mat(valorCcg_mazeSide(FSI_pos,:))))]+.1, 300, 's', 'MarkerEdgeColor', colorScatterFSI, 'MarkerFaceColor', colorCirclScatterFSI)
% Negativos
dataFSI_neg = cell2mat(valorCcg_mazeSide(FSI_neg,:));
scatter(dataFSI_neg, [1:2]-.2, 100, [.988 .263 .3], 'o', 'filled','MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
plot([mean(ccgLagFrLocomotionCellsData{1}(:,FSI_neg)) mean(ccgLagFrLocomotionCellsData{1}(:,FSI_neg))],[0 6.5], 'color', [.988 .263 .4], 'lineS', '--', 'lineW', 2)
scatter(mean(cell2mat(valorCcg_mazeSide(FSI_neg,:))',2), [1:length(mean(cell2mat(valorCcg_mazeSide(FSI_neg,:))))]-.1, 300, 'o', 'MarkerEdgeColor', colorScatterFSI, 'MarkerFaceColor', colorCirclScatterFSI)

plot([0 0], [0 6.5], 'color',[0 0 0  0.2],'lineS','-', 'lineW',2, 'HandleVisibility', 'off')

[h, p, ci, stats] = ttest(dataFSI_pos, mean(ccgLagFrLocomotionCellsData{1}(:,FSI_pos)), ...
    'alpha', pValueThreshold/2);
p < pValueThreshold/2;
errorbar(mean(dataFSI_pos), ...
    [1:length(mean(cell2mat(valorCcg_mazeSide(FSI_pos,:))))]+.1, ...
    (ci(2,1) - ci(1,:))/2, ...
    'horizontal', 'k', 'LineStyle', 'none', 'lineW', 1.5, 'HandleVisibility', 'off');

[h, p, ci, stats] = ttest(dataFSI_neg, mean(ccgLagFrLocomotionCellsData{1}(:,FSI_neg)), ...
    'alpha', pValueThreshold/2);
p < pValueThreshold/2;
errorbar(mean(dataFSI_neg), ...
    [1:length(mean(cell2mat(valorCcg_mazeSide(FSI_neg,:))))]-.1, ...
    (ci(2,1) - ci(1,:))/2, ...
    'horizontal', 'k', 'LineStyle', 'none', 'lineW', 1.5, 'HandleVisibility', 'off');

ylim([.8 2.2])
yticks(1:2)
yticklabels({'Left','Right'})
xlabel([])
ylabel('Spatial choice')
set(gca, 'FontWeight', 'bold', 'FontSize', 13);
xlim([-.6, .6])
xticks([])
text(-.25, 1.2, 'C', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

subplot(lines,columns,6)

hold on;
% Positivos
dataFSI_pos = cell2mat(valorCcg_Rew(FSI_pos,:));
scatter(dataFSI_pos, [1:2]+.2, 100, [.988 .263 .3], 's', 'filled','MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
plot([mean(ccgLagFrLocomotionCellsData{1}(:,FSI_pos)) mean(ccgLagFrLocomotionCellsData{1}(:,FSI_pos))],[0 6.5], 'color', [.988 .263 .4], 'lineS', '--', 'lineW', 2)
scatter(mean(cell2mat(valorCcg_Rew(FSI_pos,:))',2), [1:length(mean(cell2mat(valorCcg_Rew(FSI_pos,:))))]+.1, 300, 's', 'MarkerEdgeColor', colorScatterFSI, 'MarkerFaceColor', colorCirclScatterFSI)
% Negativos
dataFSI_neg = cell2mat(valorCcg_Rew(FSI_neg,:));
scatter(dataFSI_neg, [1:2]-.2, 100, [.988 .263 .3], 'o', 'filled','MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
plot([mean(ccgLagFrLocomotionCellsData{1}(:,FSI_neg)) mean(ccgLagFrLocomotionCellsData{1}(:,FSI_neg))],[0 6.5], 'color', [.988 .263 .4], 'lineS', '--', 'lineW', 2)
scatter(mean(cell2mat(valorCcg_Rew(FSI_neg,:))',2), [1:length(mean(cell2mat(valorCcg_Rew(FSI_neg,:))))]-.1, 300, 'o', 'MarkerEdgeColor', colorScatterFSI, 'MarkerFaceColor', colorCirclScatterFSI)

plot([0 0], [0 6.5], 'color',[0 0 0  0.2],'lineS','-', 'lineW',2, 'HandleVisibility', 'off')

[h, p, ci, stats] = ttest(dataFSI_pos, mean(ccgLagFrLocomotionCellsData{1}(:,FSI_pos)), ...
    'alpha', pValueThreshold/2);
p < pValueThreshold/2;
errorbar(mean(dataFSI_pos), ...
    [1:length(mean(cell2mat(valorCcg_Rew(FSI_pos,:))))]+.1, ...
    (ci(2,1) - ci(1,:))/2, ...
    'horizontal', 'k', 'LineStyle', 'none', 'lineW', 1.5, 'HandleVisibility', 'off');

[h, p, ci, stats] = ttest(dataFSI_neg, mean(ccgLagFrLocomotionCellsData{1}(:,FSI_neg)), ...
    'alpha', pValueThreshold/2);
p < pValueThreshold/2;
errorbar(mean(dataFSI_neg), ...
    [1:length(mean(cell2mat(valorCcg_Rew(FSI_neg,:))))]-.1, ...
    (ci(2,1) - ci(1,:))/2, ...
    'horizontal', 'k', 'LineStyle', 'none', 'lineW', 1.5, 'HandleVisibility', 'off');

ylim([.8 2.2])
yticks(1:2)
yticklabels({'No\newlineReward','Reward'})
xlabel([])
ylabel('Trial outcome')
set(gca, 'FontWeight', 'bold', 'FontSize', 13);
xlim([-.6, .6])
xticks([])
text(-.25, 1.2, 'D', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(lines, columns, [7, 9, 11])
hold on;
% Positivos
dataMSN_pos = cell2mat(valorCcg_portion(MSN_pos,:));
h1 = scatter(mean(cell2mat(valorCcg_portion(MSN_pos,:))',2), [1:length(mean(cell2mat(valorCcg_portion(MSN_pos,:))))]+.1, 300, 's', 'MarkerEdgeColor', colorScatterMSN, 'MarkerFaceColor', colorCirclScatterMSN);
scatter(dataMSN_pos, [1:5]+.2, 100, [.682 .698 1], 's', 'filled', 'MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
plot([mean(ccgLagFrLocomotionCellsData{1}(:,MSN_pos)) mean(ccgLagFrLocomotionCellsData{1}(:,MSN_pos))], [0 6.5],'color',[.682 .698 1],'lineS','--', 'lineW',2, 'HandleVisibility', 'off')
% Negativos
dataMSN_neg = cell2mat(valorCcg_portion(MSN_neg,:));
h2 = scatter(dataMSN_neg, [1:5]-.2, 100, [.682 .698 1], 'o', 'filled', 'MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
plot([mean(ccgLagFrLocomotionCellsData{1}(:,MSN_neg)) mean(ccgLagFrLocomotionCellsData{1}(:,MSN_neg))], [0 6.5],'color',[.682 .698 1],'lineS','--', 'lineW',2, 'HandleVisibility', 'off')
scatter(mean(cell2mat(valorCcg_portion(MSN_neg,:))',2), [1:length(mean(cell2mat(valorCcg_portion(MSN_neg,:))))]-.1, 300, 'o', 'MarkerEdgeColor', colorScatterMSN, 'MarkerFaceColor', colorCirclScatterMSN);

plot([0 0], [0 6.5], 'color',[0 0 0  0.2],'lineS','-', 'lineW',2, 'HandleVisibility', 'off')

[h, p, ci, stats] = ttest(dataMSN_pos, mean(ccgLagFrLocomotionCellsData{1}(:,MSN_pos)), 'alpha', pValueThreshold/5);
p < pValueThreshold/5;
errorbar(mean(dataMSN_pos), [1:length(mean(cell2mat(valorCcg_portion(MSN_pos,:))))]+.1, ...
    (ci(2,1) - ci(1,:))/2, ...
    'horizontal', 'k', 'LineStyle', 'none','lineW',1.5, 'HandleVisibility', 'off');

[h, p, ci, stats] = ttest(dataMSN_neg, mean(ccgLagFrLocomotionCellsData{1}(:,MSN_neg)), 'alpha', pValueThreshold/5);
p < pValueThreshold/5;
errorbar(mean(dataMSN_neg), [1:length(mean(cell2mat(valorCcg_portion(MSN_neg,:))))]-.1, ...
    (ci(2,1) - ci(1,:))/2, ...
    'horizontal', 'k', 'LineStyle', 'none','lineW',1.5, 'HandleVisibility', 'off');
ylim([.8 5.2])
yticks(1:5)
yticklabels({'1st','2nd','3rd','4th','5th'})
xlabel('Speed score (r)')
ylabel('Trial block')
xlim([-.7, .7])
xticks(-.7:.1:.7)
set(gca, 'FontWeight', 'bold', 'FontSize', 13);
ax = gca;
xticklabel = ax.XTickLabel;  % Get the tick labels
xtickangle(0);

legend([h1 h2], {'Positive MSN', 'Negative MSN'}, 'FontWeight', 'bold', 'FontSize', 15, 'NumColumns', 1)
set(legend,'Position', [0.18 0.49 0.03 0.02], 'FontWeight','bold')

text(-.1, 1.05, 'G', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

subplot(lines, columns,8)
hold on;
% Positivos
dataMSN_pos = cell2mat(valorCcg_taskTrials(MSN_pos,:));
scatter(dataMSN_pos, [1:2]+.2, 100, [.682 .698 1], 's', 'filled','MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
plot([mean(ccgLagFrLocomotionCellsData{1}(:,MSN_pos)) mean(ccgLagFrLocomotionCellsData{1}(:,MSN_pos))],[0 6.5], 'color', [.682 .698 1], 'lineS', '--', 'lineW', 2)
scatter(mean(cell2mat(valorCcg_taskTrials(MSN_pos,:))',2), [1:length(mean(cell2mat(valorCcg_taskTrials(MSN_pos,:))))]+.1, 300, 's', 'MarkerEdgeColor', colorScatterMSN, 'MarkerFaceColor', colorCirclScatterMSN)
% Negativos
dataMSN_neg = cell2mat(valorCcg_taskTrials(MSN_neg,:));
scatter(dataMSN_neg, [1:2]-.2, 100, [.682 .698 1], 'o', 'filled','MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
plot([mean(ccgLagFrLocomotionCellsData{1}(:,MSN_neg)) mean(ccgLagFrLocomotionCellsData{1}(:,MSN_neg))],[0 6.5], 'color', [.682 .698 1], 'lineS', '--', 'lineW', 2)
scatter(mean(cell2mat(valorCcg_taskTrials(MSN_neg,:))',2), [1:length(mean(cell2mat(valorCcg_taskTrials(MSN_neg,:))))]-.1, 300, 'o', 'MarkerEdgeColor', colorScatterMSN, 'MarkerFaceColor', colorCirclScatterMSN)

plot([0 0], [0 6.5], 'color',[0 0 0  0.2],'lineS','-', 'lineW',2, 'HandleVisibility', 'off')

[h, p, ci, stats] = ttest(dataMSN_pos, mean(ccgLagFrLocomotionCellsData{1}(:,MSN_pos)), ...
    'alpha', pValueThreshold/2);
p < pValueThreshold/2;
errorbar(mean(dataMSN_pos), ...
    [1:length(mean(cell2mat(valorCcg_taskTrials(MSN_pos,:))))]+.1, ...
    (ci(2,1) - ci(1,:))/2, ...
    'horizontal', 'k', 'LineStyle', 'none', 'lineW', 1.5, 'HandleVisibility', 'off');

[h, p, ci, stats] = ttest(dataMSN_neg, mean(ccgLagFrLocomotionCellsData{1}(:,MSN_neg)), ...
    'alpha', pValueThreshold/2);
p < pValueThreshold/2;
errorbar(mean(dataMSN_neg), ...
    [1:length(mean(cell2mat(valorCcg_taskTrials(MSN_neg,:))))]-.1, ...
    (ci(2,1) - ci(1,:))/2, ...
    'horizontal', 'k', 'LineStyle', 'none', 'lineW', 1.5, 'HandleVisibility', 'off');

ylim([.8 2.2])
yticks(1:2)
yticklabels({'Visual','Spatial'})
xlabel([])
ylabel('Cue')
set(gca, 'FontWeight', 'bold', 'FontSize', 13);
xlim([-.6, .6])
xticks([])
text(-.25, 1.2, 'H', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

subplot(lines, columns,10)
hold on;
% Positivos
dataMSN_pos = cell2mat(valorCcg_mazeSide(MSN_pos,:));
scatter(dataMSN_pos, [1:2]+.2, 100, [.682 .698 1], 's', 'filled','MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
plot([mean(ccgLagFrLocomotionCellsData{1}(:,MSN_pos)) mean(ccgLagFrLocomotionCellsData{1}(:,MSN_pos))],[0 6.5], 'color', [.682 .698 1], 'lineS', '--', 'lineW', 2)
scatter(mean(cell2mat(valorCcg_mazeSide(MSN_pos,:))',2), [1:length(mean(cell2mat(valorCcg_mazeSide(MSN_pos,:))))]+.1, 300, 's', 'MarkerEdgeColor', colorScatterMSN, 'MarkerFaceColor', colorCirclScatterMSN)
% Negativos
dataMSN_neg = cell2mat(valorCcg_mazeSide(MSN_neg,:));
scatter(dataMSN_neg, [1:2]-.2, 100, [.682 .698 1], 'o', 'filled','MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
plot([mean(ccgLagFrLocomotionCellsData{1}(:,MSN_neg)) mean(ccgLagFrLocomotionCellsData{1}(:,MSN_neg))],[0 6.5], 'color', [.682 .698 1], 'lineS', '--', 'lineW', 2)
scatter(mean(cell2mat(valorCcg_mazeSide(MSN_neg,:))',2), [1:length(mean(cell2mat(valorCcg_mazeSide(MSN_neg,:))))]-.1, 300, 'o', 'MarkerEdgeColor', colorScatterMSN, 'MarkerFaceColor', colorCirclScatterMSN)

plot([0 0], [0 6.5], 'color',[0 0 0  0.2],'lineS','-', 'lineW',2, 'HandleVisibility', 'off')

[h, p, ci, stats] = ttest(dataMSN_pos, mean(ccgLagFrLocomotionCellsData{1}(:,MSN_pos)), ...
    'alpha', pValueThreshold/2);
p < pValueThreshold/2;
errorbar(mean(dataMSN_pos), ...
    [1:length(mean(cell2mat(valorCcg_mazeSide(MSN_pos,:))))]+.1, ...
    (ci(2,1) - ci(1,:))/2, ...
    'horizontal', 'k', 'LineStyle', 'none', 'lineW', 1.5, 'HandleVisibility', 'off');

[h, p, ci, stats] = ttest(dataMSN_neg, mean(ccgLagFrLocomotionCellsData{1}(:,MSN_neg)), ...
    'alpha', pValueThreshold/2);
p < pValueThreshold/2;
errorbar(mean(dataMSN_neg), ...
    [1:length(mean(cell2mat(valorCcg_mazeSide(MSN_neg,:))))]-.1, ...
    (ci(2,1) - ci(1,:))/2, ...
    'horizontal', 'k', 'LineStyle', 'none', 'lineW', 1.5, 'HandleVisibility', 'off');

ylim([.8 2.2])
yticks(1:2)
yticklabels({'Left','Right'})
xlabel([])
ylabel('Spatial choice')
set(gca, 'FontWeight', 'bold', 'FontSize', 13);
xlim([-.6, .6])
xticks([])
text(-.25, 1.2, 'I', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

subplot(lines,columns,12)
hold on;
% Positivos
dataMSN_pos = cell2mat(valorCcg_Rew(MSN_pos,:));
scatter(dataMSN_pos, [1:2]+.2, 100, [.682 .698 1], 's', 'filled','MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
plot([mean(ccgLagFrLocomotionCellsData{1}(:,MSN_pos)) mean(ccgLagFrLocomotionCellsData{1}(:,MSN_pos))],[0 6.5], 'color', [.682 .698 1], 'lineS', '--', 'lineW', 2)
scatter(mean(cell2mat(valorCcg_Rew(MSN_pos,:))',2), [1:length(mean(cell2mat(valorCcg_Rew(MSN_pos,:))))]+.1, 300, 's', 'MarkerEdgeColor', colorScatterMSN, 'MarkerFaceColor', colorCirclScatterMSN)
% Negativos
dataMSN_neg = cell2mat(valorCcg_Rew(MSN_neg,:));
scatter(dataMSN_neg, [1:2]-.2, 100, [.682 .698 1], 'o', 'filled','MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
plot([mean(ccgLagFrLocomotionCellsData{1}(:,MSN_neg)) mean(ccgLagFrLocomotionCellsData{1}(:,MSN_neg))],[0 6.5], 'color', [.682 .698 1], 'lineS', '--', 'lineW', 2)
scatter(mean(cell2mat(valorCcg_Rew(MSN_neg,:))',2), [1:length(mean(cell2mat(valorCcg_Rew(MSN_neg,:))))]-.1, 300, 'o', 'MarkerEdgeColor', colorScatterMSN, 'MarkerFaceColor', colorCirclScatterMSN)

plot([0 0], [0 6.5], 'color',[0 0 0  0.2],'lineS','-', 'lineW',2, 'HandleVisibility', 'off')

[h, p, ci, stats] = ttest(dataMSN_pos, mean(ccgLagFrLocomotionCellsData{1}(:,MSN_pos)), ...
    'alpha', pValueThreshold/2);
p < pValueThreshold/2;
errorbar(mean(dataMSN_pos), ...
    [1:length(mean(cell2mat(valorCcg_Rew(MSN_pos,:))))]+.1, ...
    (ci(2,1) - ci(1,:))/2, ...
    'horizontal', 'k', 'LineStyle', 'none', 'lineW', 1.5, 'HandleVisibility', 'off');

[h, p, ci, stats] = ttest(dataMSN_neg, mean(ccgLagFrLocomotionCellsData{1}(:,MSN_neg)), ...
    'alpha', pValueThreshold/2);
p < pValueThreshold/2;
errorbar(mean(dataMSN_neg), ...
    [1:length(mean(cell2mat(valorCcg_Rew(MSN_neg,:))))]-.1, ...
    (ci(2,1) - ci(1,:))/2, ...
    'horizontal', 'k', 'LineStyle', 'none', 'lineW', 1.5, 'HandleVisibility', 'off');

ylim([.8 2.2])
yticks(1:2)
yticklabels({'No\newlineReward','Reward'})
xlabel([])
ylabel('Trial outcome')
set(gca, 'FontWeight', 'bold', 'FontSize', 13);
xlim([-.6, .6])
xticks([-.6:.1:.6])
xlabel('Speed score (r)')
text(-.25, 1.2, 'J', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

f = gcf;
%% ANOVA ONE WAY POSITIVE FSI
clc
dataFSI_pos = cell2mat(valorCcg_portion(FSI_pos,:));
% Suponha que dataFSI já esteja carregado no seu workspace
[nRows, nGroups] = size(dataFSI_pos);

% Transforma a matriz em um vetor (coluna por coluna)
data = dataFSI_pos(:);  % vetor de dados (observações)

% Cria os rótulos de grupo correspondentes
group = repelem(1:nGroups, nRows)';  % grupo 1 para a 1ª coluna, etc.

% Roda ANOVA one-way
[p, tbl, stats] = anova1(data, group);
%% ANOVA ONE WAY NEGATIE FSI
clc
dataFSI_neg = cell2mat(valorCcg_portion(FSI_neg,:));
% Suponha que dataFSI já esteja carregado no seu workspace
[nRows, nGroups] = size(dataFSI_neg);

% Transforma a matriz em um vetor (coluna por coluna)
data = dataFSI_neg(:);  % vetor de dados (observações)

% Cria os rótulos de grupo correspondentes
group = repelem(1:nGroups, nRows)';  % grupo 1 para a 1ª coluna, etc.

% Roda ANOVA one-way
[p, tbl, stats] = anova1(data, group);

%% ANOVA ONE WAY POSITIVE MSN
clc
dataMSN_pos = cell2mat(valorCcg_portion(MSN_pos,:));
% Suponha que dataFSI já esteja carregado no seu workspace
[nRows, nGroups] = size(dataMSN_pos);

% Transforma a matriz em um vetor (coluna por coluna)
data = dataMSN_pos(:);  % vetor de dados (observações)

% Cria os rótulos de grupo correspondentes
group = repelem(1:nGroups, nRows)';  % grupo 1 para a 1ª coluna, etc.

% Roda ANOVA one-way
[p, tbl, stats] = anova1(data, group);
%% ANOVA ONE WAY NEGATIVE MSN
clc
dataMSN_neg = cell2mat(valorCcg_portion(MSN_neg,:));
% Suponha que dataFSI já esteja carregado no seu workspace
[nRows, nGroups] = size(dataMSN_neg);

% Transforma a matriz em um vetor (coluna por coluna)
data = dataMSN_neg(:);  % vetor de dados (observações)

% Cria os rótulos de grupo correspondentes
group = repelem(1:nGroups, nRows)';  % grupo 1 para a 1ª coluna, etc.

% Roda ANOVA one-way
[p, tbl, stats] = anova1(data, group);

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
