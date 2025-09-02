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

% Convert lag to ms
lag_ms = valoresLag*bin_step;
bin_step = timeInfoSession-overlap;           % Bin size in ms
maxLag_seconds = lag;
maxLag_bins = round(maxLag_seconds/bin_step);  
valoresLag_ms = valoresLag*bin_step;          % Convert bin position to ms
binHisto = -maxLag_seconds:bin_step:maxLag_seconds;

%%
clf, clc
fig1 = gcf;
set(gcf,'color','w')
set(fig1,'position',[0, 0, 1300, 1050]);
pos1 = [0.13,0.489517819706499,0.775,0.435482180293501];
lines = 3;   % Number of rows of subplots
columns = 2;  % Number of columns of subplots

% Define colors
colorScatterData = '#a3a3c2';
colorCirclScatterData = '#7575a3';
colorCirclScatterSignificant = '#37a259';
colorScatterSignificant = '#37a259';

subplot(lines, columns, 1)
h1 = pie([length(valorCcg(speedCell==1)),...
    length(valorCcg(speedCell==0))],[0 1]);
h1(3).FaceColor = [.639 .639 .761];
h1(3).EdgeColor = [.639 .639 .761];
h1(1).FaceColor = [.216 .635 .349];
h1(1).EdgeColor = [.216 .635 .349];

rotate(h1, [0, 0, 1], 90);

textHandles1 = findobj(gca, 'Type', 'text');
set(textHandles1, 'FontWeight', 'bold','FontSize', 14);

% Adjust text position outward from center
for k = 1:length(textHandles1)
    pos = get(textHandles1(k), 'Position');  % Current position
    pos(1) = pos(1) * 1.2;  % Shift text along X
    pos(2) = pos(2) * 1.3;  % Shift text along Y
    set(textHandles1(k), 'Position', pos);   % Update text position
end

legend('Speed cells', 'Non-speed cells','FontWeight', 'bold','FontSize', 12);
set(legend,...
    'Position',[0.34 0.92 0.0857142842454569 0.0382323723800606],...
    'FontWeight','bold',...
    'FontSize',12);

% Add number inside each pie slice
theta = mean([h1(3).Vertices(:,1), h1(3).Vertices(:,2)], 1);
theta = atan2d(theta(2), theta(1));
angle = deg2rad(theta);
radius = 0.5;
x = radius * cos(angle);
y = radius * sin(angle);
text(x, y, num2str(length(valorCcg(speedCell==0))), 'HorizontalAlignment', 'center', 'FontWeight', 'bold','FontSize', 14);

theta = mean([h1(1).Vertices(:,1), h1(1).Vertices(:,2)], 1);
theta = atan2d(theta(2), theta(1));
angle = deg2rad(theta);
x = radius * cos(angle);
y = radius * sin(angle);
text(x, y, num2str(length(valorCcg(speedCell==1))), 'HorizontalAlignment', 'center', 'FontWeight', 'bold','FontSize', 14);

text(-.56, 1.1, 'A', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);


subplot(lines, columns,2)
hold on
s = scatter(valoresFr(speedCell==1), valorCcg(speedCell==1), ...
    'MarkerEdgeColor',colorScatterSignificant,'MarkerFaceColor',colorCirclScatterSignificant);
s.MarkerFaceAlpha = .7;
plot([0, 30], [0 0], Color=[0 0 0 .2], lineStyle='--', lineW=2)

xlabel('Firing rate (Hz)')
ylabel('Speed score (r)')
ylim([-.5 .5])
yticks(-.5:.1:.5)
xlim([0 20])
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
text(-.2, 1.1, 'B', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);


subplot(lines, columns, 3)
hold on;
[values, bins] = hist(valoresLag_ms(speedCell==0),binHisto);
[values, bins] = hist(valoresLag_ms(speedCell==1),binHisto);
bar(bins, values,'faceC',colorCirclScatterSignificant,'edgeC',colorCirclScatterSignificant,'FaceAlpha',.7, 'barW',.35)
xticks(-maxLag_seconds:.5:maxLag_seconds)
xlim([-2.5 2.5])
xlabel('Speed lag (s)')
ylabel('# Neurons')
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
text(-.18, 1.1, 'C', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);


subplot(lines, columns, 4)
hold on;
s = scatter(valoresLag_ms(speedCell==1), valoresFr(speedCell==1), ...
    'MarkerEdgeColor',colorScatterSignificant,'MarkerFaceColor',colorCirclScatterSignificant);
s.MarkerFaceAlpha = .7;
ylabel('Firing rate (Hz)')
xlabel('Speed lag (s)')
xticks(-maxLag_seconds:.5:maxLag_seconds)
xlim([-2.5 2.5])
set(get(gca, 'YAxis'),'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'),'FontWeight', 'bold','FontSize', 12);
text(-.2, 1.1, 'D', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);


subplot(lines, columns, 5)
h = pie([length(valorCcg(speedCell==1 & (cell2mat(strCells(:,6)) == 22)')),...
    length(valorCcg(speedCell==0 & (cell2mat(strCells(:,6)) == 22)'))],[0 1]);

h(1).FaceColor = [.988 .263 .4];
h(1).EdgeColor = [.988 .263 .4];
h(3).FaceColor = [.639 .639 .761];
h(3).EdgeColor = [.639 .639 .761];

textHandles1 = findobj(gca, 'Type', 'text');
set(textHandles1, 'FontWeight', 'bold','FontSize', 14);

% Push labels further from the pie chart
for i = 1:length(textHandles1)
    coords = textHandles1(i).Position;
    textHandles1(i).Position = coords * 1.2;
end

lgd = legend('Speed cells','Non-speed cells','FontWeight', 'bold','FontSize', 12);
set(legend,...
    'Position',[0.4 0.28 0.0807142843944686 0.036246275120482],...
    'FontWeight','bold');
title(lgd,'FSI')

rotate(h, [0, 0, 1], 70);

% Add numbers inside pie slices
theta = mean([h(3).Vertices(:,1), h(3).Vertices(:,2)], 1);
theta = atan2d(theta(2), theta(1));
angle = deg2rad(theta);
radius = 0.5;
x = radius * cos(angle);
y = radius * sin(angle);
text(x, y, num2str(length(valorCcg(speedCell==0 & (cell2mat(strCells(:,6)) == 22)'))), 'HorizontalAlignment', 'center', 'FontWeight', 'bold','FontSize', 14);

theta = mean([h(1).Vertices(:,1), h(1).Vertices(:,2)], 1);
theta = atan2d(theta(2), theta(1));
angle = deg2rad(theta);
x = radius * cos(angle);
y = radius * sin(angle);
text(x, y, num2str(length(valorCcg(speedCell==1 & (cell2mat(strCells(:,6)) == 22)'))), 'HorizontalAlignment', 'center', 'FontWeight', 'bold','FontSize', 14);
text(-.35, 1.1, 'E', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

subplot(lines, columns, 6)
h = pie([length(valorCcg(speedCell==1 & (cell2mat(strCells(:,6)) == 21)')),...
    length(valorCcg(speedCell==0 & (cell2mat(strCells(:,6)) == 21)'))],[0 1]);
h(1).FaceColor = [.682 .698 1];
h(1).EdgeColor = [.682 .698 1];
h(3).FaceColor = [.639 .639 .761];
h(3).EdgeColor = [.639 .639 .761];
textHandles1 = findobj(gca, 'Type', 'text');
set(textHandles1, 'FontWeight', 'bold','FontSize', 14);

lgd = legend('Speed cells','Non-speed cells','FontWeight', 'bold','FontSize', 12);
set(legend,...
    'Position',[0.85 0.28 0.0807142843944686 0.036246275120482],...
    'FontWeight','bold');
title(lgd,'MSN')
rotate(h, [0, 0, 1], 90);

% Push labels outward
for i = 1:length(textHandles1)
    coords = textHandles1(i).Position;
    textHandles1(i).Position = coords * 1.2;
end

% Add numbers inside pie slices
theta = mean([h(3).Vertices(:,1), h(3).Vertices(:,2)], 1);
theta = atan2d(theta(2), theta(1));
angle = deg2rad(theta);
radius = 0.5;
x = radius * cos(angle);
y = radius * sin(angle);
text(x, y, num2str(length(valorCcg(speedCell==0 & (cell2mat(strCells(:,6)) == 21)'))), 'HorizontalAlignment', 'center', 'FontWeight', 'bold','FontSize', 14);

theta = mean([h(1).Vertices(:,1), h(1).Vertices(:,2)], 1);
theta = atan2d(theta(2), theta(1));
angle = deg2rad(theta);
x = radius * cos(angle);
y = radius * sin(angle);
text(x, y, num2str(length(valorCcg(speedCell==1 & (cell2mat(strCells(:,6)) == 21)'))), 'HorizontalAlignment', 'center', 'FontWeight', 'bold','FontSize', 14);
text(-.2, 1.1, 'F', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

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
