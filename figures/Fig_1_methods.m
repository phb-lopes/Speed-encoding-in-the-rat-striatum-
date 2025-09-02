%%
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


%%
clc
allCells = {};
strCells = {};
for i=1:length(dataFsiAndLocomotionCells_all)
    allCells = [allCells; dataFsiAndLocomotionCells_all{1,i}.data];
end

strCells = allCells(cell2mat(allCells(:,6)) >=21,:);

%%
clc, clf
fig1 = gcf;
set(gcf,'color','w')
set(fig1,'position',[0, 0, 1400, 1050]);
pos1 = [0.13,0.489517819706499,0.775,0.435482180293501];

hist1_x_pos = [0.07, 0.11, 0.25, 0.36];
ax1_hist_x = axes('Position', hist1_x_pos);
hold on;
[f,x] = ecdf(cell2mat(strCells(:,7)));
plot(x,f,'color',[.97, .62, .47],'lineW',3)
xlabel('Firing rate (Hz)');
ylabel({'Cumulative percentage of neurons'});
xticks([0:2:max(bins)])
xlim([0 20])
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
box off
text(-.2, 1.05, 'C', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

hist1_x_pos = [0.4, 0.31, 0.3, 0.16];
ax1_hist_x = axes('Position', hist1_x_pos);
[v,bins] = hist(cell2mat(strCells(cell2mat(strCells(:,6)) == 21,7)), 0:.5:21);
bar(bins, v, 'faceC',[.682 .698 1], 'edgeC',[.682 .698 1])
xlim([0 20])
xticks([0:5:21])
xticklabels([])
ylim([0 21])
yticks([0:7:max(v)])
ylabel({'# Neurons'});
ax = gca;
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
box off
legend('MSN', 'FontWeight', 'bold', 'FontSize', 12, 'Position', [0.15, -.05, 1, 1]);
text(-.16, 1.1, 'D', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

hist1_x_pos = [0.4, 0.11, 0.3, 0.16];
ax1_hist_x = axes('Position', hist1_x_pos);
[v,bins] = hist(cell2mat(strCells(cell2mat(strCells(:,6)) == 22,7)), 0:.5:21);
bar(bins, v, 'faceC',[.988 .263 .4], 'edgeC',[.988 .263 .4])
xlim([0 20])
xticks([0:5:21])
ylim([0 21])
yticks([0:7:21])
xlabel('Firing rate (Hz)');
ylabel({'# Neurons'});
ax = gca;
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
legend('FSI', 'FontWeight', 'bold', 'FontSize', 12, 'Position', [0.15, -.25, 1, 1]);
box off
text(-.16, 1.1, 'E', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

hist1_x_pos = [0.67, 0.11, 0.35, 0.36];
ax1_hist_x = axes('Position', hist1_x_pos);
explode = [0 1];
h = pie([length(allCells(cell2mat(allCells(:,6)) == 21)'),...
    length(allCells(cell2mat(allCells(:,6)) == 22)')],[0 1]);
h(1).FaceColor = [.682 .698 1];
h(1).EdgeColor = [.682 .698 1];
h(3).FaceColor = [.988 .263 .4];
h(3).EdgeColor = [.988 .263 .4];
textHandles1 = findobj(gca, 'Type', 'text');
set(textHandles1, 'FontWeight', 'bold','FontSize',14);
theta = mean([h(3).Vertices(:,1), h(3).Vertices(:,2)], 1);
theta = atan2d(theta(2), theta(1));  
angle = deg2rad(theta); 
radius = 0.5; 
x = radius * cos(angle);
y = radius * sin(angle);
text(x, y, num2str(length(allCells(cell2mat(allCells(:,6)) == 22)')), 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 18);

theta = mean([h(1).Vertices(:,1), h(1).Vertices(:,2)], 1);
theta = atan2d(theta(2), theta(1)); 
angle = deg2rad(theta); 
radius = 0.5; 
x = radius * cos(angle);
y = radius * sin(angle);
text(x, y, num2str(length(allCells(cell2mat(allCells(:,6)) == 21)')), 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 18);

legend('MSN','FSI','FontWeight', 'bold','FontSize', 12);
set(legend,...
    'Position',[0.93 0.42 0.0578571422185217 0.0382323723800606],...
    'FontWeight','bold',...
    'FontSize',12);
text(0, 1.05, 'F', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 22);

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