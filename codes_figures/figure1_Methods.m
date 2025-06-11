%%
clc, clf, clear
infoSession = 'onlyInTrials';
timeInfoSession = '1sInterval';
% infoSession = 'wholeSession';
folder = strcat('G:\OneDrive\DocData\LFPDataMatlabFormat\spikes\spikes\msnFsiAndLocomotionCellsInfo\allCells\obertosClassification\',timeInfoSession,'\',infoSession,'\');

dataFsiAndLocomotionCells_all = readData(folder);

cd G:\OneDrive\DocData
%%
clc
allCells = {};
strCells = {};
for i=1:length(dataFsiAndLocomotionCells_all)
    allCells = [allCells; dataFsiAndLocomotionCells_all{1,i}.data];
end
% pfcCells = allCells((cell2mat(allCells(:,6)) >=11) & (cell2mat(allCells(:,6)) <=12),:);
strCells = allCells(cell2mat(allCells(:,6)) >=21,:);

%%
clc, clf
fig1 = gcf;
set(gcf,'color','w')
set(fig1,'position',[0, 0, 1400, 1050]);
pos1 = [0.13,0.489517819706499,0.775,0.435482180293501];
linhas = 2;
colunas = 2;

% subplot(linhas,colunas,[1 3])
brainFig = [0, 0.5, 0.6, .48];
axes('Position', brainFig);
% plot(1:1:10,1:1:10)
imagem = imread('C:\Users\phB\Desktop\cerebro1.png'); % Substitua pelo caminho da sua imagem
imshow(imagem); % Plota a imagem no subplot
text(-.15, .9, 'A', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

mazeFig = [.5, 0.5, 0.5, .48];
axes('Position', mazeFig);
% plot(1:1:10,1:1:10)
imagem = imread('C:\Users\phB\Desktop\mazefigura1.png'); % Substitua pelo caminho da sua imagem
imshow(imagem); % Plota a imagem no subplot
text(-.15, .9, 'A', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);


subplot(linhas,colunas,3)
hold on;
[v,bins] = hist(cell2mat(strCells(:,7)), 0:.5:21);
stairs(bins, cumsum(v)/sum(v)*100,'color',[.97, .62, .47],'lineW',3)
xlabel('Firing rate (Hz)');
ylabel({'Neurons (%)'});
xticks([0,1,3,7:7:max(bins)+1])
xlim([0 max(bins)])
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
box off
set(gca, 'XScale', 'log'); % Change x-axis to logarithmic scale
text(-.2, 1.05, 'B', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

hist1_x_pos = [0.13, 0.32, 0.35, 0.16];
ax1_hist_x = axes('Position', hist1_x_pos);
[v,bins] = hist(cell2mat(strCells(cell2mat(strCells(:,6)) == 21,7)), 0:.5:21);
bar(bins, v, 'faceC',[.682 .698 1], 'edgeC',[.682 .698 1])
xlim([0 20])
xticks([0:5:21])
xticklabels([])
ylim([0 21])
yticks([0:7:max(v)])
% xlabel('Firing rate (Hz)');
ylabel({'# Neurons'});
ax = gca;
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
box off
legend('MSN','FontWeight', 'bold','FontSize', 12);
text(-.16, 1.1, 'C', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

hist1_x_pos = [0.13, 0.11, 0.35, 0.16];
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
legend('FSI','FontWeight', 'bold','FontSize', 12);

box off
text(-.16, 1.1, 'D', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

subplot(linhas,colunas,4)
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
theta = atan2d(theta(2), theta(1));  % Convertendo para graus
angle = deg2rad(theta); % Convertendo para radianos
radius = 0.5; % Ajustar conforme necessário
x = radius * cos(angle);
y = radius * sin(angle);
text(x, y, num2str(length(allCells(cell2mat(allCells(:,6)) == 22)')), 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 18);

theta = mean([h(1).Vertices(:,1), h(1).Vertices(:,2)], 1);
theta = atan2d(theta(2), theta(1));  % Convertendo para graus
angle = deg2rad(theta); % Convertendo para radianos
radius = 0.5; % Ajustar conforme necessário
x = radius * cos(angle);
y = radius * sin(angle);
text(x, y, num2str(length(allCells(cell2mat(allCells(:,6)) == 21)')), 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 18);

legend('MSN','FSI','FontWeight', 'bold','FontSize', 12);
set(legend,...
    'Position',[0.85 0.408970547330099 0.0578571422185217 0.0382323723800606],...
    'FontWeight','bold',...
    'FontSize',12);
text(-.16, 1.1, 'E', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 22);


f = gcf;
%%
clc, clf
fig1 = gcf;
set(gcf,'color','w')
set(fig1,'position',[0, 0, 1400, 1050]);
pos1 = [0.13,0.489517819706499,0.775,0.435482180293501];
linhas = 2;
colunas = 3;


hist1_x_pos = [0.07, 0.11, 0.25, 0.36];
ax1_hist_x = axes('Position', hist1_x_pos);
hold on;
% [v,bins] = hist(cell2mat(strCells(:,7)), 0:.5:21);
% stairs(bins, cumsum(v)/sum(v)*100,'color',[.97, .62, .47],'lineW',3)
[f,x] = ecdf(cell2mat(strCells(:,7)));
plot(x,f,'color',[.97, .62, .47],'lineW',3)
xlabel('Firing rate (Hz)');
ylabel({'Cumulative percentage of neurons'});
% xticks([0,.5,1,3,7:7:max(bins)+1])
xticks([0:2:max(bins)])
xlim([0 20])
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
box off
% set(gca, 'XScale', 'log'); % Change x-axis to logarithmic scale
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
% xlabel('Firing rate (Hz)');
ylabel({'# Neurons'});
ax = gca;
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
box off
% legend('MSN','FontWeight', 'bold','FontSize', 12);
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
% legend('FSI','FontWeight', 'bold','FontSize', 12);
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
theta = atan2d(theta(2), theta(1));  % Convertendo para graus
angle = deg2rad(theta); % Convertendo para radianos
radius = 0.5; % Ajustar conforme necessário
x = radius * cos(angle);
y = radius * sin(angle);
text(x, y, num2str(length(allCells(cell2mat(allCells(:,6)) == 22)')), 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 18);

theta = mean([h(1).Vertices(:,1), h(1).Vertices(:,2)], 1);
theta = atan2d(theta(2), theta(1));  % Convertendo para graus
angle = deg2rad(theta); % Convertendo para radianos
radius = 0.5; % Ajustar conforme necessário
x = radius * cos(angle);
y = radius * sin(angle);
text(x, y, num2str(length(allCells(cell2mat(allCells(:,6)) == 21)')), 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 18);

legend('MSN','FSI','FontWeight', 'bold','FontSize', 12);
set(legend,...
    'Position',[0.93 0.42 0.0578571422185217 0.0382323723800606],...
    'FontWeight','bold',...
    'FontSize',12);
text(0, 1.05, 'F', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 22);

f = gcf;
%% SAVE IMAGE
clc
name = 'figure1_Methods_OC_';
saveName = strcat(name, '.tiff');
filename = fullfile('./codes/paperFigures_v1.0/figures/', saveName);
exportgraphics(f, filename, 'Resolution',300)

'done'
%%
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

function [fr] = getFiringRate(data)
   fr = sum(data)/length(data);
end