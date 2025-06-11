%%
clc, clf, clear
infoSession = 'onlyInTrials';
% infoSession = 'wholeSession';

timeInfoSession = '1sInterval';
% timeInfoSession = '0.5sInterval';
% timeInfoSession = '0.25sInterval';
folder = strcat('E:\OneDrive\DocData\LFPDataMatlabFormat\spikes\spikes\msnFsiAndLocomotionCellsInfo\allCells\obertosClassification\',timeInfoSession,'\',infoSession,'\');
dataFsiAndLocomotionCells_all = readData(folder);

cd E:\OneDrive\DocData

%%
clc
allCells = {};
strCells = {};
for i=1:length(dataFsiAndLocomotionCells_all)
    allCells = [allCells; dataFsiAndLocomotionCells_all{1,i}.data];
end
strCells = allCells(cell2mat(allCells(:,6)) >=21,:);

% for i=1:length(allCells)
%     if cell2mat(allCells(i,6)) >= 21
% %         cell2mat(allCells(i,6))
%         strCells = [strCells; allCells(i,:)];
% %         firingRate_all(i,1) = getFiringRate(strCells{i,4}(:,1));
%     end
% end
% % 
load('.\LFPDataMatlabFormat\spikes\spikes\msnFsiAndLocomotionCellsInfo\ccgLagFrLocomotionCellsData_OC_str_surNumber_10000.mat')
valorCcg = ccgLagFrLocomotionCellsData{1};
valoresLag = ccgLagFrLocomotionCellsData{2};
valoresFr = ccgLagFrLocomotionCellsData{3};
locomotioncell = ccgLagFrLocomotionCellsData{4};
%%
clear valorCcg valoresLag valoresFr locomotioncell
clc
maxLag = 10;
surNumber = 10000;
std_gaussiana = .5;

tic
for neuron = 1:length(strCells)
    data = strCells(neuron,:);
    dataFR = data{4}(:,1);
    dataVel = data{5};
    
    tempo_firing_rate = 1:1:length(dataFR);
%     smooth_fr = suavizarDisparos(tempo_firing_rate, dataFR, std_gaussiana);
%     smooth_vel = suavizarDisparos(tempo_firing_rate, dataVel, std_gaussiana);
    smooth_fr = dataFR;
    smooth_vel = dataVel;
    
    for surIdx=1:surNumber
        shiftAmount = randi([20,50],1)/100;
        circularShift = round(length(smooth_vel)*shiftAmount,0);
        [ccg, ~] = xcorr(zscore(circshift(smooth_fr',circularShift)), zscore(smooth_vel), maxLag, 'normalized');
        allSurCCG(:,surIdx) = ccg;
    end
    
    [ccg, lag] = xcorr(zscore(smooth_fr), zscore(smooth_vel), maxLag, 'coeff');
    Threshold(neuron) = max(mean(allSurCCG, 2)+2*std(allSurCCG,1, 2));
    [m i] = max(ccg);
    valorCcg(neuron) = m;
    valoresLag(neuron) = lag(i);
    valoresFr(neuron) = cell2mat(data(7));
    locomotioncell(neuron) =  m > Threshold(neuron);
%     neuron
    m, Threshold(neuron)
    neuron
%     pause
end
toc
%%
% locomotioncellInfo = strCells((locomotioncell == 1),1:3);
locomotioncellInfo = strCells(:,[1:3]);
ccgLagFrLocomotionCellsData = {valorCcg, valoresLag, valoresFr, locomotioncell, locomotioncellInfo};
%%
jsonText = jsonencode(ccgLagFrLocomotionCellsData);
fid = fopen('./dataCCGLagFRLocomotionCellAnimalInfo.json', 'w'); % abrir arquivo para escrita
if fid == -1
    error('Erro ao abrir o arquivo para escrita.');
end

fwrite(fid, jsonText, 'char'); % escrever texto JSON
fclose(fid); % fechar o arquivo
%%
clc
if ~isempty(ccgLagFrLocomotionCellsData)
    saveName = strcat('ccgLagFrLocomotionCellsData_OC_str_surNumber_',num2str(surNumber));
    saveName = strcat('./LFPDataMatlabFormat/spikes/spikes/msnFsiAndLocomotionCellsInfo/',saveName);
    save(saveName, 'ccgLagFrLocomotionCellsData', '-v7.3');
    disp('saved')
end

%%
clc, clf
fig1 = gcf;
set(gcf,'color','w')
set(fig1,'position',[0, 0, 1400, 1050]);
pos1 = [0.13,0.489517819706499,0.775,0.435482180293501];
linhas = 3;
colunas = 9;
colorScatterData = '#a3a3c2';
colorCirclScatterData = '#7575a3';

colorCirclScatterSignificant = '#37a259';
colorScatterSignificant = '#37a259';

% subplot(linhas, colunas,[1 3])
% Create axes
axes1 = axes('Parent',fig1,...
    'Position',[0.13 0.709264705882353 0.235 0.215735294117647]);
% hold(axes1,'on');

h1 = pie([length(valorCcg(locomotioncell==1)),...
    length(valorCcg(locomotioncell==0))],[0 1]);
h1(3).FaceColor = [.639 .639 .761];
h1(3).EdgeColor = [.639 .639 .761];
h1(1).FaceColor = [.216 .635 .349];
h1(1).EdgeColor = [.216 .635 .349];

rotate(h1, [0, 0, 1], 90);

textHandles1 = findobj(gca, 'Type', 'text');
set(textHandles1, 'FontWeight', 'bold','FontSize', 14);
% Ajustando a posição do texto
for k = 1:length(textHandles1)
    % Pega a posição atual
    pos = get(textHandles1(k), 'Position');
    
    % Desloca o texto para mais longe do centro (ajustar os valores conforme necessário)
    pos(1) = pos(1) * 1.2;  % Deslocar ao longo do eixo X
    pos(2) = pos(2) * 1.3;  % Deslocar ao longo do eixo Y
    
    % Atualiza a posição do texto
    set(textHandles1(k), 'Position', pos);
end

legend('Speed cells', 'Non-speed cells','FontWeight', 'bold','FontSize', 12);
set(legend,...
    'Position',[0.24 0.92 0.0857142842454569 0.0382323723800606],...
    'FontWeight','bold',...
    'FontSize',12);

theta = mean([h1(3).Vertices(:,1), h1(3).Vertices(:,2)], 1);
theta = atan2d(theta(2), theta(1));  % Convertendo para graus
angle = deg2rad(theta); % Convertendo para radianos
radius = 0.5; % Ajustar conforme necessário
x = radius * cos(angle);
y = radius * sin(angle);
text(x, y, num2str(length(valorCcg(locomotioncell==0))), 'HorizontalAlignment', 'center', 'FontWeight', 'bold','FontSize', 14);

theta = mean([h1(1).Vertices(:,1), h1(1).Vertices(:,2)], 1);
theta = atan2d(theta(2), theta(1));  % Convertendo para graus
angle = deg2rad(theta); % Convertendo para radianos
radius = 0.5; % Ajustar conforme necessário
x = radius * cos(angle);
y = radius * sin(angle);
text(x, y, num2str(length(valorCcg(locomotioncell==1))), 'HorizontalAlignment', 'center', 'FontWeight', 'bold','FontSize', 14);
text(-.56, 1.1, 'A', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);


% subplot(linhas, colunas,[4 6])
axes2 = axes('Parent',fig1,...
    'Position',[0.41 0.709264705882353 0.232857142857143 0.215735294117647]);
hold(axes2,'on');
hold on
scatter(valoresFr(locomotioncell==0), valorCcg(locomotioncell==0), 'MarkerEdgeColor',colorScatterData,'MarkerFaceColor',colorCirclScatterData)
s = scatter(valoresFr(locomotioncell==1), valorCcg(locomotioncell==1), 'MarkerEdgeColor',colorScatterSignificant,'MarkerFaceColor',colorCirclScatterSignificant);
s.MarkerFaceAlpha = .7;
% s.MarkerEdgeAlpha = .6;

xlabel('Firing rate (Hz)')
ylabel('Speed score (r)')
ylim([0 .51])
xlim([0 20])
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
text(-.2, 1.1, 'B', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);


% subplot(linhas, colunas,[7 9])
axes3 = axes('Parent',fig1,...
    'Position',[0.682857142857143 0.709264705882353 0.222142857142857 0.215735294117647]);
hold(axes3,'on');
hold on;
[values, bins] = hist(valoresLag(locomotioncell==0), -10:1:10);
bar(bins,values,'faceC',colorScatterData,'edgeC',colorScatterData)
[values, bins] = hist(valoresLag(locomotioncell==1), -10:1:10);
bar(bins, values,'faceC',colorCirclScatterSignificant,'edgeC',colorCirclScatterSignificant,'FaceAlpha',.7)%,'FaceAlpha',1,'barW',.5)
% ylim([-25 25])
% yticks([-25:5:25])
% yticklabels([flip(0:5:25),5:5:25])
xlabel('Speed lag (s)')
ylabel('# Neurons')
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
% legend('Data','Significant','FontWeight', 'bold','FontSize', 12);
text(-.18, 1.1, 'C', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);



subplot(linhas, colunas,[10 12])
hold on;
scatter(valoresLag(locomotioncell==0), valoresFr(locomotioncell==0), 'MarkerEdgeColor',colorScatterData,'MarkerFaceColor',colorCirclScatterData)
s = scatter(valoresLag(locomotioncell==1), valoresFr(locomotioncell==1),  'MarkerEdgeColor',colorScatterSignificant,'MarkerFaceColor',colorCirclScatterSignificant);
s.MarkerFaceAlpha = .7;
% s.MarkerEdgeAlpha = .6;
ylabel('Firing rate (Hz)')
xlabel('Speed lag (s)')
xticks([-10:5:10])
set(get(gca, 'YAxis'),'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'),'FontWeight', 'bold','FontSize', 12);
text(-.2, 1.1, 'D', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

subplot(linhas, colunas, [13 15])

h = pie([length(valorCcg(locomotioncell==1 & (cell2mat(strCells(:,6)) == 22)')),...
    length(valorCcg(locomotioncell==0 & (cell2mat(strCells(:,6)) == 22)'))],[0 1]);

h(1).FaceColor = [.988 .263 .4];
h(1).EdgeColor = [.988 .263 .4];
h(3).FaceColor = [.639 .639 .761];
h(3).EdgeColor = [.639 .639 .761];

textHandles1 = findobj(gca, 'Type', 'text');
set(textHandles1, 'FontWeight', 'bold','FontSize', 14);
% Ajustar as posições dos textos das porcentagens

for i = 1:length(textHandles1)
    % Obter as coordenadas do texto atual
    coords = textHandles1(i).Position;
    
    % Aumentar a distância (afastar o texto)
    textHandles1(i).Position = coords * 1.2;  % Aumente o fator (1.1) para mais afastamento
end


lgd = legend('Speed cells','Non-speed cells','FontWeight', 'bold','FontSize', 12);
set(legend,...
    'Position',[0.605 0.42 0.0807142843944686 0.036246275120482],...
    'FontWeight','bold');
title(lgd,'FSI')

rotate(h, [0, 0, 1], 70);

theta = mean([h(3).Vertices(:,1), h(3).Vertices(:,2)], 1);
theta = atan2d(theta(2), theta(1));  % Convertendo para graus
angle = deg2rad(theta); % Convertendo para radianos
radius = 0.5; % Ajustar conforme necessário
x = radius * cos(angle);
y = radius * sin(angle);
text(x, y, num2str(length(valorCcg(locomotioncell==0 & (cell2mat(strCells(:,6)) == 22)'))), 'HorizontalAlignment', 'center', 'FontWeight', 'bold','FontSize', 14);

theta = mean([h(1).Vertices(:,1), h(1).Vertices(:,2)], 1);
theta = atan2d(theta(2), theta(1));  % Convertendo para graus
angle = deg2rad(theta); % Convertendo para radianos
radius = 0.5; % Ajustar conforme necessário
x = radius * cos(angle);
y = radius * sin(angle);
text(x, y, num2str(length(valorCcg(locomotioncell==1 & (cell2mat(strCells(:,6)) == 22)'))), 'HorizontalAlignment', 'center', 'FontWeight', 'bold','FontSize', 14);
text(-.35, 1.1, 'E', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

subplot(linhas, colunas,[16 18])
h = pie([length(valorCcg(locomotioncell==1 & (cell2mat(strCells(:,6)) == 21)')),...
    length(valorCcg(locomotioncell==0 & (cell2mat(strCells(:,6)) == 21)'))],[0 1]);
h(1).FaceColor = [.682 .698 1];
h(1).EdgeColor = [.682 .698 1];
h(3).FaceColor = [.639 .639 .761];
h(3).EdgeColor = [.639 .639 .761];
textHandles1 = findobj(gca, 'Type', 'text');
set(textHandles1, 'FontWeight', 'bold','FontSize', 14);

lgd = legend('Speed cells','Non-speed cells','FontWeight', 'bold','FontSize', 12);
set(legend,...
    'Position',[0.875 0.42 0.0807142843944686 0.036246275120482],...
    'FontWeight','bold');
title(lgd,'MSN')
rotate(h, [0, 0, 1], 90);

for i = 1:length(textHandles1)
    % Obter as coordenadas do texto atual
    coords = textHandles1(i).Position;
    
    % Aumentar a distância (afastar o texto)
    textHandles1(i).Position = coords * 1.2;  % Aumente o fator (1.1) para mais afastamento
end

theta = mean([h(3).Vertices(:,1), h(3).Vertices(:,2)], 1);
theta = atan2d(theta(2), theta(1));  % Convertendo para graus
angle = deg2rad(theta); % Convertendo para radianos
radius = 0.5; % Ajustar conforme necessário
x = radius * cos(angle);
y = radius * sin(angle);
text(x, y, num2str(length(valorCcg(locomotioncell==0 & (cell2mat(strCells(:,6)) == 21)'))), 'HorizontalAlignment', 'center', 'FontWeight', 'bold','FontSize', 14);

theta = mean([h(1).Vertices(:,1), h(1).Vertices(:,2)], 1);
theta = atan2d(theta(2), theta(1));  % Convertendo para graus
angle = deg2rad(theta); % Convertendo para radianos
radius = 0.5; % Ajustar conforme necessário
x = radius * cos(angle);
y = radius * sin(angle);
text(x, y, num2str(length(valorCcg(locomotioncell==1 & (cell2mat(strCells(:,6)) == 21)'))), 'HorizontalAlignment', 'center', 'FontWeight', 'bold','FontSize', 14);
text(-.2, 1.1, 'F', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

subplot(linhas, colunas,[19 21])
hold on;
[values, bins] = hist(valorCcg(locomotioncell==1 & (cell2mat(strCells(:,6)) == 22)'), -.1:.01:.6);
bar(bins, values,'faceC',[.988 .263 .4],'edgeC',[.988 .263 .4])
[values, bins] = hist(valorCcg(locomotioncell==1 & (cell2mat(strCells(:,6)) == 21)'), -.1:.01:.6);
bar(bins,values*-1,'faceC',[.682 .698 1],'edgeC',[.682 .698 1])
xlabel('Speed score (r)')
ylabel('# Neurons')
xlim([0 .5])
yticks(-5:1:5)
ylim([-5, 5])
yticklabels([flip(1:1:5),0:1:5])
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
lgd = legend('FSI','MSN','FontWeight', 'bold','FontSize', 12);
% legend boxoff
title(lgd,'Speed cells')
set(legend,...
    'Position',[0.28 0.28 0.0621428563765117 0.0362462751204823],...
    'FontWeight','bold');
text(-.2, 1.1, 'G', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

% subplot(linhas, colunas,[22])
axes1 = axes('Parent',fig1,...
    'Position',[0.39541095890411 0.11 0.0824461839530333 0.215735294117647]);
hold(axes1,'on');
b = bar([1,10],[mean(valorCcg(locomotioncell==1 & (cell2mat(strCells(:,6)) == 22)')), mean(valorCcg(locomotioncell==1 & (cell2mat(strCells(:,6)) == 21)'))]);
rgbColors = [
    .988 .263 .4;
    .682 .698 1;
];
b.FaceColor = 'flat';
b.EdgeColor = 'flat';
b.CData = rgbColors;

hold on
errorbar([1],mean(valorCcg(locomotioncell==1 & (cell2mat(strCells(:,6)) == 22)')),std(valorCcg(locomotioncell==1 & (cell2mat(strCells(:,6)) == 22)'))/sqrt(length(valorCcg(locomotioncell==1 & (cell2mat(strCells(:,6)) == 22)'))),'k.','markersize',.1,'lineW',1.2)
errorbar([10],mean(valorCcg(locomotioncell==1 & (cell2mat(strCells(:,6)) == 21)')),std(valorCcg(locomotioncell==1 & (cell2mat(strCells(:,6)) == 21)'))/sqrt(length(valorCcg(locomotioncell==1 & (cell2mat(strCells(:,6)) == 21)'))),'k.','markersize',.1,'lineW',1.2)
xticklabels({'FSI','MSN'})
set(gca, 'YAxisLocation', 'right');
% ylabel('Mean speed score (r)')
h = ylabel('\fontsize{13}Mean speed score (r)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [10, 0, 0]); % Ajuste a posição conforme necessário

set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
box off
text(-.1, 1.1, 'H', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);


subplot(linhas, colunas,[24 26])
hold on;
[values, bins] = hist(valoresLag(locomotioncell==1 & (cell2mat(strCells(:,6)) == 22)'), -10:1:10);
bar(bins, values,'faceC',[.988 .263 .4],'edgeC',[.988 .263 .4])
[values, bins] = hist(valoresLag(locomotioncell==1 & (cell2mat(strCells(:,6)) == 21)'), -10:1:10);
bar(bins,values*-1,'faceC',[.682 .698 1],'edgeC',[.682 .698 1])
xlabel('Speed lag (s)')
ylabel('# Neurons')
yticks(-15:5:15)
ylim([-15, 15])
yticklabels([flip(5:5:15),0:5:15])
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
% legend('FSI','MSN','FontWeight', 'bold','FontSize', 12);
text(-.13, 1.1, 'I', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);


% subplot(linhas, colunas,27)
axes1 = axes('Parent',fig1,...
    'Position',[0.837762557077625 0.11 0.0793803000652317 0.215735294117647]);
hold(axes1,'on');
b = bar([1,10],[mean(valoresLag(locomotioncell==1 & (cell2mat(strCells(:,6)) == 22)')), mean(valoresLag(locomotioncell==1 & (cell2mat(strCells(:,6)) == 21)'))]);
rgbColors = [
    .988 .263 .4;
    .682 .698 1;
];
b.FaceColor = 'flat';
b.EdgeColor = 'flat';
b.CData = rgbColors;
hold on
errorbar([1],mean(valoresLag(locomotioncell==1 & (cell2mat(strCells(:,6)) == 22)')),std(valoresLag(locomotioncell==1 & (cell2mat(strCells(:,6)) == 22)'))/sqrt(length(valoresLag(locomotioncell==1 & (cell2mat(strCells(:,6)) == 22)'))),'k.','markersize',.1,'lineW',1.2)
errorbar([10],mean(valoresLag(locomotioncell==1 & (cell2mat(strCells(:,6)) == 21)')),std(valoresLag(locomotioncell==1 & (cell2mat(strCells(:,6)) == 21)'))/sqrt(length(valoresLag(locomotioncell==1 & (cell2mat(strCells(:,6)) == 21)'))),'k.','markersize',.1,'lineW',1.2)
xticklabels({'FSI','MSN'})
set(gca, 'YAxisLocation', 'right');
% ylabel('Mean speed lag (s)')
h = ylabel('\fontsize{13}Mean speed lag (s)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [6, 0, 0]); % Ajuste a posição conforme necessário

set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
box off
text(-.1, 1.1, 'J', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);
ylim([-3 3])

f = gcf;

%% SAVE IMAGE
clc
name = strcat('figure3_ACG_Stats_OC_str_2SD_',timeInfoSession,'_v1');
saveName = strcat(name, '.tiff');
filename = fullfile('G:\OneDrive\DocData\codes\paperFigures_v1.0\figures\', saveName);
% filename = fullfile('./codes/paperFigures_v1.0/figures/', saveName);
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

function taxa_suavizada = suavizarDisparos(vetor_tempo, vetor_disparos, std_gaussiana)
    % std_gaussiana é o desvio padrão da função gaussiana em segundos

    % Calcular o passo de tempo (dt) do vetor de tempo
    dt = vetor_tempo(2) - vetor_tempo(1);

    % Criar vetor de tempo com espaçamento igual ao vetor original, mas que abranja um intervalo maior
    tempo_extenso = min(vetor_tempo):dt:max(vetor_tempo);

    % Criar função gaussiana
    gaussiana = exp(-(tempo_extenso - mean(tempo_extenso)).^2 / (2 * std_gaussiana^2));

    % Normalizar a gaussiana para que a área seja igual a 1
    gaussiana = gaussiana / sum(gaussiana);

    % Criar vetor de taxa suavizada
    taxa_suavizada = conv(vetor_disparos, gaussiana, 'same');

    % Ajustar para o tamanho correto
    taxa_suavizada = taxa_suavizada(1:length(vetor_tempo));
end