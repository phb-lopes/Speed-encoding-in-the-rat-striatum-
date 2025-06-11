%%
clc, clf, clear
cd E:\OneDrive\DocData
folder = 'E:\OneDrive\DocData\LFPDataMatlabFormat\spikes\spikes\msnFsiAndLocomotionCellsInfo\startEndSpksAndVelocity\7sAnalise\spks\';
allSpks = readData(folder);

folder = 'E:\OneDrive\DocData\LFPDataMatlabFormat\spikes\spikes\msnFsiAndLocomotionCellsInfo\startEndSpksAndVelocity\7sAnalise\meanVelocity\start\';
dataMeanVelocityStart = readData(folder);

folder = 'E:\OneDrive\DocData\LFPDataMatlabFormat\spikes\spikes\msnFsiAndLocomotionCellsInfo\startEndSpksAndVelocity\7sAnalise\meanVelocity\end\';
dataMeanVelocityEnd = readData(folder);

% folder = 'G:\OneDrive\DocData\LFPDataMatlabFormat\spikes\spikes\msnFsiAndLocomotionCellsInfo\startEndSpksAndVelocity\7sAnalise\meanVelocity\running\';
% dataMeanVelocityRunning = readData(folder);

infoSession = 'onlyInTrials';
timeInfoSession = '1sInterval';
folder = strcat('E:\OneDrive\DocData\LFPDataMatlabFormat\spikes\spikes\msnFsiAndLocomotionCellsInfo\allCells\obertosClassification\',timeInfoSession,'\',infoSession,'\');
dataFsiAndLocomotionCells_all = readData(folder);

load('.\LFPDataMatlabFormat\spikes\spikes\msnFsiAndLocomotionCellsInfo\ccgLagFrLocomotionCellsData_OC_str_surNumber_10000.mat')
% valorCcg = ccgLagFrLocomotionCellsData{1};
% valoresLag = ccgLagFrLocomotionCellsData{2};
% valoresFr = ccgLagFrLocomotionCellsData{3};
locomotioncell = ccgLagFrLocomotionCellsData{4};

%%
clc
allMeanVelocityStart = {};
allMeanVelocityEnd = {};
allMeanVelocityRunning = {};

for i=1:length(dataMeanVelocityStart)
    allMeanVelocityStart(i) = {dataMeanVelocityStart{i}.mediaStartMovment};
    allMeanVelocityEnd(i) = {dataMeanVelocityEnd{i}.mediaEndMovment};
%     allMeanVelocityRunning(i) = {dataMeanVelocityRunning{i}.mediaRunningMovment};
end

somaS = [];
somaE = [];
% somaR = zeros(size(dataMeanVelocityRunning{1}.mediaRunningMovment)); 
% somaR = [];
for i=1:length(allMeanVelocityStart)
    somaS = [somaS; {allMeanVelocityStart{i}}];
    somaE = [somaE; {allMeanVelocityEnd{i}}];
%     somaR = [somaR; {allMeanVelocityStart{i}}];
%     somaR = somaR + allMeanVelocityRunning{i};
end

% mediaFinalS = somaS / i;
% mediaFinalE = somaE / i; 
% mediaFinalR = somaR / i; 

%%
matrizS = cell2mat(somaS');
mediaFinalS = mean(matrizS, 2);

matrizE = cell2mat(somaE');
mediaFinalE = mean(matrizE, 2);

%%
clc
allCells = {};
strCells = {};
for i=1:length(dataFsiAndLocomotionCells_all)
    allCells = [allCells; dataFsiAndLocomotionCells_all{1,i}.data];
end
% allCells = allCells(~strcmp(allCells(:,1), 'R268 20140930'), :);
strCells = allCells(cell2mat(allCells(:,6)) >=21,:);
strCells = strCells(locomotioncell,:);

spkNeurons = {};
for j=1:length(allSpks)
    spkNeurons = [spkNeurons; allSpks{1,j}.spks];
end

indicesLC_MSN = [];
indicesLC_FSI = [];
for i = 1:size(strCells, 1)
    for j = 1:size(spkNeurons, 1)
        if isequal(strCells(i, [1:3]), spkNeurons(j, [1,2,3]))
            strCells(i, [1:3,6])
            if cell2mat(strCells(i,6)) == 21
                indicesLC_MSN = [indicesLC_MSN; j]; % Adiciona o índice
            else
                indicesLC_FSI = [indicesLC_FSI; j]; % Adiciona o índice
            end
            break; % Sai do loop se encontrar uma correspondência
        end
    end
end

indicesNoLC_MSN = [];
indicesNoLC_FSI = [];
% indicesNoLC = setdiff(1:length(spkNeurons), vertcat(indicesLC_MSN,indicesLC_FSI));

strCells = allCells(cell2mat(allCells(:,6)) >=21,:);
strCells = strCells(~locomotioncell,:);
for i = 1:size(strCells, 1)
    for j = 1:size(spkNeurons, 1)
        if isequal(strCells(i, [1:3]), spkNeurons(j, [1,2,3]))
            if cell2mat(strCells(i,6)) == 21
                indicesNoLC_MSN = [indicesNoLC_MSN; j]; % Adiciona o índice
            else
                indicesNoLC_FSI = [indicesNoLC_FSI; j]; % Adiciona o índice
            end
            break; % Sai do loop se encontrar uma correspondência
        end
    end
end

%%
clc
strCells = allCells(cell2mat(allCells(:,6)) >=21,:);
valorCcg = ccgLagFrLocomotionCellsData{1}';
[m, I] = maxk(valorCcg,10);
I
strCells(12,[1:3,6])
%%
clc

binSize = .1;
edges7s = 0:binSize:7;


spksLcS_MSN = [];
allN_spksLcS_MSN = [];
spksLcS_FSI = [];
allN_spksLcS_FSI = [];

spksLcE_MSN = [];
allN_spksLcE_MSN = [];
spksLcE_FSI = [];
allN_spksLcE_FSI = [];

aux_exempleLC_MSN_S = [];
aux_exempleLC_FSI_S = [];
aux_exempleLC_MSN_E = [];
aux_exempleLC_FSI_E = [];
idxLC_MSN = 25;%25
c = 0;
for i=1:length(indicesLC_MSN)

    aux = spkNeurons(indicesLC_MSN(i),5);
    for trial=1:length(aux{1,1})
        spikeCounts = histcounts(cell2mat(aux{1,1}(trial)), edges7s);
        firingRate = spikeCounts / binSize;
        spksLcS_MSN = [spksLcS_MSN;firingRate];
    end
%     allN_spksLcS_MSN = [allN_spksLcS_MSN; normalize(mean(spksLcS_MSN),'range')];
%     allN_spksLcS_MSN = [allN_spksLcS_MSN; mean(zscore(spksLcS_MSN))];
    allN_spksLcS_MSN = [allN_spksLcS_MSN; mean(spksLcS_MSN)/mean(mean(spksLcS_MSN))];
    
    if indicesLC_MSN(i) == idxLC_MSN;%25   % Bom
        for trial=1:length(aux{1,1})
            spikeCounts = histcounts(cell2mat(aux{1,1}(trial)), edges7s);
            firingRate = spikeCounts / binSize;
            aux_exempleLC_MSN_S = [aux_exempleLC_MSN_S;firingRate];
        end
        exempleLC_MSN_S = normalize(mean(aux_exempleLC_MSN_S),'range');
    end

    aux = spkNeurons(indicesLC_MSN(i),7);
    for trial=1:length(aux{1,1})
        spikeCounts = histcounts(cell2mat(aux{1,1}(trial)), edges7s);
        firingRate = spikeCounts / binSize;
        spksLcE_MSN = [spksLcE_MSN;firingRate];
    end
%     allN_spksLcE_MSN = [allN_spksLcE_MSN; normalize(mean(spksLcE_MSN),'range')];
%     allN_spksLcE_MSN = [allN_spksLcE_MSN; mean(zscore(spksLcE_MSN))];
    allN_spksLcE_MSN = [allN_spksLcE_MSN; mean(spksLcE_MSN)/mean(mean(spksLcE_MSN))];
     
    if indicesLC_MSN(i) == idxLC_MSN;%25   % Bom
        for trial=1:length(aux{1,1})
            spikeCounts = histcounts(cell2mat(aux{1,1}(trial)), edges7s);
            firingRate = spikeCounts / binSize;
            aux_exempleLC_MSN_E = [aux_exempleLC_MSN_E;firingRate];
        end
        exempleLC_MSN_E = normalize(mean(aux_exempleLC_MSN_E),'range');
    end    
end

idxLC_FSI = 28;%28
for i=1:length(indicesLC_FSI)
 
    aux = spkNeurons(indicesLC_FSI(i),5);
    for trial=1:length(aux{1,1})
        spikeCounts = histcounts(cell2mat(aux{1,1}(trial)), edges7s);
        firingRate = spikeCounts / binSize;
        spksLcS_FSI = [spksLcS_FSI;firingRate];
    end
%     allN_spksLcS_FSI = [allN_spksLcS_FSI; normalize(mean(spksLcS_FSI),'range')];
%     allN_spksLcS_FSI = [allN_spksLcS_FSI; mean(zscore(spksLcS_FSI))];
    allN_spksLcS_FSI = [allN_spksLcS_FSI; mean(spksLcS_FSI)/mean(mean(spksLcS_FSI))];
    
    if indicesLC_FSI(i) == idxLC_FSI;%28
        for trial=1:length(aux{1,1})
            spikeCounts = histcounts(cell2mat(aux{1,1}(trial)), edges7s);
            firingRate = spikeCounts / binSize;
            aux_exempleLC_FSI_S = [aux_exempleLC_FSI_S;firingRate];
        end
        exempleLC_FSI_S = normalize(mean(aux_exempleLC_FSI_S),'range');
    end
    
    aux = spkNeurons(indicesLC_FSI(i),7);
    for trial=1:length(aux{1,1})
        spikeCounts = histcounts(cell2mat(aux{1,1}(trial)), edges7s);
        firingRate = spikeCounts / binSize;
        spksLcE_FSI = [spksLcE_FSI;firingRate];
    end    
%     allN_spksLcE_FSI = [allN_spksLcE_FSI; normalize(mean(spksLcE_FSI),'range')];
%     allN_spksLcE_FSI = [allN_spksLcE_FSI; mean(zscore(spksLcE_FSI))];
    allN_spksLcE_FSI = [allN_spksLcE_FSI; mean(spksLcE_FSI)/mean(mean(spksLcE_FSI))];

    if indicesLC_FSI(i) == idxLC_FSI;%28
        for trial=1:length(aux{1,1})
            spikeCounts = histcounts(cell2mat(aux{1,1}(trial)), edges7s);
            firingRate = spikeCounts / binSize;
            aux_exempleLC_FSI_E = [aux_exempleLC_FSI_E;firingRate];
        end
        exempleLC_FSI_E = normalize(mean(aux_exempleLC_FSI_E),'range');
    end
end

spksNoLcS_MSN = [];
allN_spksNoLcS_MSN = [];

spksNoLcS_FSI = [];
allN_spksNoLcS_FSI = [];

spksNoLcE_MSN = [];
allN_spksNoLcE_MSN = [];

spksNoLcE_FSI = [];
allN_spksNoLcE_FSI = [];

aux_exempleNoLCMSN_S = [];
aux_exempleNoLCFSI_S = [];
aux_exempleNoLCMSN_E = [];
aux_exempleNoLCFSI_E = [];

for i=1:length(indicesNoLC_MSN)
 
    aux = spkNeurons(indicesNoLC_MSN(i),5);
    for trial=1:length(aux{1,1})
        spikeCounts = histcounts(cell2mat(aux{1,1}(trial)), edges7s);
        firingRate = spikeCounts / binSize;
        spksNoLcS_MSN = [spksNoLcS_MSN;firingRate];
    end
%     allN_spksNoLcS_MSN = [allN_spksNoLcS_MSN; normalize(mean(spksNoLcS_MSN),'range')];
%     allN_spksNoLcS_MSN = [allN_spksNoLcS_MSN; mean(zscore(spksNoLcS_MSN))];
    allN_spksNoLcS_MSN = [allN_spksNoLcS_MSN; mean(spksNoLcS_MSN)/mean(mean(spksNoLcS_MSN))];
    
    if indicesNoLC_MSN(i) == 22
        for trial=1:length(aux{1,1})
            spikeCounts = histcounts(cell2mat(aux{1,1}(trial)), edges7s);
            firingRate = spikeCounts / binSize;
            aux_exempleNoLCMSN_S = [aux_exempleNoLCMSN_S;firingRate];
        end
        exempleNoLCMSN_S = normalize(mean(aux_exempleNoLCMSN_S),'range');
    end

    aux = spkNeurons(indicesNoLC_MSN(i),7);
    for trial=1:length(aux{1,1})
        spikeCounts = histcounts(cell2mat(aux{1,1}(trial)), edges7s);
        firingRate = spikeCounts / binSize;
        spksNoLcE_MSN = [spksNoLcE_MSN;firingRate];
    end
%     allN_spksNoLcE_MSN = [allN_spksNoLcE_MSN; normalize(mean(spksNoLcE_MSN),'range')];
%     allN_spksNoLcE_MSN = [allN_spksNoLcE_MSN; mean(zscore(spksNoLcE_MSN))];
    allN_spksNoLcE_MSN = [allN_spksNoLcE_MSN; mean(spksNoLcE_MSN)/mean(mean(spksNoLcE_MSN))];
    
    if indicesNoLC_MSN(i) == 22
        for trial=1:length(aux{1,1})
            spikeCounts = histcounts(cell2mat(aux{1,1}(trial)), edges7s);
            firingRate = spikeCounts / binSize;
            aux_exempleNoLCMSN_E = [aux_exempleNoLCMSN_E;firingRate];
        end
        exempleNoLCMSN_E = normalize(mean(aux_exempleNoLCMSN_E),'range');
    end
end

idxNoLC_FSI = 23;%23
for i=1:length(indicesNoLC_FSI)
    
    aux = spkNeurons(indicesNoLC_FSI(i),5);
    for trial=1:length(aux{1,1})
        spikeCounts = histcounts(cell2mat(aux{1,1}(trial)), edges7s);
        firingRate = spikeCounts / binSize;
        spksNoLcS_FSI = [spksNoLcS_FSI;firingRate];
    end
%     allN_spksNoLcS_FSI = [allN_spksNoLcS_FSI; normalize(mean(spksNoLcS_FSI),'range')];
%     allN_spksNoLcS_FSI = [allN_spksNoLcS_FSI; mean(zscore(spksNoLcS_FSI))];
    allN_spksNoLcS_FSI = [allN_spksNoLcS_FSI; mean(spksNoLcS_FSI)/mean(mean(spksNoLcS_FSI))];
    
    if indicesNoLC_FSI(i) == idxNoLC_FSI;%23 % Bom
        for trial=1:length(aux{1,1})
            spikeCounts = histcounts(cell2mat(aux{1,1}(trial)), edges7s);
            firingRate = spikeCounts / binSize;
            aux_exempleNoLCFSI_S = [aux_exempleNoLCFSI_S;firingRate];
        end
         exempleNoLCFSI_S = normalize(mean(aux_exempleNoLCFSI_S),'range');
    end
    
    aux = spkNeurons(indicesNoLC_FSI(i),7);
    for trial=1:length(aux{1,1})
        spikeCounts = histcounts(cell2mat(aux{1,1}(trial)), edges7s);
        firingRate = spikeCounts / binSize;
        spksNoLcE_FSI = [spksNoLcE_FSI;firingRate];
    end
%     allN_spksNoLcE_FSI = [allN_spksNoLcE_FSI; normalize(mean(spksNoLcE_FSI),'range')];
%     allN_spksNoLcE_FSI = [allN_spksNoLcE_FSI; mean(zscore(spksNoLcE_FSI))];
    allN_spksNoLcE_FSI = [allN_spksNoLcE_FSI; mean(spksNoLcE_FSI)/mean(mean(spksNoLcE_FSI))];
    
    if indicesNoLC_FSI(i) == idxNoLC_FSI;%23 % Bom
        for trial=1:length(aux{1,1})
            spikeCounts = histcounts(cell2mat(aux{1,1}(trial)), edges7s);
            firingRate = spikeCounts / binSize;
            aux_exempleNoLCFSI_E = [aux_exempleNoLCFSI_E;firingRate];
        end
         exempleNoLCFSI_E = normalize(mean(aux_exempleNoLCFSI_S),'range');
    end
end
%%
clf, clc
fig1 = gcf;
linhas = 2;
colunas = 2;
set(gcf,'color','w')
set(fig1,'position',[0, 0, 1400, 1050]);
pos1 = [0.13,0.489517819706499,0.775,0.435482180293501];
yMaxSpks = 1.6;
yMaxVel = 50;
yticksVel = [0:5:yMaxVel];
xLim = [0 7];

colorAllNeurons = '#37a259';

figWidth = 0.35;   % largura de cada gráfico
figHeight = 0.4;   % altura de cada gráfico
hSpacing = 0.12;   % espaço horizontal entre colunas
vSpacing = 0.08;   % espaço vertical entre linhas

x1 = 0.1;
y1 = 0.55;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(linhas, colunas,1)

ax1 = axes('Position', [x1, y1, figWidth, figHeight]);
% nexttile
title('Speed cells', 'FontWeight', 'bold','FontSize', 18,'color','k');
hold on;
yyaxis left;
plot(edges7s(1:end-1),smooth(mean(vertcat(allN_spksLcS_FSI,allN_spksLcS_MSN)),10),'color',colorAllNeurons,'lineW',2, 'lineS','-')
plot(edges7s(1:end-1),smooth(mean(allN_spksLcS_FSI),10),'color',[.988 .263 .4],'lineW',2, 'lineS','-')
plot(edges7s(1:end-1),smooth(mean(allN_spksLcS_MSN),10),'color',[.682 .698 1],'lineW',2, 'lineS','-')

A = smooth(mean(vertcat(allN_spksLcS_FSI,allN_spksLcS_MSN)),10)+...
    std(mean(vertcat(allN_spksLcS_FSI,allN_spksLcS_MSN)),1)...
    /sqrt(size(vertcat(allN_spksLcS_FSI,allN_spksLcS_MSN),1));
B = smooth(mean(vertcat(allN_spksLcS_FSI,allN_spksLcS_MSN)),10)-...
    std(mean(vertcat(allN_spksLcS_FSI,allN_spksLcS_MSN)),1)...
    /sqrt(size(vertcat(allN_spksLcS_FSI,allN_spksLcS_MSN),1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.216 .635 .349],'facealpha',0.5,'EdgeColor','none', 'HandleVisibility', 'off')
 

A = smooth(mean(allN_spksLcS_FSI),10)+...
    std(mean(allN_spksLcS_FSI),1)...
    /sqrt(size(allN_spksLcS_FSI,1));
B = smooth(mean(allN_spksLcS_FSI),10)-...
    std(mean(allN_spksLcS_FSI),1)...
    /sqrt(size(allN_spksLcS_FSI,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.988 .263 .4],'facealpha',0.2,'EdgeColor','none', 'HandleVisibility', 'off')
 
A = smooth(mean(allN_spksLcS_MSN),10)+...
    std(mean(allN_spksLcS_MSN),1)...
    /sqrt(size(allN_spksLcS_MSN,1));
B = smooth(mean(allN_spksLcS_MSN),10)-...
    std(mean(allN_spksLcS_MSN),1)...
    /sqrt(size(allN_spksLcS_MSN,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.682 .698 1],'facealpha',0.2,'EdgeColor','none', 'HandleVisibility', 'off')
 

xlim(xLim)
ylim([0.6 yMaxSpks])
xticks(0:1:7)
xticklabels([])
yticks(0:.2:1.6)
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
ylabel({'\fontsize{18}Locomotion onset', '\fontsize{13}Normalized firing rate'}, 'FontWeight', 'bold');
% ylabel({'\fontsize{18}Locomotion onset', '\fontsize{13}Normalized firing rate - \color[rgb]{0.216, 0.635, 0.349}All\color{black},\color[rgb]{.988 .263 .4}FSI\color{black},\color[rgb]{.682 .698 1}MSN'}, 'FontWeight', 'bold');

yyaxis right;
x = linspace(0, 7, length(mediaFinalS));
plot(x, smooth(mediaFinalS),'k','lineW',2); 
ylim([0 35]);
h = ylabel('\fontsize{13}Mean speed (cm/s)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [0.4, 0, 0]); % Ajuste a posição conforme necessário
xt = xticks; % Obter os valores dos xticks
for i = 2:length(xt)-1
    xline(xt(i), '--', 'color',[.8 .8 .8 .2], 'HandleVisibility', 'off'); % Adiciona uma linha vertical em cada xtick
end

A = smooth(mean(matrizS, 2),10)+...
    std(mean(matrizS, 2),1)...
    /sqrt(size(somaS,1));
B = smooth(mean(matrizS, 2),10)-...
    std(mean(matrizS, 2),1)...
    /sqrt(size(somaS,1));
fill([x fliplr(x)], ...
     [A' fliplr(B')],...
     [.51 .51 .51],'facealpha',0.2,'EdgeColor','none', 'HandleVisibility', 'off')
% lgd = legend({'All (n = 102)','FSI (n = 44)','MSN (n = 58)','Mean speed'});
lgd = legend({'All','FSI','MSN','Speed'});
set(legend,...
    'Position',[0.072 0.89 0.16428571156093 0.0526315775268009],...
    'NumColumns',1,...
    'FontWeight','bold',...
    'FontSize',13);
legend boxoff
% title(lgd,'Speed cells')

text(-.15, 1.05, 'A', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(linhas, colunas,2)
% nexttile

ax2 = axes('Position', [x1 + figWidth + hSpacing, y1, figWidth, figHeight]);

hold on
plot(edges7s(1:end-1),smooth(mean(vertcat(allN_spksNoLcS_FSI,allN_spksNoLcS_MSN)),10),'color',colorAllNeurons,'lineW',2, 'lineS',':')
% plot(edges7s(1:end-1)-10,smooth(mean(allN_spksNoLcS_FSI),10),'color',[.988 .263 .4],'lineW',2, 'lineS',':')
% plot(edges7s(1:end-1)-10,smooth(mean(allN_spksNoLcS_MSN),10),'color',[.682 .698 1],'lineW',2, 'lineS',':')
% plot(edges7s(1:end-1)-10,smooth(mean(vertcat(allN_spksLcS_FSI,allN_spksLcS_MSN)),10),'color',colorAllNeurons,'lineW',2, 'lineS','-')
plot(edges7s(1:end-1),smooth(mean(allN_spksNoLcS_FSI),10),'color',[.988 .263 .4],'lineW',2, 'lineS',':')
plot(edges7s(1:end-1),smooth(mean(allN_spksNoLcS_MSN),10),'color',[.682 .698 1],'lineW',2, 'lineS',':')
% plot(edges7s(1:end-1)-10,smooth(mean(allN_spksLcS_FSI),10),'color',[.988 .263 .4],'lineW',2, 'lineS','-')
% plot(edges7s(1:end-1)-10,smooth(mean(allN_spksLcS_MSN),10),'color',[.682 .698 1],'lineW',2, 'lineS','-')

A = smooth(mean(allN_spksNoLcS_FSI),10)+...
    std(mean(allN_spksNoLcS_FSI),1)...
    /sqrt(size(allN_spksNoLcS_FSI,1));
B = smooth(mean(allN_spksNoLcS_FSI),10)-...
    std(mean(allN_spksNoLcS_FSI),1)...
    /sqrt(size(allN_spksNoLcS_FSI,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.988 .263 .4],'facealpha',0.2,'EdgeColor','none', 'HandleVisibility', 'off')
 
A = smooth(mean(allN_spksNoLcS_MSN),10)+...
    std(mean(allN_spksNoLcS_MSN),1)...
    /sqrt(size(allN_spksNoLcS_MSN,1));
B = smooth(mean(allN_spksNoLcS_MSN),10)-...
    std(mean(allN_spksNoLcS_MSN),1)...
    /sqrt(size(allN_spksNoLcS_MSN,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.682 .698 1],'facealpha',0.2,'EdgeColor','none', 'HandleVisibility', 'off')
 
A = smooth(mean(vertcat(allN_spksNoLcS_FSI,allN_spksNoLcS_MSN)),10)+...
    std(mean(vertcat(allN_spksNoLcS_FSI,allN_spksNoLcS_MSN)),1)...
    /sqrt(size(vertcat(allN_spksNoLcS_FSI,allN_spksNoLcS_MSN),1));
B = smooth(mean(vertcat(allN_spksNoLcS_FSI,allN_spksNoLcS_MSN)),10)-...
    std(mean(vertcat(allN_spksNoLcS_FSI,allN_spksNoLcS_MSN)),1)...
    /sqrt(size(vertcat(allN_spksNoLcS_FSI,allN_spksNoLcS_MSN),1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.216 .635 .349],'facealpha',0.5,'EdgeColor','none', 'HandleVisibility', 'off')


title('Non-speed cells', 'FontWeight', 'bold','FontSize', 18,'color','k');
hold on;
yyaxis left;

xlim(xLim)
ylim([0.6 yMaxSpks])
xticks(0:1:7)
xticklabels([])
yticks(0:.2:1.6)
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
ylabel('Normalized firing rate', 'FontWeight', 'bold','FontSize', 13);
% ylabel({'\fontsize{13}Normalized firing rate - \color[rgb]{0.216, 0.635, 0.349}All\color{black},\color[rgb]{.988 .263 .4}FSI\color{black},\color[rgb]{.682 .698 1}MSN'}, 'FontWeight', 'bold');

yyaxis right;
x = linspace(0, 7, length(mediaFinalS));
plot(x, smooth(mediaFinalS),'k','lineW',2); 
ylim([0 35]);
h = ylabel('\fontsize{13}Mean speed (cm/s)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [0.4, 0, 0]); % Ajuste a posição conforme necessário
xt = xticks; % Obter os valores dos xticks
for i = 2:length(xt)-1
    xline(xt(i), '--', 'color',[.8 .8 .8 .2], 'HandleVisibility', 'off'); % Adiciona uma linha vertical em cada xtick
end

A = smooth(mean(matrizS, 2),10)+...
    std(mean(matrizS, 2),1)...
    /sqrt(size(somaS,1));
B = smooth(mean(matrizS, 2),10)-...
    std(mean(matrizS, 2),1)...
    /sqrt(size(somaS,1));
fill([x fliplr(x)], ...
     [A' fliplr(B')],...
     [.51 .51 .51],'facealpha',0.2,'EdgeColor','none', 'HandleVisibility', 'off')

% size(somaE,1)
% lgd = legend({'All (n = 23)','FSI (n = 7)','MSN (n = 16)','Mean speed'});
% set(legend,...
%     'Position',[0.58 0.865 0.16428571156093 0.0526315775268009],...
%     'NumColumns',1,...
%     'FontWeight','bold',...
%     'FontSize',13);
% title(lgd,'Non-speed cells')

text(-.16, 1.05, 'C', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(linhas, colunas,3)
% nexttile

ax3 = axes('Position', [x1, y1 - figHeight - vSpacing, figWidth, figHeight]);
hold on;
yyaxis left;
plot(edges7s(1:end-1),smooth(mean(allN_spksLcE_FSI),10),'color',[.988 .263 .4],'lineW',2, 'lineS','-');
plot(edges7s(1:end-1),smooth(mean(allN_spksLcE_MSN),10),'color',[.682 .698 1],'lineW',2, 'lineS','-');

plot(edges7s(1:end-1),smooth(mean(vertcat(allN_spksLcE_FSI,allN_spksLcE_MSN)),10),'color',colorAllNeurons,'lineW',2, 'lineS','-')

A = smooth(mean(allN_spksLcE_FSI),10)+...
    std(mean(allN_spksLcE_FSI),1)...
    /sqrt(size(allN_spksLcE_FSI,1));
B = smooth(mean(allN_spksLcE_FSI),10)-...
    std(mean(allN_spksLcE_FSI),1)...
    /sqrt(size(allN_spksLcE_FSI,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.988 .263 .4],'facealpha',0.2,'EdgeColor','none')
 
A = smooth(mean(allN_spksLcE_MSN),10)+...
    std(mean(allN_spksLcE_MSN),1)...
    /sqrt(size(allN_spksLcE_MSN,1));
B = smooth(mean(allN_spksLcE_MSN),10)-...
    std(mean(allN_spksLcE_MSN),1)...
    /sqrt(size(allN_spksLcE_MSN,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.682 .698 1],'facealpha',0.2,'EdgeColor','none')
 
 
 A = smooth(mean(vertcat(allN_spksLcE_FSI,allN_spksLcE_MSN)),10)+...
    std(mean(vertcat(allN_spksLcE_FSI,allN_spksLcE_MSN)),1)...
    /sqrt(size(vertcat(allN_spksLcE_FSI,allN_spksLcE_MSN),1));
B = smooth(mean(vertcat(allN_spksLcE_FSI,allN_spksLcE_MSN)),10)-...
    std(mean(vertcat(allN_spksLcE_FSI,allN_spksLcE_MSN)),1)...
    /sqrt(size(vertcat(allN_spksLcE_FSI,allN_spksLcE_MSN),1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.216 .635 .349],'facealpha',0.5,'EdgeColor','none')
 
xlim(xLim)
ylim([0.6 yMaxSpks])
xticks(0:1:7)
yticks(0:.2:1.6)
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
ylabel({'\fontsize{18}Locomotion offset', '\fontsize{13}Normalized firing rate'}, 'FontWeight', 'bold');
% ylabel({'\fontsize{18}Locomotion offset', '\fontsize{13}Normalized firing rate - \color[rgb]{0.216, 0.635, 0.349}All\color{black},\color[rgb]{.988 .263 .4}FSI\color{black},\color[rgb]{.682 .698 1}MSN'}, 'FontWeight', 'bold');
yyaxis right;
x = linspace(0, 7, length(mediaFinalE));
plot(x, smooth(mediaFinalE),'k','lineW',2);
ylim([0 35])
h = ylabel('\fontsize{13}Mean speed (cm/s)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [0.4, 0, 0]); % Ajuste a posição conforme necessário
xlabel('Time (s)','FontWeight', 'bold','FontSize', 12);

xt = xticks; % Obter os valores dos xticks
for i = 2:length(xt)-1
    xline(xt(i), '--', 'color',[.8 .8 .8 .2]); % Adiciona uma linha vertical em cada xtick
end
A = smooth(mean(matrizE, 2),10)+...
    std(mean(matrizE, 2),1)...
    /sqrt(size(somaE,1));
B = smooth(mean(matrizE, 2),10)-...
    std(mean(matrizE, 2),1)...
    /sqrt(size(somaE,1));
fill([x fliplr(x)], ...
     [A' fliplr(B')],...
     [.51 .51 .51],'facealpha',0.2,'EdgeColor','none', 'HandleVisibility', 'off')
 
text(-.15, 1.05, 'B', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(linhas, colunas,4)
% nexttile
ax4 = axes('Position', [x1 + figWidth + hSpacing, y1 - figHeight - vSpacing, figWidth, figHeight]);

hold on;
yyaxis left;
plot(edges7s(1:end-1),smooth(mean(allN_spksNoLcE_FSI),10),'color',[.988 .263 .4],'lineW',2, 'lineS',':');
plot(edges7s(1:end-1),smooth(mean(allN_spksNoLcE_MSN),10),'color',[.682 .698 1],'lineW',2, 'lineS',':');

plot(edges7s(1:end-1),smooth(mean(vertcat(allN_spksNoLcE_FSI,allN_spksNoLcE_MSN)),10),'color',[.216 .635 .349],'lineW',2, 'lineS',':')

A = smooth(mean(allN_spksNoLcE_FSI),10)+...
    std(mean(allN_spksNoLcE_FSI),1)...
    /sqrt(size(allN_spksNoLcE_FSI,1));
B = smooth(mean(allN_spksNoLcE_FSI),10)-...
    std(mean(allN_spksNoLcE_FSI),1)...
    /sqrt(size(allN_spksNoLcE_FSI,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.988 .263 .4],'facealpha',0.2,'EdgeColor','none')
 
A = smooth(mean(allN_spksNoLcE_MSN),10)+...
    std(mean(allN_spksNoLcE_MSN),1)...
    /sqrt(size(allN_spksNoLcE_MSN,1));
B = smooth(mean(allN_spksNoLcE_MSN),10)-...
    std(mean(allN_spksNoLcE_MSN),1)...
    /sqrt(size(allN_spksNoLcE_MSN,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.682 .698 1],'facealpha',0.2,'EdgeColor','none')

A = smooth(mean(vertcat(allN_spksNoLcE_FSI,allN_spksNoLcE_MSN)),10)+...
    std(mean(vertcat(allN_spksNoLcE_FSI,allN_spksNoLcE_MSN)),1)...
    /sqrt(size(vertcat(allN_spksNoLcE_FSI,allN_spksNoLcE_MSN),1));
B = smooth(mean(vertcat(allN_spksNoLcE_FSI,allN_spksNoLcE_MSN)),10)-...
    std(mean(vertcat(allN_spksNoLcE_FSI,allN_spksNoLcE_MSN)),1)...
    /sqrt(size(vertcat(allN_spksNoLcE_FSI,allN_spksNoLcE_MSN),1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.216 .635 .349],'facealpha',0.5,'EdgeColor','none')
 
xlim(xLim)
ylim([0.6 yMaxSpks])
xticks(0:1:7)
yticks(0:.2:1.6)
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
ylabel('Normalized firing rate', 'FontWeight', 'bold','FontSize', 13);
% ylabel({'\fontsize{13}Normalized firing rate - \color[rgb]{0.216, 0.635, 0.349}All\color{black},\color[rgb]{.988 .263 .4}FSI\color{black},\color[rgb]{.682 .698 1}MSN'}, 'FontWeight', 'bold');
yyaxis right;
x = linspace(0, 7, length(mediaFinalE));
plot(x, smooth(mediaFinalE),'k','lineW',2);
ylim([0 35])
h = ylabel('\fontsize{13}Mean speed (cm/s)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [0.4, 0, 0]); % Ajuste a posição conforme necessário
xlabel('Time (s)','FontWeight', 'bold','FontSize', 12);

xt = xticks; % Obter os valores dos xticks
for i = 2:length(xt)-1
    xline(xt(i), '--', 'color',[.8 .8 .8 .2]); % Adiciona uma linha vertical em cada xtick
end
A = smooth(mean(matrizE, 2),10)+...
    std(mean(matrizE, 2),1)...
    /sqrt(size(somaE,1));
B = smooth(mean(matrizE, 2),10)-...
    std(mean(matrizE, 2),1)...
    /sqrt(size(somaE,1));
fill([x fliplr(x)], ...
     [A' fliplr(B')],...
     [.51 .51 .51],'facealpha',0.2,'EdgeColor','none', 'HandleVisibility', 'off')
 
text(-.16, 1.05, 'D', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

f = gcf;
%% SAVE IMAGE
clc
name = 'figure5_SpkStartEndMovment_v1';
saveName = strcat(name, '.tiff');
filename = fullfile('./codes/paperFigures_v1.0/figures/', saveName);
exportgraphics(f, filename, 'Resolution',300)

'done'
%%
clf, clc
fig1 = gcf;
linhas = 2;
colunas = 2;
set(gcf,'color','w')
set(fig1,'position',[0, 0, 1400, 1050]);
pos1 = [0.13,0.489517819706499,0.775,0.435482180293501];
yMaxSpks = 1.5;
yMaxVel = 50;
yticksVel = [0:5:yMaxVel];
xLim = [0 7];

subplot(linhas, colunas,1)
title('Neurons data', 'FontWeight', 'bold','FontSize', 18,'color','k');
hold on;
yyaxis left;
plot(edges7s(1:end-1),smooth(exempleLC_FSI_S,10),'color',[.988 .263 .4],'lineW',2,'lineS','-')
plot(edges7s(1:end-1),smooth(exempleLC_MSN_S,10),'color',[.682 .698 1],'lineW',2,'lineS','-')
plot(edges7s(1:end-1),smooth(exempleNoLCFSI_S,10),'color',[.988 .263 .4],'lineW',2, 'lineS',':')
plot(edges7s(1:end-1),smooth(exempleNoLCMSN_S,10),'color',[.682 .698 1],'lineW',2, 'lineS',':')
xlim(xLim)
ylim([0 1])
xticks(0:1:7)
xticklabels([])
yticks(0:.2:1)
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
ylabel({'\fontsize{18}Starting locomotion', '\fontsize{13}Normalized firing rate'}, 'FontWeight', 'bold');
yyaxis right;
x = linspace(0, 7, length(allMeanVelocityStart{1, 4}));
plot(x, smooth(allMeanVelocityStart{1, 4},10),'k','lineW',2)
ylim([0 35])
h = ylabel('\fontsize{13}Mean speed (cm/s)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [0.4, 0, 0]); % Ajuste a posição conforme necessário

xt = xticks; % Obter os valores dos xticks
for i = 1:length(xt)-1
    xline(xt(i)+1, '--', 'color',[.8 .8 .8 .2], 'HandleVisibility', 'off'); % Adiciona uma linha vertical em cada xtick
end
legend({'LC FSI','LC MSN','NoLC FSI','NoLC MSN','Mean speed'})
set(legend,...
    'Position',[0.44 0.5 0.16428571156093 0.0526315775268009],...
    'NumColumns',3,...
    'FontWeight','bold',...
    'FontSize',12);

subplot(linhas, colunas,2)
title('Group data', 'FontWeight', 'bold','FontSize', 18,'color','k');
hold on;
yyaxis left;
plot(edges7s(1:end-1),smooth(mean(allN_spksLcS_FSI),10),'color',[.988 .263 .4],'lineW',2, 'lineS','-')
plot(edges7s(1:end-1),smooth(mean(allN_spksLcS_MSN),10),'color',[.682 .698 1],'lineW',2, 'lineS','-')
plot(edges7s(1:end-1),smooth(mean(allN_spksNoLcS_FSI),10),'color',[.988 .263 .4],'lineW',2, 'lineS',':')
plot(edges7s(1:end-1),smooth(mean(allN_spksNoLcS_MSN),10),'color',[.682 .698 1],'lineW',2, 'lineS',':')

A = smooth(mean(allN_spksLcS_FSI),10)+...
    std(mean(allN_spksLcS_FSI),1)...
    /sqrt(size(allN_spksLcS_FSI,1));
B = smooth(mean(allN_spksLcS_FSI),10)-...
    std(mean(allN_spksLcS_FSI),1)...
    /sqrt(size(allN_spksLcS_FSI,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.988 .263 .4],'facealpha',0.2,'EdgeColor','none')
 
A = smooth(mean(allN_spksLcS_MSN),10)+...
    std(mean(allN_spksLcS_MSN),1)...
    /sqrt(size(allN_spksLcS_MSN,1));
B = smooth(mean(allN_spksLcS_MSN),10)-...
    std(mean(allN_spksLcS_MSN),1)...
    /sqrt(size(allN_spksLcS_MSN,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.682 .698 1],'facealpha',0.2,'EdgeColor','none')
 
A = smooth(mean(allN_spksNoLcS_FSI),10)+...
    std(mean(allN_spksNoLcS_FSI),1)...
    /sqrt(size(allN_spksNoLcS_FSI,1));
B = smooth(mean(allN_spksNoLcS_FSI),10)-...
    std(mean(allN_spksNoLcS_FSI),1)...
    /sqrt(size(allN_spksNoLcS_FSI,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.988 .263 .4],'facealpha',0.2,'EdgeColor','none')
 
A = smooth(mean(allN_spksNoLcS_MSN),10)+...
    std(mean(allN_spksNoLcS_MSN),1)...
    /sqrt(size(allN_spksNoLcS_MSN,1));
B = smooth(mean(allN_spksNoLcS_MSN),10)-...
    std(mean(allN_spksNoLcS_MSN),1)...
    /sqrt(size(allN_spksNoLcS_MSN,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.682 .698 1],'facealpha',0.2,'EdgeColor','none')
 
xlim(xLim)
ylim([.5 yMaxSpks])
xticks(0:1:7)
xticklabels([])
yticks(0:.2:1)
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
ylabel('Mean firing rate', 'FontWeight', 'bold','FontSize', 13);

yyaxis right;
x = linspace(0, 7, length(mediaFinalS));
plot(x, smooth(mediaFinalS),'k','lineW',2); 
ylim([0 35]);
h = ylabel('\fontsize{13}Mean speed (cm/s)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [0.4, 0, 0]); % Ajuste a posição conforme necessário
xt = xticks; % Obter os valores dos xticks
for i = 1:length(xt)-1
    xline(xt(i)+1, '--', 'color',[.8 .8 .8 .2], 'HandleVisibility', 'off'); % Adiciona uma linha vertical em cada xtick
end


subplot(linhas, colunas,3)
hold on;
yyaxis left;
plot(edges7s(1:end-1),smooth((exempleLC_FSI_E),10),'color',[.988 .263 .4],'lineW',2,'lineS','-')
plot(edges7s(1:end-1),smooth((exempleLC_MSN_E),10),'color',[.682 .698 1],'lineW',2,'lineS','-')
plot(edges7s(1:end-1),smooth((exempleNoLCFSI_E),10),'color',[.988 .263 .4],'lineW',2, 'lineS',':')
plot(edges7s(1:end-1),smooth((exempleNoLCMSN_E),10),'color',[.682 .698 1],'lineW',2, 'lineS',':')

xlim(xLim)
ylim([0 1])
xticks(0:1:7)
yticks(0:.2:1)
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
ylabel({'\fontsize{18}Finishing locomotion', '\fontsize{13}Normalized firing rate'}, 'FontWeight', 'bold');

yyaxis right;
x = linspace(0, 7, length(allMeanVelocityEnd{1, 4}));
plot(x, smooth(allMeanVelocityEnd{1, 4},10),'k','lineW',2)
ylim([0 35])
h = ylabel('\fontsize{13}Mean speed (cm/s)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [0.4, 0, 0]); % Ajuste a posição conforme necessário
xt = xticks; % Obter os valores dos xticks
for i = 1:length(xt)-1
    xline(xt(i)+1, '--', 'color',[.8 .8 .8 .2], 'HandleVisibility', 'off'); % Adiciona uma linha vertical em cada xtick
end
xlim(xLim)
xlabel('Time (s)','FontWeight', 'bold','FontSize', 12);


% xlabel('Time (s)','FontWeight', 'bold','FontSize', 12);
% 

subplot(linhas, colunas,4)
hold on;
yyaxis left;
plot(edges7s(1:end-1),smooth(mean(allN_spksLcE_FSI),10),'color',[.988 .263 .4],'lineW',2, 'lineS','-');
plot(edges7s(1:end-1),smooth(mean(allN_spksLcE_MSN),10),'color',[.682 .698 1],'lineW',2, 'lineS','-');
plot(edges7s(1:end-1),smooth(mean(allN_spksNoLcE_FSI),10),'color',[.988 .263 .4],'lineW',2, 'lineS',':');
plot(edges7s(1:end-1),smooth(mean(allN_spksNoLcE_MSN),10),'color',[.682 .698 1],'lineW',2, 'lineS',':');

A = smooth(mean(allN_spksLcE_FSI),10)+...
    std(mean(allN_spksLcE_FSI),1)...
    /sqrt(size(allN_spksLcE_FSI,1));
B = smooth(mean(allN_spksLcE_FSI),10)-...
    std(mean(allN_spksLcE_FSI),1)...
    /sqrt(size(allN_spksLcE_FSI,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.988 .263 .4],'facealpha',0.2,'EdgeColor','none')
 
A = smooth(mean(allN_spksLcE_MSN),10)+...
    std(mean(allN_spksLcE_MSN),1)...
    /sqrt(size(allN_spksLcE_MSN,1));
B = smooth(mean(allN_spksLcE_MSN),10)-...
    std(mean(allN_spksLcE_MSN),1)...
    /sqrt(size(allN_spksLcE_MSN,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.682 .698 1],'facealpha',0.2,'EdgeColor','none')
 
A = smooth(mean(allN_spksNoLcE_FSI),10)+...
    std(mean(allN_spksNoLcE_FSI),1)...
    /sqrt(size(allN_spksNoLcE_FSI,1));
B = smooth(mean(allN_spksNoLcE_FSI),10)-...
    std(mean(allN_spksNoLcE_FSI),1)...
    /sqrt(size(allN_spksNoLcE_FSI,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.988 .263 .4],'facealpha',0.2,'EdgeColor','none')
 
A = smooth(mean(allN_spksNoLcE_MSN),10)+...
    std(mean(allN_spksNoLcE_MSN),1)...
    /sqrt(size(allN_spksNoLcE_MSN,1));
B = smooth(mean(allN_spksNoLcE_MSN),10)-...
    std(mean(allN_spksNoLcE_MSN),1)...
    /sqrt(size(allN_spksNoLcE_MSN,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.682 .698 1],'facealpha',0.2,'EdgeColor','none')
 
ylim([0 yMaxSpks]);
xlim(xLim);
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
ylabel('Mean firing rate', 'FontWeight', 'bold','FontSize', 13);
xticks(0:1:7);
yticks(0:.2:1)
yticklabels();

yyaxis right;
x = linspace(0, 7, length(mediaFinalE));
plot(x, smooth(mediaFinalE),'k','lineW',2);
ylim([0 35])
h = ylabel('\fontsize{13}Mean speed (cm/s)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [0.4, 0, 0]); % Ajuste a posição conforme necessário
xlabel('Time (s)','FontWeight', 'bold','FontSize', 12);

xt = xticks; % Obter os valores dos xticks
for i = 1:length(xt)-1
    xline(xt(i)+1, '--', 'color',[.8 .8 .8 .2]); % Adiciona uma linha vertical em cada xtick
end

f = gcf;
%% SAVE IMAGE
clc
name = 'figure5_SpkStartEndMovment_v1_';
saveName = strcat(name, '.tiff');
filename = fullfile('./codes/paperFigures/figures/', saveName);
exportgraphics(f, filename, 'Resolution',300)

'done'
%%
clf, clc
fig1 = gcf;
linhas = 2;
colunas = 4;
set(gcf,'color','w')
set(fig1,'position',[0, 0, 1400, 1050]);
pos1 = [0.13,0.489517819706499,0.775,0.435482180293501];
yMaxSpks = 1;
yMaxVel = 50;
yticksVel = [0:5:yMaxVel];
xLim = [0 7];

subplot(linhas, colunas,1)
title('FSI neurons data', 'FontWeight', 'bold','FontSize', 18,'color','k');
hold on;
yyaxis left;
plot(edges7s(1:end-1),smooth(exempleLC_FSI_S,10),'color',[.988 .263 .4],'lineW',2,'lineS','-')
plot(edges7s(1:end-1)-10,smooth(exempleLC_MSN_S,10),'color',[.682 .698 1],'lineW',2,'lineS','-')
plot(edges7s(1:end-1),smooth(exempleNoLCFSI_S,10),'color',[.988 .263 .4],'lineW',2, 'lineS',':')
plot(edges7s(1:end-1)-10,smooth(exempleNoLCMSN_S,10),'color',[.682 .698 1],'lineW',2, 'lineS',':')
xlim(xLim)
ylim([0 1])
xticks(0:1:7)
xticklabels([])
yticks(0:.2:1)
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
ylabel({'\fontsize{18}Starting locomotion', '\fontsize{13}Normalized firing rate'}, 'FontWeight', 'bold');
yyaxis right;
x = linspace(0, 7, length(allMeanVelocityStart{1, 4}));
plot(x, smooth(allMeanVelocityStart{1, 4},10),'k','lineW',2)
ylim([0 35])
yticks(0:5:35)
yticklabels([])
% h = ylabel('\fontsize{13}Mean speed (cm/s)', 'FontWeight', 'bold','Rotation', -90);
% set(h, 'Position', get(h, 'Position') + [0.4, 0, 0]); % Ajuste a posição conforme necessário

xt = xticks; % Obter os valores dos xticks
for i = 1:length(xt)-1
    xline(xt(i)+1, '--', 'color',[.8 .8 .8 .2], 'HandleVisibility', 'off'); % Adiciona uma linha vertical em cada xtick
end

legend({'LC FSI','LC MSN','NoLC FSI','NoLC MSN','Mean speed'})
set(legend,...
    'Position',[0.44 0.5 0.16428571156093 0.0526315775268009],...
    'NumColumns',5,...
    'FontWeight','bold',...
    'FontSize',12);

subplot(linhas, colunas,2)
title('MSN neurons data', 'FontWeight', 'bold','FontSize', 18,'color','k');
hold on;
yyaxis left;
% plot(edges7s(1:end-1),smooth(exempleLC_FSI_S,10),'color',[.988 .263 .4],'lineW',2,'lineS','-')
plot(edges7s(1:end-1),smooth(exempleLC_MSN_S,10),'color',[.682 .698 1],'lineW',2,'lineS','-')
% plot(edges7s(1:end-1),smooth(exempleNoLCFSI_S,10),'color',[.988 .263 .4],'lineW',2, 'lineS',':')
plot(edges7s(1:end-1),smooth(exempleNoLCMSN_S,10),'color',[.682 .698 1],'lineW',2, 'lineS',':')
xlim(xLim)
ylim([0 1])
xticks(0:1:7)
xticklabels([])
yticks(0:.2:1)
yticklabels([])
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
yyaxis right;
x = linspace(0, 7, length(allMeanVelocityStart{1, 4}));
plot(x, smooth(allMeanVelocityStart{1, 4},10),'k','lineW',2)
ylim([0 35])
yticks(0:5:35)
yticklabels([])

xt = xticks; % Obter os valores dos xticks
for i = 1:length(xt)-1
    xline(xt(i)+1, '--', 'color',[.8 .8 .8 .2], 'HandleVisibility', 'off'); % Adiciona uma linha vertical em cada xtick
end

subplot(linhas, colunas,3)
title('FSI group data', 'FontWeight', 'bold','FontSize', 18,'color','k');
hold on;
yyaxis left;
plot(edges7s(1:end-1),smooth(mean(allN_spksLcS_FSI),10),'color',[.988 .263 .4],'lineW',2, 'lineS','-')
plot(edges7s(1:end-1),smooth(mean(allN_spksNoLcS_FSI),10),'color',[.988 .263 .4],'lineW',2, 'lineS',':')

A = smooth(mean(allN_spksLcS_FSI),10)+...
    std(mean(allN_spksLcS_FSI),1)...
    /sqrt(size(allN_spksLcS_FSI,1));
B = smooth(mean(allN_spksLcS_FSI),10)-...
    std(mean(allN_spksLcS_FSI),1)...
    /sqrt(size(allN_spksLcS_FSI,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.988 .263 .4],'facealpha',0.2,'EdgeColor','none')

A = smooth(mean(allN_spksNoLcS_FSI),10)+...
    std(mean(allN_spksNoLcS_FSI),1)...
    /sqrt(size(allN_spksNoLcS_FSI,1));
B = smooth(mean(allN_spksNoLcS_FSI),10)-...
    std(mean(allN_spksNoLcS_FSI),1)...
    /sqrt(size(allN_spksNoLcS_FSI,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.988 .263 .4],'facealpha',0.2,'EdgeColor','none')
 
xlim(xLim)
ylim([0 yMaxSpks])
xticks(0:1:7)
xticklabels([])
yticks(0:.2:1)
yticklabels([])
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
yyaxis right;
x = linspace(0, 7, length(mediaFinalS));
plot(x, smooth(mediaFinalS),'k','lineW',2); 
ylim([0 35]);
yticks([0:5:35]);
yticklabels([])
xt = xticks; % Obter os valores dos xticks
for i = 1:length(xt)-1
    xline(xt(i)+1, '--', 'color',[.8 .8 .8 .2], 'HandleVisibility', 'off'); % Adiciona uma linha vertical em cada xtick
end

subplot(linhas, colunas,4)
title('MSN group data', 'FontWeight', 'bold','FontSize', 18,'color','k');
hold on;
yyaxis left;
plot(edges7s(1:end-1),smooth(mean(allN_spksLcS_MSN),10),'color',[.682 .698 1],'lineW',2, 'lineS','-')
plot(edges7s(1:end-1),smooth(mean(allN_spksNoLcS_MSN),10),'color',[.682 .698 1],'lineW',2, 'lineS',':')

A = smooth(mean(allN_spksLcS_MSN),10)+...
    std(mean(allN_spksLcS_MSN),1)...
    /sqrt(size(allN_spksLcS_MSN,1));
B = smooth(mean(allN_spksLcS_MSN),10)-...
    std(mean(allN_spksLcS_MSN),1)...
    /sqrt(size(allN_spksLcS_MSN,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.682 .698 1],'facealpha',0.2,'EdgeColor','none')

A = smooth(mean(allN_spksNoLcS_MSN),10)+...
    std(mean(allN_spksNoLcS_MSN),1)...
    /sqrt(size(allN_spksNoLcS_MSN,1));
B = smooth(mean(allN_spksNoLcS_MSN),10)-...
    std(mean(allN_spksNoLcS_MSN),1)...
    /sqrt(size(allN_spksNoLcS_MSN,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.682 .698 1],'facealpha',0.2,'EdgeColor','none')
 
xlim(xLim)
ylim([0 yMaxSpks])
xticks(0:1:7)
xticklabels([])
yticks(0:.2:1)
yticklabels([])
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
% ylabel('Mean firing rate', 'FontWeight', 'bold','FontSize', 13);

yyaxis right;
x = linspace(0, 7, length(mediaFinalS));
plot(x, smooth(mediaFinalS),'k','lineW',2); 
ylim([0 35]);
h = ylabel('\fontsize{13}Mean speed (cm/s)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [0.8, 0, 0]); % Ajuste a posição conforme necessário
xt = xticks; % Obter os valores dos xticks
for i = 1:length(xt)-1
    xline(xt(i)+1, '--', 'color',[.8 .8 .8 .2], 'HandleVisibility', 'off'); % Adiciona uma linha vertical em cada xtick
end


subplot(linhas, colunas,5)
hold on;
yyaxis left;
plot(edges7s(1:end-1),smooth(exempleLC_FSI_E,10),'color',[.988 .263 .4],'lineW',2,'lineS','-')
% plot(edges7s(1:end-1),smooth(exempleLC_MSN_E,10),'color',[.682 .698 1],'lineW',2,'lineS','-')
plot(edges7s(1:end-1),smooth(exempleNoLCFSI_E,10),'color',[.988 .263 .4],'lineW',2, 'lineS',':')
% plot(edges7s(1:end-1),smooth(exempleNoLCMSN_E,10),'color',[.682 .698 1],'lineW',2, 'lineS',':')
xlim(xLim)
ylim([0 1])
xticks(0:1:7)
yticks(0:.2:1)
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
ylabel({'\fontsize{18}Finishing locomotion', '\fontsize{13}Normalized firing rate'}, 'FontWeight', 'bold');
yyaxis right;
x = linspace(0, 7, length(allMeanVelocityEnd{1, 4}));
plot(x, smooth(allMeanVelocityEnd{1, 4},10),'k','lineW',2)
ylim([0 35])
yticks(0:5:35)
yticklabels([])
xlabel('Time (s)','FontWeight', 'bold','FontSize', 12);
xt = xticks; % Obter os valores dos xticks
for i = 1:length(xt)-1
    xline(xt(i)+1, '--', 'color',[.8 .8 .8 .2], 'HandleVisibility', 'off'); % Adiciona uma linha vertical em cada xtick
end

subplot(linhas, colunas,6)
hold on;
yyaxis left;
plot(edges7s(1:end-1),smooth(exempleLC_MSN_E,10),'color',[.682 .698 1],'lineW',2,'lineS','-')
plot(edges7s(1:end-1),smooth(exempleNoLCMSN_E,10),'color',[.682 .698 1],'lineW',2, 'lineS',':')
xlim(xLim)
ylim([0 1])
xticks(0:1:7)
yticks(0:.2:1)
yticklabels([])
xlabel('Time (s)','FontWeight', 'bold','FontSize', 12);
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
yyaxis right;
x = linspace(0, 7, length(allMeanVelocityEnd{1, 4}));
plot(x, smooth(allMeanVelocityEnd{1, 4},10),'k','lineW',2)
ylim([0 35])
yticks(0:5:35)
yticklabels([])

xt = xticks; % Obter os valores dos xticks
for i = 1:length(xt)-1
    xline(xt(i)+1, '--', 'color',[.8 .8 .8 .2], 'HandleVisibility', 'off'); % Adiciona uma linha vertical em cada xtick
end

subplot(linhas, colunas,7)
hold on;
yyaxis left;
plot(edges7s(1:end-1),smooth(mean(allN_spksLcE_FSI),10),'color',[.988 .263 .4],'lineW',2, 'lineS','-')
plot(edges7s(1:end-1),smooth(mean(allN_spksNoLcE_FSI),10),'color',[.988 .263 .4],'lineW',2, 'lineS',':')

A = smooth(mean(allN_spksLcE_FSI),10)+...
    std(mean(allN_spksLcE_FSI),1)...
    /sqrt(size(allN_spksLcE_FSI,1));
B = smooth(mean(allN_spksLcE_FSI),10)-...
    std(mean(allN_spksLcE_FSI),1)...
    /sqrt(size(allN_spksLcE_FSI,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.988 .263 .4],'facealpha',0.2,'EdgeColor','none')

A = smooth(mean(allN_spksNoLcE_FSI),10)+...
    std(mean(allN_spksNoLcE_FSI),1)...
    /sqrt(size(allN_spksNoLcE_FSI,1));
B = smooth(mean(allN_spksNoLcE_FSI),10)-...
    std(mean(allN_spksNoLcE_FSI),1)...
    /sqrt(size(allN_spksNoLcE_FSI,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.988 .263 .4],'facealpha',0.2,'EdgeColor','none')
 
xlim(xLim)
ylim([0 yMaxSpks])
xticks(0:1:7)
yticks(0:.2:1)
yticklabels([])
xlabel('Time (s)','FontWeight', 'bold','FontSize', 12);
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
yyaxis right;
x = linspace(0, 7, length(mediaFinalE));
plot(x, smooth(mediaFinalE),'k','lineW',2); 
ylim([0 35]);
yticks([0:5:35]);
yticklabels([])
xt = xticks; % Obter os valores dos xticks
for i = 1:length(xt)-1
    xline(xt(i)+1, '--', 'color',[.8 .8 .8 .2], 'HandleVisibility', 'off'); % Adiciona uma linha vertical em cada xtick
end

subplot(linhas, colunas,8)
hold on;
yyaxis left;
plot(edges7s(1:end-1),smooth(mean(allN_spksLcE_MSN),10),'color',[.682 .698 1],'lineW',2, 'lineS','-')
plot(edges7s(1:end-1),smooth(mean(allN_spksNoLcE_MSN),10),'color',[.682 .698 1],'lineW',2, 'lineS',':')

A = smooth(mean(allN_spksLcE_MSN),10)+...
    std(mean(allN_spksLcE_MSN),1)...
    /sqrt(size(allN_spksLcE_MSN,1));
B = smooth(mean(allN_spksLcE_MSN),10)-...
    std(mean(allN_spksLcE_MSN),1)...
    /sqrt(size(allN_spksLcE_MSN,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.682 .698 1],'facealpha',0.2,'EdgeColor','none')

A = smooth(mean(allN_spksNoLcE_MSN),10)+...
    std(mean(allN_spksNoLcE_MSN),1)...
    /sqrt(size(allN_spksNoLcE_MSN,1));
B = smooth(mean(allN_spksNoLcE_MSN),10)-...
    std(mean(allN_spksNoLcE_MSN),1)...
    /sqrt(size(allN_spksNoLcE_MSN,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.682 .698 1],'facealpha',0.2,'EdgeColor','none')
 
xlim(xLim)
ylim([0 yMaxSpks])
xticks(0:1:7)
yticks(0:.2:1)
yticklabels([])
xlabel('Time (s)','FontWeight', 'bold','FontSize', 12);
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
% ylabel('Mean firing rate', 'FontWeight', 'bold','FontSize', 13);

yyaxis right;
x = linspace(0, 7, length(mediaFinalE));
plot(x, smooth(mediaFinalE),'k','lineW',2); 
ylim([0 35]);
h = ylabel('\fontsize{13}Mean speed (cm/s)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [0.8, 0, 0]); % Ajuste a posição conforme necessário
xt = xticks; % Obter os valores dos xticks
for i = 1:length(xt)-1
    xline(xt(i)+1, '--', 'color',[.8 .8 .8 .2], 'HandleVisibility', 'off'); % Adiciona uma linha vertical em cada xtick
end

% 
f = gcf;
%% SAVE IMAGE
clc
name = 'figure5_SpkStartEndMovment_v2_';
saveName = strcat(name, '.tiff');
filename = fullfile('./codes/paperFigures/figures/', saveName);
exportgraphics(f, filename, 'Resolution',300)

'done'
%%
clf, clc
fig1 = gcf;
linhas = 2;
colunas = 2;
set(gcf,'color','w')
set(fig1,'position',[0, 0, 1400, 1050]);
pos1 = [0.13,0.489517819706499,0.775,0.435482180293501];
yMaxSpks = 1;
yMaxVel = 50;
yticksVel = [0:5:yMaxVel];
xLim = [0 7];

subplot(linhas, colunas,1)
title('FSI', 'FontWeight', 'bold','FontSize', 18,'color','k');
hold on;
yyaxis left;
plot(edges7s(1:end-1),smooth(exempleLC_FSI_S,10),'color',[.988 .263 .4],'lineW',2,'lineS','-')
plot(edges7s(1:end-1),smooth(exempleNoLCFSI_S,10),'color',[.988 .263 .4],'lineW',2, 'lineS',':')
plot(edges7s(1:end-1),smooth(mean(allN_spksLcS_FSI),10),'color',[.988 .263 .4],'lineW',2, 'lineS','-', 'HandleVisibility', 'off')
plot(edges7s(1:end-1),smooth(mean(allN_spksNoLcS_FSI),10),'color',[.988 .263 .4],'lineW',2, 'lineS',':', 'HandleVisibility', 'off')

A = smooth(mean(allN_spksLcS_FSI),10)+...
    std(mean(allN_spksLcS_FSI),1)...
    /sqrt(size(allN_spksLcS_FSI,1));
B = smooth(mean(allN_spksLcS_FSI),10)-...
    std(mean(allN_spksLcS_FSI),1)...
    /sqrt(size(allN_spksLcS_FSI,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.988 .263 .4],'facealpha',0.2,'EdgeColor','none')
 
A = smooth(mean(allN_spksNoLcS_FSI),10)+...
    std(mean(allN_spksNoLcS_FSI),1)...
    /sqrt(size(allN_spksNoLcS_FSI,1));
B = smooth(mean(allN_spksNoLcS_FSI),10)-...
    std(mean(allN_spksNoLcS_FSI),1)...
    /sqrt(size(allN_spksNoLcS_FSI,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.988 .263 .4],'facealpha',0.2,'EdgeColor','none')
xlim(xLim)
ylim([0 1])
xticks(0:1:7)
xticklabels([])
yticks(0:.2:1)
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
ylabel({'\fontsize{18}Starting locomotion', '\fontsize{13}Normalized firing rate'}, 'FontWeight', 'bold');
yyaxis right;
x = linspace(0, 7, length(allMeanVelocityStart{1, 4}));
plot(x, smooth(allMeanVelocityStart{1, 4},10),'k','lineW',2)
ylim([0 35])
h = ylabel('\fontsize{13}Mean speed (cm/s)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [0.4, 0, 0]); % Ajuste a posição conforme necessário

xt = xticks; % Obter os valores dos xticks
for i = 1:length(xt)-2
    xline(xt(i)+1, '--', 'color',[.8 .8 .8 .2], 'HandleVisibility', 'off'); % Adiciona uma linha vertical em cada xtick
end
legend({'LC neuron','NoLC neuron','LC group data','NoLC group data','Mean speed'})
set(legend,...
    'Position',[0.44 0.5 0.16428571156093 0.0526315775268009],...
    'NumColumns',3,...
    'FontWeight','bold',...
    'FontSize',12);

subplot(linhas, colunas,2)
title('MSN', 'FontWeight', 'bold','FontSize', 18,'color','k');
hold on;
yyaxis left;
plot(edges7s(1:end-1),smooth(exempleLC_MSN_S,10),'color',[.949 .875 0],'lineW',2,'lineS','-')
plot(edges7s(1:end-1),smooth(exempleNoLCMSN_S,10),'color',[.949 .875 0],'lineW',2, 'lineS',':')
plot(edges7s(1:end-1),smooth(mean(allN_spksLcS_MSN),10),'color',[.682 .698 1],'lineW',2, 'lineS','-')
plot(edges7s(1:end-1),smooth(mean(allN_spksNoLcS_MSN),10),'color',[.682 .698 1],'lineW',2, 'lineS',':')

A = smooth(mean(allN_spksLcS_MSN),10)+...
    std(mean(allN_spksLcS_MSN),1)...
    /sqrt(size(allN_spksLcS_MSN,1));
B = smooth(mean(allN_spksLcS_MSN),10)-...
    std(mean(allN_spksLcS_MSN),1)...
    /sqrt(size(allN_spksLcS_MSN,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.682 .698 1],'facealpha',0.2,'EdgeColor','none')
 
A = smooth(mean(allN_spksNoLcS_MSN),10)+...
    std(mean(allN_spksNoLcS_MSN),1)...
    /sqrt(size(allN_spksNoLcS_MSN,1));
B = smooth(mean(allN_spksNoLcS_MSN),10)-...
    std(mean(allN_spksNoLcS_MSN),1)...
    /sqrt(size(allN_spksNoLcS_MSN,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.682 .698 1],'facealpha',0.2,'EdgeColor','none')
 
xlim(xLim)
ylim([0 yMaxSpks])
xticks(0:1:7)
xticklabels([])
yticks(0:.2:1)
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
ylabel('Mean firing rate', 'FontWeight', 'bold','FontSize', 13);

yyaxis right;
x = linspace(0, 7, length(mediaFinalS));
plot(x, smooth(mediaFinalS),'k','lineW',2); 
ylim([0 35]);
h = ylabel('\fontsize{13}Mean speed (cm/s)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [0.4, 0, 0]); % Ajuste a posição conforme necessário
xt = xticks; % Obter os valores dos xticks
for i = 1:length(xt)-2
    xline(xt(i)+1, '--', 'color',[.8 .8 .8 .2], 'HandleVisibility', 'off'); % Adiciona uma linha vertical em cada xtick
end


subplot(linhas, colunas,3)
hold on;
yyaxis left;
plot(edges7s(1:end-1),smooth((exempleLC_FSI_E),10),'color',[.949 .875 0],'lineW',2,'lineS','-')
plot(edges7s(1:end-1),smooth((exempleNoLCFSI_E),10),'color',[.949 .875 0],'lineW',2, 'lineS',':')
plot(edges7s(1:end-1),smooth(mean(allN_spksLcE_FSI),10),'color',[.988 .263 .4],'lineW',2, 'lineS','-');
plot(edges7s(1:end-1),smooth(mean(allN_spksNoLcE_FSI),10),'color',[.988 .263 .4],'lineW',2, 'lineS',':');

A = smooth(mean(allN_spksLcE_FSI),10)+...
    std(mean(allN_spksLcE_FSI),1)...
    /sqrt(size(allN_spksLcE_FSI,1));
B = smooth(mean(allN_spksLcE_FSI),10)-...
    std(mean(allN_spksLcE_FSI),1)...
    /sqrt(size(allN_spksLcE_FSI,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.988 .263 .4],'facealpha',0.2,'EdgeColor','none')
 
A = smooth(mean(allN_spksNoLcE_FSI),10)+...
    std(mean(allN_spksNoLcE_FSI),1)...
    /sqrt(size(allN_spksNoLcE_FSI,1));
B = smooth(mean(allN_spksNoLcE_FSI),10)-...
    std(mean(allN_spksNoLcE_FSI),1)...
    /sqrt(size(allN_spksNoLcE_FSI,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.988 .263 .4],'facealpha',0.2,'EdgeColor','none')
 
xlim(xLim)
ylim([0 1])
xticks(0:1:7)
yticks(0:.2:1)
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
ylabel({'\fontsize{18}Finishing locomotion', '\fontsize{13}Normalized firing rate'}, 'FontWeight', 'bold');

yyaxis right;
x = linspace(0, 7, length(allMeanVelocityEnd{1, 4}));
plot(x, smooth(allMeanVelocityEnd{1, 4},10),'k','lineW',2)
ylim([0 35])
h = ylabel('\fontsize{13}Mean speed (cm/s)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [0.4, 0, 0]); % Ajuste a posição conforme necessário
xt = xticks; % Obter os valores dos xticks
for i = 1:length(xt)-2
    xline(xt(i)+1, '--', 'color',[.8 .8 .8 .2], 'HandleVisibility', 'off'); % Adiciona uma linha vertical em cada xtick
end
xlim(xLim)
xlabel('Time (s)','FontWeight', 'bold','FontSize', 12);


% xlabel('Time (s)','FontWeight', 'bold','FontSize', 12);
% 

subplot(linhas, colunas,4)
hold on;
yyaxis left;
plot(edges7s(1:end-1),smooth(mean(allN_spksLcE_MSN),10),'color',[.682 .698 1],'lineW',2, 'lineS','-');
plot(edges7s(1:end-1),smooth((exempleLC_MSN_E),10),'color',[.949 .875 0],'lineW',2,'lineS','-')
plot(edges7s(1:end-1),smooth((exempleNoLCMSN_E),10),'color',[.949 .875 0],'lineW',2, 'lineS',':')
plot(edges7s(1:end-1),smooth(mean(allN_spksNoLcE_MSN),10),'color',[.682 .698 1],'lineW',2, 'lineS',':');


A = smooth(mean(allN_spksLcE_MSN),10)+...
    std(mean(allN_spksLcE_MSN),1)...
    /sqrt(size(allN_spksLcE_MSN,1));
B = smooth(mean(allN_spksLcE_MSN),10)-...
    std(mean(allN_spksLcE_MSN),1)...
    /sqrt(size(allN_spksLcE_MSN,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.682 .698 1],'facealpha',0.2,'EdgeColor','none')
  
A = smooth(mean(allN_spksNoLcE_MSN),10)+...
    std(mean(allN_spksNoLcE_MSN),1)...
    /sqrt(size(allN_spksNoLcE_MSN,1));
B = smooth(mean(allN_spksNoLcE_MSN),10)-...
    std(mean(allN_spksNoLcE_MSN),1)...
    /sqrt(size(allN_spksNoLcE_MSN,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.682 .698 1],'facealpha',0.2,'EdgeColor','none')
 
ylim([0 yMaxSpks]);
xlim(xLim);
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
ylabel('Mean firing rate', 'FontWeight', 'bold','FontSize', 13);
xticks(0:1:7);
yticks(0:.2:1)
yticklabels();

yyaxis right;
x = linspace(0, 7, length(mediaFinalE));
plot(x, smooth(mediaFinalE),'k','lineW',2);
ylim([0 35])
h = ylabel('\fontsize{13}Mean speed (cm/s)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [0.4, 0, 0]); % Ajuste a posição conforme necessário
xlabel('Time (s)','FontWeight', 'bold','FontSize', 12);

xt = xticks; % Obter os valores dos xticks
for i = 1:length(xt)-2
    xline(xt(i)+1, '--', 'color',[.8 .8 .8 .2]); % Adiciona uma linha vertical em cada xtick
end

f = gcf;
%% SAVE IMAGE
clc
name = 'figure5_SpkStartEndMovment_v3_';
saveName = strcat(name, '.tiff');
filename = fullfile('./codes/paperFigures/figures/', saveName);
exportgraphics(f, filename, 'Resolution',300)

'done'
%% PEARSONS CORRELATION
clc

a = smooth(mean(vertcat(allN_spksLcS_FSI,allN_spksLcS_MSN)));
b = smooth(mean(matrizS, 2));
indices = round(linspace(1, length(b), 70));
b_sub = b(indices);
[r, p] = corr(b_sub(:), a(:))


%
a = smooth(mean(vertcat(allN_spksNoLcS_FSI,allN_spksNoLcS_MSN)));
b = smooth(mean(matrizS, 2));
indices = round(linspace(1, length(b), 70));
b_sub = b(indices);
[r, p] = corr(b_sub(:), a(:))

%
a = smooth(mean(vertcat(allN_spksLcE_FSI,allN_spksLcE_MSN)));
b = smooth(mean(matrizE, 2));
indices = round(linspace(1, length(b), 70));
b_sub = b(indices);
[r, p] = corr(b_sub(:), a(:))

%
a = smooth(mean(vertcat(allN_spksNoLcE_FSI,allN_spksNoLcE_MSN)));
b = smooth(mean(matrizE, 2));
indices = round(linspace(1, length(b), 70));
b_sub = b(indices);
[r, p] = corr(b_sub(:), a(:))
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
