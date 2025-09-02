%%
clc, clf, clear
cd /Users/phb_lopes/Library/CloudStorage/OneDrive-Pessoal/projetoEletrofisiologiaMacbook
folder = './startEndSpksAndVelocity/7sAnalise/spks/';
allSpks = readData(folder);

folder = './startEndSpksAndVelocity/7sAnalise/meanVelocity/start/';
dataMeanVelocityStart = readData(folder);

folder = './startEndSpksAndVelocity/7sAnalise/meanVelocity/end/';
dataMeanVelocityEnd = readData(folder);

timeInfoSession = 1;   % Size of time window (s)
lag = 3;               % Maximum lag for CCG analysis (s)
overlap = .9;          % Fraction of overlap between windows


folder = strcat('/Users/phb_lopes/Library/CloudStorage/OneDrive-Pessoal/projetoEletrofisiologiaMacbook/',num2str(timeInfoSession),'sInterval/', num2str(lag),'sLag/',num2str(overlap),'_overlap/');
dataFsiAndLocomotionCells_all = readData(folder);

load(strcat('./ccgLagFrSpeedCellsInfoPosAndNeg_1sInterval_0.9overlap_3sLagNoSCwith2.5Lag.mat'))
valorCcg = ccgLagFrLocomotionCellsData{1};
speedCells = ccgLagFrLocomotionCellsData{4};

%%
clc
allMeanVelocityStart = {};
allMeanVelocityEnd = {};

for i=1:length(dataMeanVelocityStart)
    allMeanVelocityStart(i) = {dataMeanVelocityStart{i}.mediaStartMovment};
    allMeanVelocityEnd(i) = {dataMeanVelocityEnd{i}.mediaEndMovment};
end

sumStart = [];
sumEndind = [];
for i=1:length(allMeanVelocityStart)
    sumStart = [sumStart; {allMeanVelocityStart{i}}];
    sumEndind = [sumEndind; {allMeanVelocityEnd{i}}];
end

%%
matrizS = cell2mat(sumStart');
meanVelocityStarting = mean(matrizS, 2);

matrizE = cell2mat(sumEndind');
meanVelocityEnding = mean(matrizE, 2);

%%
clc
allCells = {};
strCells = {};
for i=1:length(dataFsiAndLocomotionCells_all)
    allCells = [allCells; dataFsiAndLocomotionCells_all{1,i}.data];
end

strCells = allCells(cell2mat(allCells(:,6)) >=21,:);
strCells = strCells(speedCells,:);

spkNeurons = {};
for j=1:length(allSpks)
    spkNeurons = [spkNeurons; allSpks{1,j}.spks];
end

indexSC_MSN = [];
indexSC_FSI = [];


for i = 1:size(strCells, 1)
    for j = 1:size(spkNeurons, 1)
        if isequal(strCells(i, [1:3]), spkNeurons(j, [1,2,3]))
            if cell2mat(strCells(i,6)) == 21
                indexSC_MSN = [indexSC_MSN; j]; 
            else
                indexSC_FSI = [indexSC_FSI; j]; 
            end
            break; 
        end
    end
end

%% SC Normalization
clc
binSize = .1;
edges7s = 0:binSize:7;
indexEndVelocity = 7;

spksScS_MSN_pos = [];
spksScS_MSN_neg = [];

spksScS_FSI_pos = [];
spksScS_FSI_neg = [];

spksScE_MSN_pos = [];
spksScE_MSN_neg = [];

spksScE_FSI_pos = [];
spksScE_FSI_neg = [];


allN_spksScS_MSN_pos = [];
allN_spksScS_MSN_neg = [];
allN_spksScS_FSI_pos = [];
allN_spksScS_FSI_neg = [];

allN_spksScE_MSN_pos = [];
allN_spksScE_MSN_neg = [];
allN_spksScE_FSI_pos = [];
allN_spksScE_FSI_neg = [];


for i=1:length(indexSC_MSN)
    
    if valorCcg(indexSC_MSN(i)) > 0 % Positive MSN
        auxS = spkNeurons(indexSC_MSN(i),5);
        auxE = spkNeurons(indexSC_MSN(i),indexEndVelocity);
        
        if isempty(auxE{1,1}) || isempty(auxS{1,1})
            continue
        end

        for trial=1:length(auxS{1,1}) % Locomotion onset
            spikeCounts = histcounts(cell2mat(auxS{1,1}(trial)), edges7s);
            firingRate = spikeCounts / binSize;
            spksScS_MSN_pos = [spksScS_MSN_pos; firingRate];
        end

        
        for trial=1:length(auxE{1,1}) % Locomotion offset
            spikeCounts = histcounts(cell2mat(auxE{1,1}(trial)), edges7s);
            firingRate = spikeCounts / binSize;
            spksScE_MSN_pos = [spksScE_MSN_pos; firingRate];
        end
        
        allN_spksScS_MSN_pos = [allN_spksScS_MSN_pos; mean(spksScS_MSN_pos)/mean(mean([spksScS_MSN_pos; spksScE_MSN_pos]))];
        allN_spksScE_MSN_pos = [allN_spksScE_MSN_pos; mean(spksScE_MSN_pos)/mean(mean([spksScS_MSN_pos; spksScE_MSN_pos]))];
        
        spksScS_MSN_pos = [];
        spksScE_MSN_pos = [];
        
    else % Negative MSN
        
        auxS = spkNeurons(indexSC_MSN(i),5);
        auxE = spkNeurons(indexSC_MSN(i),indexEndVelocity);
        

        if isempty(auxE{1,1}) || isempty(auxS{1,1})
            continue
        end

        for trial=1:length(auxS{1,1}) % Locomotion onset
            spikeCounts = histcounts(cell2mat(auxS{1,1}(trial)), edges7s);
            firingRate = spikeCounts / binSize;
            spksScS_MSN_neg = [spksScS_MSN_neg;firingRate];
        end

        for trial=1:length(auxE{1,1}) % Locomotion offset
            spikeCounts = histcounts(cell2mat(auxE{1,1}(trial)), edges7s);
            firingRate = spikeCounts / binSize;
            spksScE_MSN_neg = [spksScE_MSN_neg;firingRate];
        end

        allN_spksScS_MSN_neg = [allN_spksScS_MSN_neg; mean(spksScS_MSN_neg)/mean(mean([spksScS_MSN_neg; spksScE_MSN_neg]))];
        allN_spksScE_MSN_neg = [allN_spksScE_MSN_neg; mean(spksScE_MSN_neg)/mean(mean([spksScS_MSN_neg; spksScE_MSN_neg]))];
        
        spksScS_MSN_neg = [];
        spksScE_MSN_neg = [];
    end

end


for i=1:length(indexSC_FSI)
    if valorCcg(indexSC_FSI(i)) > 0 % Positive FSI
        auxS = spkNeurons(indexSC_FSI(i),5);
        auxE = spkNeurons(indexSC_FSI(i),indexEndVelocity);
                       
        if isempty(auxE{1,1}) || isempty(auxS{1,1})
            continue
        end
        for trial=1:length(auxS{1,1}) % Locomotion onset
            spikeCounts = histcounts(cell2mat(auxS{1,1}(trial)), edges7s);
            firingRate = spikeCounts / binSize;
            spksScS_FSI_pos = [spksScS_FSI_pos;firingRate];
        end

        for trial=1:length(auxE{1,1}) % Locomotion offset
            spikeCounts = histcounts(cell2mat(auxE{1,1}(trial)), edges7s);
            firingRate = spikeCounts / binSize;
            spksScE_FSI_pos = [spksScE_FSI_pos;firingRate];
        end
        
        allN_spksScS_FSI_pos = [allN_spksScS_FSI_pos; mean(spksScS_FSI_pos)/mean(mean([spksScS_FSI_pos; spksScE_FSI_pos]))];
        allN_spksScE_FSI_pos = [allN_spksScE_FSI_pos; mean(spksScE_FSI_pos)/mean(mean([spksScS_FSI_pos; spksScE_FSI_pos]))];
        
        spksScS_FSI_pos = [];
        spksScE_FSI_pos = [];
        
    else % Negative FSI
        auxS = spkNeurons(indexSC_FSI(i),5);
        auxE = spkNeurons(indexSC_FSI(i),indexEndVelocity);

        if isempty(auxE{1,1}) || isempty(auxS{1,1})
            continue
        end
        for trial=1:length(auxS{1,1}) % Locomotion onset
            spikeCounts = histcounts(cell2mat(auxS{1,1}(trial)), edges7s);
            firingRate = spikeCounts / binSize;
            spksScS_FSI_neg = [spksScS_FSI_neg;firingRate];
        end

        for trial=1:length(auxE{1,1}) % Locomotion offset
            spikeCounts = histcounts(cell2mat(auxE{1,1}(trial)), edges7s);
            firingRate = spikeCounts / binSize;
            spksScE_FSI_neg = [spksScE_FSI_neg;firingRate];
        end

        allN_spksScS_FSI_neg = [allN_spksScS_FSI_neg; mean(spksScS_FSI_neg)/mean(mean([spksScS_FSI_neg; spksScE_FSI_neg]))];
        allN_spksScE_FSI_neg = [allN_spksScE_FSI_neg; mean(spksScE_FSI_neg)/mean(mean([spksScS_FSI_neg; spksScE_FSI_neg]))];
        
        spksScS_FSI_neg = [];
        spksScE_FSI_neg = [];
    end
end
%% 
clf, clc
fig1 = gcf;
lines = 2;
columns = 2;
set(gcf,'color','w')
set(fig1,'position',[0, 0, 1400, 1050]);
pos1 = [0.13,0.489517819706499,0.775,0.435482180293501];
yMaxSpks = 1.8;
yMinSpks = .6;
yMaxVel = 50;
yticksVel = [0:5:yMaxVel];
xLim = [0 max(edges7s)];
smoothFactor = 10;
colorAllNeurons = '#37a259';

figWidth = 0.35;   % widht
figHeight = 0.4;   % hight
hSpacing = 0.12;   % W-space 
vSpacing = 0.08;   % H-space

x1 = 0.1;
y1 = 0.55;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1 = axes('Position', [x1, y1, figWidth, figHeight]);
title('Positive speed cells', 'FontWeight', 'bold','FontSize', 18,'color','k');
hold on;
yyaxis left;
plot(edges7s(1:end-1),smooth(mean(vertcat(allN_spksScS_FSI_pos,allN_spksScS_MSN_pos)),smoothFactor),'color',colorAllNeurons,'lineW',2, 'lineS','-')
plot(edges7s(1:end-1),smooth(mean(allN_spksScS_FSI_pos),smoothFactor),'color',[.988 .263 .4],'lineW',2, 'lineS','-')
plot(edges7s(1:end-1),smooth(mean(allN_spksScS_MSN_pos),smoothFactor),'color',[.682 .698 1],'lineW',2, 'lineS','-')

A = smooth(mean(vertcat(allN_spksScS_FSI_pos,allN_spksScS_MSN_pos)),smoothFactor)+...
    std(mean(vertcat(allN_spksScS_FSI_pos,allN_spksScS_MSN_pos)),1)...
    /sqrt(size(vertcat(allN_spksScS_FSI_pos,allN_spksScS_MSN_pos),1));
B = smooth(mean(vertcat(allN_spksScS_FSI_pos,allN_spksScS_MSN_pos)),smoothFactor)-...
    std(mean(vertcat(allN_spksScS_FSI_pos,allN_spksScS_MSN_pos)),1)...
    /sqrt(size(vertcat(allN_spksScS_FSI_pos,allN_spksScS_MSN_pos),1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.216 .635 .349],'facealpha',0.5,'EdgeColor','none', 'HandleVisibility', 'off')
 

A = smooth(mean(allN_spksScS_FSI_pos),smoothFactor)+...
    std(mean(allN_spksScS_FSI_pos),1)...
    /sqrt(size(allN_spksScS_FSI_pos,1));
B = smooth(mean(allN_spksScS_FSI_pos),smoothFactor)-...
    std(mean(allN_spksScS_FSI_pos),1)...
    /sqrt(size(allN_spksScS_FSI_pos,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.988 .263 .4],'facealpha',0.2,'EdgeColor','none', 'HandleVisibility', 'off')
 
A = smooth(mean(allN_spksScS_MSN_pos),smoothFactor)+...
    std(mean(allN_spksScS_MSN_pos),1)...
    /sqrt(size(allN_spksScS_MSN_pos,1));
B = smooth(mean(allN_spksScS_MSN_pos),smoothFactor)-...
    std(mean(allN_spksScS_MSN_pos),1)...
    /sqrt(size(allN_spksScS_MSN_pos,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.682 .698 1],'facealpha',0.2,'EdgeColor','none', 'HandleVisibility', 'off')
 

xlim(xLim)
ylim([yMinSpks yMaxSpks])
xticks(0:1:7)
xticklabels([])
yticks(yMinSpks:.2:yMaxSpks)
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
ylabel({'\fontsize{18}Locomotion onset', '\fontsize{13}Normalized firing rate'}, 'FontWeight', 'bold');

yyaxis right;
x = linspace(0, max(edges7s), length(meanVelocityStarting));
plot(x, smooth(meanVelocityStarting, smoothFactor),'k','lineW',2); 
ylim([0 35]);
h = ylabel('\fontsize{13}Mean speed (cm/s)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [0.4, 0, 0]); 
xt = xticks; 
for i = 2:length(xt)-1
    xline(xt(i), '--', 'color',[.8 .8 .8 .2], 'HandleVisibility', 'off'); 
end

A = smooth(mean(matrizS, 2),smoothFactor)+...
    std(mean(matrizS, 2),1)...
    /sqrt(size(sumStart,1));
B = smooth(mean(matrizS, 2),smoothFactor)-...
    std(mean(matrizS, 2),1)...
    /sqrt(size(sumStart,1));
fill([x fliplr(x)], ...
     [A' fliplr(B')],...
     [.51 .51 .51],'facealpha',0.2,'EdgeColor','none', 'HandleVisibility', 'off')

lgd = legend({'All','FSI','MSN','Speed'});

set(legend,...
    'Position',[0.072 0.89 0.16428571156093 0.0526315775268009],...
    'NumColumns',1,...
    'FontWeight','bold',...
    'FontSize',13);
legend boxoff
text(-.15, 1.05, 'A', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax3 = axes('Position', [x1, y1 - figHeight - vSpacing, figWidth, figHeight]);
hold on
plot(edges7s(1:end-1),smooth(mean(vertcat(allN_spksScE_FSI_pos,allN_spksScE_MSN_pos)),smoothFactor),'color',colorAllNeurons,'lineW',2, 'lineS','-')
plot(edges7s(1:end-1),smooth(mean(allN_spksScE_FSI_pos),smoothFactor),'color',[.988 .263 .4],'lineW',2, 'lineS','-')
plot(edges7s(1:end-1),smooth(mean(allN_spksScE_MSN_pos),smoothFactor),'color',[.682 .698 1],'lineW',2, 'lineS','-')

A = smooth(mean(allN_spksScE_FSI_pos),smoothFactor)+...
    std(mean(allN_spksScE_FSI_pos),1)...
    /sqrt(size(allN_spksScE_FSI_pos,1));
B = smooth(mean(allN_spksScE_FSI_pos),smoothFactor)-...
    std(mean(allN_spksScE_FSI_pos),1)...
    /sqrt(size(allN_spksScE_FSI_pos,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.988 .263 .4],'facealpha',0.2,'EdgeColor','none', 'HandleVisibility', 'off')
 
A = smooth(mean(allN_spksScE_MSN_pos),smoothFactor)+...
    std(mean(allN_spksScE_MSN_pos),1)...
    /sqrt(size(allN_spksScE_MSN_pos,1));
B = smooth(mean(allN_spksScE_MSN_pos),smoothFactor)-...
    std(mean(allN_spksScE_MSN_pos),1)...
    /sqrt(size(allN_spksScE_MSN_pos,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.682 .698 1],'facealpha',0.2,'EdgeColor','none', 'HandleVisibility', 'off')
 
A = smooth(mean(vertcat(allN_spksScE_FSI_pos,allN_spksScE_MSN_pos)),smoothFactor)+...
    std(mean(vertcat(allN_spksScE_FSI_pos,allN_spksScE_MSN_pos)),1)...
    /sqrt(size(vertcat(allN_spksScE_FSI_pos,allN_spksScE_MSN_pos),1));
B = smooth(mean(vertcat(allN_spksScE_FSI_pos,allN_spksScE_MSN_pos)),smoothFactor)-...
    std(mean(vertcat(allN_spksScE_FSI_pos,allN_spksScE_MSN_pos)),1)...
    /sqrt(size(vertcat(allN_spksScE_FSI_pos,allN_spksScE_MSN_pos),1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.216 .635 .349],'facealpha',0.5,'EdgeColor','none', 'HandleVisibility', 'off')


hold on;
yyaxis left;

xlim(xLim)
ylim([yMinSpks yMaxSpks])
xticks(0:1:7)
yticks(yMinSpks:.2:yMaxSpks)
xlabel('Time (s)')
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
ylabel({'\fontsize{18}Locomotion offset', '\fontsize{13}Normalized firing rate'}, 'FontWeight', 'bold');

yyaxis right;
x = linspace(0, max(edges7s), length(meanVelocityEnding));
plot(x, smooth(meanVelocityEnding, smoothFactor),'k','lineW',2); 
ylim([0 35]);
h = ylabel('\fontsize{13}Mean speed (cm/s)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [0.4, 0, 0]); 
xt = xticks; 
for i = 2:length(xt)-1
    xline(xt(i), '--', 'color',[.8 .8 .8 .2], 'HandleVisibility', 'off'); 
end

A = smooth(mean(matrizE, 2),smoothFactor)+...
    std(mean(matrizE, 2),1)...
    /sqrt(size(sumEndind,1));
B = smooth(mean(matrizE, 2),smoothFactor)-...
    std(mean(matrizE, 2),1)...
    /sqrt(size(sumEndind,1));
fill([x fliplr(x)], ...
     [A' fliplr(B')],...
     [.51 .51 .51],'facealpha',0.2,'EdgeColor','none', 'HandleVisibility', 'off')

text(-.16, 1.05, 'B', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax2 = axes('Position', [x1 + figWidth + hSpacing, y1, figWidth, figHeight]);

hold on;
title('Negative speed cells', 'FontWeight', 'bold','FontSize', 18,'color','k');
yyaxis left;
plot(edges7s(1:end-1),smooth(mean(vertcat(allN_spksScS_FSI_neg,allN_spksScS_MSN_neg)),smoothFactor),'color',colorAllNeurons,'lineW',2, 'lineS','-')
plot(edges7s(1:end-1),smooth(mean(allN_spksScS_FSI_neg),smoothFactor),'color',[.988 .263 .4],'lineW',2, 'lineS','-')
plot(edges7s(1:end-1),smooth(mean(allN_spksScS_MSN_neg),smoothFactor),'color',[.682 .698 1],'lineW',2, 'lineS','-')


A = smooth(mean(vertcat(allN_spksScS_FSI_neg,allN_spksScS_MSN_neg)),smoothFactor)+...
    std(mean(vertcat(allN_spksScS_FSI_neg,allN_spksScS_MSN_neg)),1)...
    /sqrt(size(vertcat(allN_spksScS_FSI_neg,allN_spksScS_MSN_neg),1));
B = smooth(mean(vertcat(allN_spksScS_FSI_neg,allN_spksScS_MSN_neg)),smoothFactor)-...
    std(mean(vertcat(allN_spksScS_FSI_neg,allN_spksScS_MSN_neg)),1)...
    /sqrt(size(vertcat(allN_spksScS_FSI_neg,allN_spksScS_MSN_neg),1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.216 .635 .349],'facealpha',0.5,'EdgeColor','none', 'HandleVisibility', 'off')
 

A = smooth(mean(allN_spksScS_FSI_neg),smoothFactor)+...
    std(mean(allN_spksScS_FSI_neg),1)...
    /sqrt(size(allN_spksScS_FSI_neg,1));
B = smooth(mean(allN_spksScS_FSI_neg),smoothFactor)-...
    std(mean(allN_spksScS_FSI_neg),1)...
    /sqrt(size(allN_spksScS_FSI_neg,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.988 .263 .4],'facealpha',0.2,'EdgeColor','none', 'HandleVisibility', 'off')
 
A = smooth(mean(allN_spksScS_MSN_neg),smoothFactor)+...
    std(mean(allN_spksScS_MSN_neg),1)...
    /sqrt(size(allN_spksScS_MSN_neg,1));
B = smooth(mean(allN_spksScS_MSN_neg),smoothFactor)-...
    std(mean(allN_spksScS_MSN_neg),1)...
    /sqrt(size(allN_spksScS_MSN_neg,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.682 .698 1],'facealpha',0.2,'EdgeColor','none', 'HandleVisibility', 'off')
 
xlim(xLim)
ylim([yMinSpks yMaxSpks])
xticks(0:1:7)
xticklabels([])
yticks(yMinSpks:.2:yMaxSpks)
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
ylabel('\fontsize{13}Normalized firing rate', 'FontWeight', 'bold');

yyaxis right;
x = linspace(0, max(edges7s), length(meanVelocityStarting));
plot(x, smooth(meanVelocityStarting, smoothFactor),'k','lineW',2); 
ylim([0 35]);
h = ylabel('\fontsize{13}Mean speed (cm/s)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [0.4, 0, 0]); 
xt = xticks; 
for i = 2:length(xt)-1
    xline(xt(i), '--', 'color',[.8 .8 .8 .2], 'HandleVisibility', 'off'); 
end

A = smooth(mean(matrizS, 2),smoothFactor)+...
    std(mean(matrizS, 2),1)...
    /sqrt(size(sumStart,1));
B = smooth(mean(matrizS, 2),smoothFactor)-...
    std(mean(matrizS, 2),1)...
    /sqrt(size(sumStart,1));
fill([x fliplr(x)], ...
     [A' fliplr(B')],...
     [.51 .51 .51],'facealpha',0.2,'EdgeColor','none', 'HandleVisibility', 'off')

text(-.15, 1.05, 'C', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax4 = axes('Position', [x1 + figWidth + hSpacing, y1 - figHeight - vSpacing, figWidth, figHeight]);

hold on;
yyaxis left;
plot(edges7s(1:end-1),smooth(mean(allN_spksScE_FSI_neg),smoothFactor),'color',[.988 .263 .4],'lineW',2, 'lineS','-');
plot(edges7s(1:end-1),smooth(mean(allN_spksScE_MSN_neg),smoothFactor),'color',[.682 .698 1],'lineW',2, 'lineS','-');
plot(edges7s(1:end-1),smooth(mean(vertcat(allN_spksScE_FSI_neg,allN_spksScE_MSN_neg)),smoothFactor),'color',[.216 .635 .349],'lineW',2, 'lineS','-')

A = smooth(mean(allN_spksScE_FSI_neg),smoothFactor)+...
    std(mean(allN_spksScE_FSI_neg),1)...
    /sqrt(size(allN_spksScE_FSI_neg,1));
B = smooth(mean(allN_spksScE_FSI_neg),smoothFactor)-...
    std(mean(allN_spksScE_FSI_neg),1)...
    /sqrt(size(allN_spksScE_FSI_neg,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.988 .263 .4],'facealpha',0.2,'EdgeColor','none')
 
A = smooth(mean(allN_spksScE_MSN_neg),smoothFactor)+...
    std(mean(allN_spksScE_MSN_neg),1)...
    /sqrt(size(allN_spksScE_MSN_neg,1));
B = smooth(mean(allN_spksScE_MSN_neg),smoothFactor)-...
    std(mean(allN_spksScE_MSN_neg),1)...
    /sqrt(size(allN_spksScE_MSN_neg,1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.682 .698 1],'facealpha',0.2,'EdgeColor','none')

A = smooth(mean(vertcat(allN_spksScE_FSI_neg,allN_spksScE_MSN_neg)),smoothFactor)+...
    std(mean(vertcat(allN_spksScE_FSI_neg,allN_spksScE_MSN_neg)),1)...
    /sqrt(size(vertcat(allN_spksScE_FSI_neg,allN_spksScE_MSN_neg),1));
B = smooth(mean(vertcat(allN_spksScE_FSI_neg,allN_spksScE_MSN_neg)),smoothFactor)-...
    std(mean(vertcat(allN_spksScE_FSI_neg,allN_spksScE_MSN_neg)),1)...
    /sqrt(size(vertcat(allN_spksScE_FSI_neg,allN_spksScE_MSN_neg),1));
fill([edges7s(1:end-1) fliplr(edges7s(1:end-1))], ...
     [A' fliplr(B')],...
     [.216 .635 .349],'facealpha',0.5,'EdgeColor','none')
 
xlim(xLim)
ylim([yMinSpks yMaxSpks])
xticks(0:1:7)
yticks(yMinSpks:.2:yMaxSpks)
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13,'color','k');
ylabel('Normalized firing rate', 'FontWeight', 'bold','FontSize', 13);

yyaxis right;
x = linspace(0, max(edges7s), length(meanVelocityEnding));
plot(x, smooth(meanVelocityEnding, smoothFactor),'k','lineW',2);
ylim([0 35])
h = ylabel('\fontsize{13}Mean speed (cm/s)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [0.4, 0, 0]); 
xlabel('Time (s)','FontWeight', 'bold','FontSize', 12);

xt = xticks; 
for i = 2:length(xt)-1
    xline(xt(i), '--', 'color',[.8 .8 .8 .2]); 
end
A = smooth(mean(matrizE, 2),smoothFactor)+...
    std(mean(matrizE, 2),1)...
    /sqrt(size(sumEndind,1));
B = smooth(mean(matrizE, 2),smoothFactor)-...
    std(mean(matrizE, 2),1)...
    /sqrt(size(sumEndind,1));
fill([x fliplr(x)], ...
     [A' fliplr(B')],...
     [.51 .51 .51],'facealpha',0.2,'EdgeColor','none', 'HandleVisibility', 'off')
 
text(-.16, 1.05, 'D', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);
f = gcf;
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
