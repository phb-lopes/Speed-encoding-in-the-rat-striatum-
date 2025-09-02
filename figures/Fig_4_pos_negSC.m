%%
clc, clf, clear
timeInfoSession = 1;
lag = 3;
overlap = 0.9;
folder = strcat('/Users/phb_lopes/Library/CloudStorage/OneDrive-Pessoal/projetoEletrofisiologiaMacbook/',num2str(timeInfoSession),'sInterval/', num2str(lag),'sLag/',num2str(overlap),'_overlap/'); % Change to whatever pattern you need
dataFsiAndLocomotionCells_all = readData(folder);
bin_step = timeInfoSession-overlap;
maxLag_seconds = lag;
cd /Users/phb_lopes/Library/CloudStorage/OneDrive-Pessoal/projetoEletrofisiologiaMacbook/ 
%%
allCells = {};
for i=1:length(dataFsiAndLocomotionCells_all)
    allCells = [allCells; dataFsiAndLocomotionCells_all{1,i}.data];
end
strCells = allCells(cell2mat(allCells(:,6)) >=21,:);

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
%%
clf, clc
fig1 = gcf;
set(gcf,'color','w')
set(fig1,'position',[0, 0, 1400, 1050]);
pos1 = [0.13,0.489517819706499,0.775,0.435482180293501];
linhas = 2;
colunas = 9;

FSI_pos_SC_ccg = valorCcg(valorCcg>0 & speedCell==1 & (cell2mat(strCells(:,6)) == 22)');
MSN_pos_SC_ccg = valorCcg(valorCcg>0 & speedCell==1 & (cell2mat(strCells(:,6)) == 21)');
FSI_pos_SC_lag_ms = valoresLag_ms(valorCcg>0 & speedCell==1 & (cell2mat(strCells(:,6)) == 22)');
MSN_pos_SC_lag_ms = valoresLag_ms(valorCcg>0 & speedCell==1 & (cell2mat(strCells(:,6)) == 21)');

[h,p] = ttest2(FSI_pos_SC_ccg, MSN_pos_SC_ccg)
[h,p] = ttest2(FSI_pos_SC_lag_ms, MSN_pos_SC_lag_ms)

subplot(linhas, colunas,[1 3])
hold on;
[values, bins] = hist(FSI_pos_SC_ccg, 0:.01:.6);
bar(bins, values,'faceC',[.988 .263 .4],'edgeC',[.988 .263 .4], 'barW', .8)
[values, bins] = hist(MSN_pos_SC_ccg, 0:.01:.6);
bar(bins,values*-1,'faceC',[.682 .698 1],'edgeC',[.682 .698 1], 'barW', .8)

% Weak
plot([0.205, 0.205], [-7 8], Color='k', lineStyle='--', lineW=2)
text(.095, 5, 'Weak', 'FontSize', 17, 'FontWeight', 'bold','Color', 'k')
w = sum(FSI_pos_SC_ccg<0.205);
text(.1, 4.2, strcat('n=',num2str(w)), 'FontSize', 16, 'FontWeight', 'bold','Color', [.988 .263 .4])
w = sum(MSN_pos_SC_ccg<0.205);
text(.1, -4.4, strcat('n=',num2str(w)), 'FontSize', 16, 'FontWeight', 'bold', 'Color', [.682 .698 1])

% Moderate
plot([0.405, 0.405], [-7 8], Color='k', lineStyle='--', lineW=2)
text(.25, 5, 'Moderate', 'FontSize', 17, 'FontWeight', 'bold','Color', 'k')
m = sum(FSI_pos_SC_ccg > 0.205 & FSI_pos_SC_ccg < 0.4);
text(.28, 4.2, strcat('n=',num2str(m)), 'FontSize', 16, 'FontWeight', 'bold','Color', [.988 .263 .4])
m = sum(MSN_pos_SC_ccg > 0.205 & MSN_pos_SC_ccg < 0.4);
text(.28, -4.3, strcat('n=',num2str(m)), 'FontSize', 16, 'FontWeight', 'bold', 'Color', [.682 .698 1])

% Strong
text(.45, 5, 'Strong', 'FontSize', 17, 'FontWeight', 'bold','Color', 'k')
s = sum(FSI_pos_SC_ccg > 0.4);
text(.47, 4.2, strcat('n=',num2str(s)), 'FontSize', 16, 'FontWeight', 'bold','Color', [.988 .263 .4])
s = sum(MSN_pos_SC_ccg > 0.4);
xlabel('Speed score (r)')
ylabel({'# Neurons'})
xlim([0 .5])
yticks(-5:1:5)
ylim([-5, 5])
yticklabels([flip(1:1:5),0:1:5])
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 15);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 15);
lgd = legend('FSI','MSN','FontWeight', 'bold','FontSize', 15);

title(lgd,'Speed cells')
set(legend,...
    'Position',[0.335 0.64 0.0621428563765117 0.0362462751204823],...
    'FontWeight','bold');
text(.2, 6.8,'Positive speed cells', 'FontWeight', 'bold', 'FontSize', 30);
text(-.1, 1.05, 'A', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 25);

axes1 = axes('Parent',fig1,...
    'Position',[0.4 0.583 0.1 0.345]);
hold(axes1,'on');
b = bar([1,10],[mean(FSI_pos_SC_ccg), mean(MSN_pos_SC_ccg)]);
rgbColors = [
    .988 .263 .4;
    .682 .698 1;
];
b.FaceColor = 'flat';
b.EdgeColor = 'flat';
b.CData = rgbColors;

hold on
errorbar([1],mean(FSI_pos_SC_ccg),std(FSI_pos_SC_ccg)/sqrt(length(FSI_pos_SC_ccg)),'k.','markersize',.1,'lineW',1.2)
errorbar([10],mean(MSN_pos_SC_ccg),std(MSN_pos_SC_ccg)/sqrt(length(MSN_pos_SC_ccg)),'k.','markersize',.1,'lineW',1.2)
xticks([1, 10])
xticklabels({'FSI','MSN'})
set(gca, 'YAxisLocation', 'right');
h = ylabel('\fontsize{16}Mean speed score (r)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [8, 0, 0]); % Ajuste a posição conforme necessário

set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 16);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 16);
box off
text(-.1, 1.05, 'B', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 25);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FSI_neg_SC_ccg = valorCcg(valorCcg<0 & speedCell==1 & (cell2mat(strCells(:,6)) == 22)');
MSN_neg_SC_ccg = valorCcg(valorCcg<0 & speedCell==1 & (cell2mat(strCells(:,6)) == 21)');

FSI_neg_SC_lag_ms = valoresLag_ms(valorCcg<0 & speedCell==1 & (cell2mat(strCells(:,6)) == 22)');
MSN_neg_SC_lag_ms = valoresLag_ms(valorCcg<0 & speedCell==1 & (cell2mat(strCells(:,6)) == 21)');

subplot(linhas, colunas,[6 8])
hold on;
[values, bins] = hist(FSI_neg_SC_ccg, -.6:.01:0);
bar(bins, values,'faceC',[.988 .263 .4],'edgeC',[.988 .263 .4], 'barW', .8)
[values, bins] = hist(MSN_neg_SC_ccg, -.6:.01:0);
bar(bins,values*-1,'faceC',[.682 .698 1],'edgeC',[.682 .698 1], 'barW', .8)

% Weak
plot([-0.205, -0.205], [-7 8], Color='k', lineStyle='--', lineW=2)
text(-.13, 5, 'Weak', 'FontSize', 17, 'FontWeight', 'bold','Color', 'k')
w = sum(FSI_neg_SC_ccg<0 & FSI_neg_SC_ccg>-0.205);
text(-.12, 4.2, strcat('n=',num2str(w)), 'FontSize', 16, 'FontWeight', 'bold','Color', [.988 .263 .4])
w = sum(MSN_neg_SC_ccg<0 & MSN_neg_SC_ccg>-0.205);
text(-.12, -4.4, strcat('n=',num2str(w)), 'FontSize', 16, 'FontWeight', 'bold', 'Color', [.682 .698 1])

% Moderate
plot([-0.405, -0.405], [-7 8], Color='k', lineStyle='--', lineW=2)
text(-.35, 5, 'Moderate', 'FontSize', 17, 'FontWeight', 'bold','Color', 'k')
m = sum(FSI_neg_SC_ccg < -0.205 & FSI_neg_SC_ccg > -0.4);
text(-.32, 4.2, strcat('n=',num2str(m)), 'FontSize', 16, 'FontWeight', 'bold','Color', [.988 .263 .4])
m = sum(MSN_neg_SC_ccg < -0.205 & MSN_neg_SC_ccg > -0.4);
text(-.32, -4.3, strcat('n=',num2str(m)), 'FontSize', 16, 'FontWeight', 'bold', 'Color', [.682 .698 1])

xlabel('Speed score (r)')
ylabel('# Neurons')
xlim([-.5 0])
yticks(-5:1:5)
ylim([-5, 5])
yticklabels([flip(1:1:5),0:1:5])
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 15);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 15);

text(-.32, 6.8,'Negative speed cells', 'FontWeight', 'bold', 'FontSize', 30);
text(-.1, 1.05, 'E', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 25);

axes1 = axes('Parent',fig1,...
    'Position',[0.85 0.583 0.1 0.345]);
hold(axes1,'on');
b = bar([1,10],[mean(FSI_neg_SC_ccg), mean(MSN_neg_SC_ccg)]);
rgbColors = [
    .988 .263 .4;
    .682 .698 1;
];
b.FaceColor = 'flat';
b.EdgeColor = 'flat';
b.CData = rgbColors;

hold on
errorbar([1],mean(FSI_neg_SC_ccg),std(FSI_neg_SC_ccg)/sqrt(length(FSI_neg_SC_ccg)),'k.','markersize',.1,'lineW',1.2)
errorbar([10],mean(MSN_neg_SC_ccg),std(MSN_neg_SC_ccg)/sqrt(length(MSN_neg_SC_ccg)),'k.','markersize',.1,'lineW',1.2)
xticks([1, 10])
xticklabels({'FSI','MSN'})
set(gca, 'YAxisLocation', 'right');
h = ylabel('\fontsize{16}Mean speed score (r)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [8, 0, 0]); % Ajuste a posição conforme necessário

set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 16);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 16);
box off
text(-.1, 1.05, 'F', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 25);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(linhas, colunas,[10 12])
hold(axes1,'on');
hold on;
[values, bins] = hist(FSI_pos_SC_lag_ms,binHisto);
bar(bins, values,'faceC',[.988 .263 .4],'edgeC',[.988 .263 .4], 'barW',.8)
[values, bins] = hist(MSN_pos_SC_lag_ms,binHisto);
bar(bins,values*-1,'faceC',[.682 .698 1],'edgeC',[.682 .698 1], 'barW',.8)
xlabel('Speed lag (s)')
ylabel('# Neurons')
xticks(-maxLag_seconds:.5:maxLag_seconds)
xlim([-2.5 2.5])
yticks(-5:1:5)
ylim([-5, 5])
yticklabels([flip(0:1:5),1:1:5])
ax = gca;
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 16);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 16);
xticklabel = ax.XTickLabel;  % Get the tick labels
xtickangle(0);
% legend('FSI','MSN','FontWeight', 'bold','FontSize', 12);
text(-.13, 1.05, 'C', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 25);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes1 = axes('Parent',fig1,...
    'Position',[0.4 0.105 0.1 0.345]);

hold(axes1,'on');
b = bar([1,10],[mean(FSI_pos_SC_lag_ms), mean(MSN_pos_SC_lag_ms)]);
rgbColors = [
    .988 .263 .4;
    .682 .698 1;
];
b.FaceColor = 'flat';
b.EdgeColor = 'flat';
b.CData = rgbColors;
hold on
errorbar([1],mean(FSI_pos_SC_lag_ms),std(FSI_pos_SC_lag_ms)/sqrt(length(FSI_pos_SC_lag_ms)),'k.','markersize',.1,'lineW',1.2)
errorbar([10],mean(MSN_pos_SC_lag_ms),std(MSN_pos_SC_lag_ms)/sqrt(length(MSN_pos_SC_lag_ms)),'k.','markersize',.1,'lineW',1.2)
xticklabels({'FSI','MSN'})
set(gca, 'YAxisLocation', 'right');
h = ylabel('\fontsize{16}Mean speed lag (s)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [8, 0, 0]); % Ajuste a posição conforme necessário

set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 16);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 16);
box off

ylim([-2.5 2.5])
yticks(-2.5:.5:2.5)

text(-.1, 1.05, 'D', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 25);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FSI_neg_SC_ccg = valorCcg(valorCcg<0 & speedCell==1 & (cell2mat(strCells(:,6)) == 22)');
MSN_neg_SC_ccg = valorCcg(valorCcg<0 & speedCell==1 & (cell2mat(strCells(:,6)) == 21)');

FSI_neg_SC_lag_ms = valoresLag_ms(valorCcg<0 & speedCell==1 & (cell2mat(strCells(:,6)) == 22)');
MSN_neg_SC_lag_ms = valoresLag_ms(valorCcg<0 & speedCell==1 & (cell2mat(strCells(:,6)) == 21)');

[h,p] = ttest2(FSI_neg_SC_ccg, MSN_neg_SC_ccg)
[h,p] = ttest2(FSI_neg_SC_lag_ms, MSN_neg_SC_lag_ms)


subplot(linhas, colunas,[15 17])
hold(axes1,'on');
hold on;
[values, bins] = hist(FSI_neg_SC_lag_ms,binHisto);
bar(bins, values,'faceC',[.988 .263 .4],'edgeC',[.988 .263 .4], 'barW',.8)
[values, bins] = hist(MSN_neg_SC_lag_ms,binHisto);
bar(bins,values*-1,'faceC',[.682 .698 1],'edgeC',[.682 .698 1], 'barW',.8)
xlabel('Speed lag (s)')
ylabel('# Neurons')
xticks(-maxLag_seconds:.5:maxLag_seconds)
xlim([-2.5 2.5])
yticks(-5:1:5)
ylim([-5, 5])
yticklabels([flip(0:1:5),1:1:5])
ax = gca;
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 16);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 16);
xticklabel = ax.XTickLabel;  % Get the tick labels
xtickangle(0);
text(-.13, 1.05, 'G', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 25);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes1 = axes('Parent',fig1,...
    'Position',[0.85 0.105 0.1 0.345]);
hold(axes1,'on');
b = bar([1,10],[mean(FSI_neg_SC_lag_ms), mean(MSN_neg_SC_lag_ms)]);
rgbColors = [
    .988 .263 .4;
    .682 .698 1;
];
b.FaceColor = 'flat';
b.EdgeColor = 'flat';
b.CData = rgbColors;
hold on
errorbar([1],mean(FSI_neg_SC_lag_ms),std(FSI_neg_SC_lag_ms)/sqrt(length(FSI_neg_SC_lag_ms)),'k.','markersize',.1,'lineW',1.2)
errorbar([10],mean(MSN_neg_SC_lag_ms),std(MSN_neg_SC_lag_ms)/sqrt(length(MSN_neg_SC_lag_ms)),'k.','markersize',.1,'lineW',1.2)
xticklabels({'FSI','MSN'})
set(gca, 'YAxisLocation', 'right');
h = ylabel('\fontsize{16}Mean speed lag (s)', 'FontWeight', 'bold','Rotation', -90);
set(h, 'Position', get(h, 'Position') + [8, 0, 0]); % Ajuste a posição conforme necessário

set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 16);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 16);
box off
text(-.1, 1.05, 'H', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 25);
ylim([-2.5 2.5])
yticks(-2.5:.5:2.5)

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
