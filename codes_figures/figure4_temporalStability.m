clc, clf, clear
infoSession = 'onlyInTrials';
% infoSession = 'wholeSession';

timeInfoSession = '1sInterval';
% timeInfoSession = '0.5sInterval';
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
strCells = allCells(cell2mat(allCells(:,6)) >=21,:);

load('.\LFPDataMatlabFormat\spikes\spikes\msnFsiAndLocomotionCellsInfo\ccgLagFrLocomotionCellsData_OC_str_surNumber_10000.mat')
locomotioncell = ccgLagFrLocomotionCellsData{4};

%%
clc
maxLag = 10;
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

    parte = floor(maxTrial / numPortion);
    resto = mod(maxTrial, numPortion);
    matriz_vetores = cell(numPortion, 1);
    
    inicio = 1;
    for i = 1:numPortion
        if i <= resto
            fim = inicio + parte;
        else
            fim = inicio + parte - 1;
        end
        matriz_vetores{i} = inicio:fim;
        inicio = fim + 1;
    end
    
    for portion = 1:numPortion
        dataFR = data{4}(ismember(data{4}(:,6), matriz_vetores{portion}),1);
        dataVel = data{5}(ismember(data{4}(:,6), matriz_vetores{portion}),1);
        tempo_firing_rate = 1:1:length(dataFR);
%         smooth_fr = suavizarDisparos(tempo_firing_rate, dataFR, std_gaussiana);
%         smooth_vel = suavizarDisparos(tempo_firing_rate, dataVel, std_gaussiana);
        smooth_fr = dataFR;
        smooth_vel = dataVel;

        [ccg, lags] = xcorr(zscore(smooth_fr), zscore(smooth_vel), maxLag, 'coeff');
        [m i] = max(ccg);
        valorCcg_portion{neuron,portion} = m;
%         valoresLag_portion{neuron,portion} = lags(i);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Correct choice
    I = find(data{4}(:,3) == 0); % Wrong
    dataFR = data{4}(I);
    dataVel = data{5}(I);
    tempo_firing_rate = 1:1:length(dataFR);
%     smooth_fr = suavizarDisparos(tempo_firing_rate, dataFR, std_gaussiana);
%     smooth_vel = suavizarDisparos(tempo_firing_rate, dataVel, std_gaussiana);
    smooth_fr = dataFR;
    smooth_vel = dataVel;

    [ccg, lags] = xcorr(zscore(smooth_fr), zscore(smooth_vel), maxLag, 'coeff');
    [m i] = max(ccg);
    valorCcg_Rew{neuron,1} = m;
%     valoresLag_Rew{neuron,1} = lags(i);
    
    I = find(data{4}(:,3) == 1); % Right
    dataFR = data{4}(I);
    dataVel = data{5}(I);
    tempo_firing_rate = 1:1:length(dataFR);
%     smooth_fr = suavizarDisparos(tempo_firing_rate, dataFR, std_gaussiana);
%     smooth_vel = suavizarDisparos(tempo_firing_rate, dataVel, std_gaussiana);
    smooth_fr = dataFR;
    smooth_vel = dataVel;

    [ccg, lags] = xcorr(zscore(smooth_fr), zscore(smooth_vel), maxLag, 'coeff');
    [m i] = max(ccg);
    valorCcg_Rew{neuron,2} = m;
%     valoresLag_Rew{neuron,2} = lags(i);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Visual and Spatial tasks
    I = find(data{4}(:,4) == 0);    % Visual
    dataFR = data{4}(I);
    dataVel = data{5}(I);
    tempo_firing_rate = 1:1:length(dataFR);
%     smooth_fr = suavizarDisparos(tempo_firing_rate, dataFR, std_gaussiana);
%     smooth_vel = suavizarDisparos(tempo_firing_rate, dataVel, std_gaussiana);
    smooth_fr = dataFR;
    smooth_vel = dataVel;

    [ccg, lags] = xcorr(zscore(smooth_fr), zscore(smooth_vel), maxLag, 'coeff');
    [m i] = max(ccg);
    valorCcg_taskTrials{neuron,1} = m;
%     valoresLag_taskTrials{neuron,1} = lags(i);
    
    I = find(data{4}(:,4) == 1);    % Spatial
    dataFR = data{4}(I);
    dataVel = data{5}(I);
    tempo_firing_rate = 1:1:length(dataFR);
%     smooth_fr = suavizarDisparos(tempo_firing_rate, dataFR, std_gaussiana);
%     smooth_vel = suavizarDisparos(tempo_firing_rate, dataVel, std_gaussiana);
    smooth_fr = dataFR;
    smooth_vel = dataVel;

    [ccg, lags] = xcorr(zscore(smooth_fr), zscore(smooth_vel), maxLag, 'coeff');
    [m i] = max(ccg);
    valorCcg_taskTrials{neuron,2} = m;
%     valoresLag_taskTrials{neuron,2} = lags(i);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Maze side
    I = find(data{4}(:,5) == 0); % left
    dataFR = data{4}(I);
    dataVel = data{5}(I);
    tempo_firing_rate = 1:1:length(dataFR);
%     smooth_fr = suavizarDisparos(tempo_firing_rate, dataFR, std_gaussiana);
%     smooth_vel = suavizarDisparos(tempo_firing_rate, dataVel, std_gaussiana);
    smooth_fr = dataFR;
    smooth_vel = dataVel;

    [ccg, lags] = xcorr(zscore(smooth_fr), zscore(smooth_vel), maxLag, 'coeff');
    [m i] = max(ccg);
    valorCcg_mazeSide{neuron,1} = m;
%     valoresLag_mazeSide{neuron,1} = lags(i);
    
    I = find(data{4}(:,5) == 1); % right
    dataFR = data{4}(I);
    dataVel = data{5}(I);
    tempo_firing_rate = 1:1:length(dataFR);
%     smooth_fr = suavizarDisparos(tempo_firing_rate, dataFR, std_gaussiana);
%     smooth_vel = suavizarDisparos(tempo_firing_rate, dataVel, std_gaussiana);
    smooth_fr = dataFR;
    smooth_vel = dataVel;

    [ccg, lags] = xcorr(zscore(smooth_fr), zscore(smooth_vel), maxLag, 'coeff');
    [m i] = max(ccg);
    valorCcg_mazeSide{neuron,2} = m;
%     valoresLag_mazeSide{neuron,2} = lags(i);
    
    neuron
end
%%
clc, clf
fig1 = gcf;
set(gcf,'color','w')
set(fig1,'position',[0, 0, 1400, 1050]);
pos1 = [0.13,0.489517819706499,0.775,0.435482180293501];
linhas = 4;
colunas = 3;
pValueThreshold = .01;
colorCirclScatterFSI = '#FC4366';
colorScatterFSI = '#FC4366';

colorCirclScatterMSN = '#AEB2FF';
colorScatterMSN = '#AEB2FF';

FSI = (cell2mat(strCells(:,6)) == 22) & (locomotioncell==1)';
MSN = (cell2mat(strCells(:,6)) == 21) & (locomotioncell==1)';
% ccgLagFrLocomotionCellsData{1}(:,MSN)
% ccgLagFrLocomotionCellsData{1}(:,FSI)
subplot(linhas, colunas,[1 3])
hold on;
dataFSI = cell2mat(valorCcg_portion(FSI,:));
plot(1:5, dataFSI,'color',[.988 .263 .4 0.3], 'HandleVisibility', 'off')

plot([0 6.5], [mean(ccgLagFrLocomotionCellsData{1}(:,FSI)) mean(ccgLagFrLocomotionCellsData{1}(:,FSI))], 'color',[.988 .263 .4],'lineS','--', 'lineW',2, 'HandleVisibility', 'off')

[h p ci stats] = ttest(dataFSI, mean(ccgLagFrLocomotionCellsData{1}(:,FSI)),...
    'alpha',pValueThreshold/5);
p<pValueThreshold/5
errorbar(1:length(mean(cell2mat(valorCcg_portion(FSI,:)))),...
    mean(dataFSI),...
    (ci(2,1) - ci(1,:))/2,...
    'k', 'LineStyle', 'none','lineW',1.5, 'HandleVisibility', 'off');
scatter(1:length(mean(cell2mat(valorCcg_portion(FSI,:)))),mean(cell2mat(valorCcg_portion(FSI,:))',2), 'MarkerEdgeColor',colorScatterFSI,'MarkerFaceColor',colorCirclScatterFSI)

% text(1.98, .46,'*','FontSize', 21)
xlim([.8 5.2])
xticks(1:5)
xticklabels({'1st','2nd','3rd','4th','5th'})
ylabel('Speed score (r)')
xlabel('Trial block')
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13);
ylim([0, .7])
yticks(0:.1:.7)
legend({'FSI'},'FontWeight', 'bold','FontSize', 13,'NumColumns',2);
text(-.07, 1.1, 'A', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

subplot(linhas, colunas,4)
hold on;
dataFSI = cell2mat(valorCcg_taskTrials(FSI,:));
plot(1:2, dataFSI,'color',[.988 .263 .4 0.3])
plot([0 6.5], [mean(ccgLagFrLocomotionCellsData{1}(:,FSI)) mean(ccgLagFrLocomotionCellsData{1}(:,FSI))], 'color',[.988 .263 .4],'lineS','--', 'lineW',2)

[h p ci stats] = ttest(dataFSI, mean(ccgLagFrLocomotionCellsData{1}(:,FSI)),...
    'alpha',pValueThreshold/2);
p<pValueThreshold/2
errorbar(1:length(mean(cell2mat(valorCcg_taskTrials(FSI,:)))),...
    mean(dataFSI),...
    (ci(2,1) - ci(1,:))/2,...
    'k', 'LineStyle', 'none','lineW',1.5, 'HandleVisibility', 'off');
scatter(1:length(mean(cell2mat(valorCcg_taskTrials(FSI,:)))),mean(cell2mat(valorCcg_taskTrials(FSI,:))',2), 'MarkerEdgeColor',colorScatterFSI,'MarkerFaceColor',colorCirclScatterFSI)

xlim([.8 2.2])
xticks(1:2)
xticklabels({'Visual','Spatial'})
ylabel('Speed score (r)')
xlabel('Cue')
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13);
ylim([0, .6])
yticks(0:.1:.6)
text(-.26, 1.1, 'B', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

subplot(linhas, colunas,5)
hold on;
dataFSI = cell2mat(valorCcg_mazeSide(FSI,:));
plot(1:2, dataFSI,'color',[.988 .263 .4 0.3])
plot([0 6.5], [mean(ccgLagFrLocomotionCellsData{1}(:,FSI)) mean(ccgLagFrLocomotionCellsData{1}(:,FSI))], 'color',[.988 .263 .4],'lineS','--', 'lineW',2)

[h p ci stats] = ttest(dataFSI, mean(ccgLagFrLocomotionCellsData{1}(:,FSI)),...
    'alpha',pValueThreshold/2);
p<pValueThreshold/2

errorbar(1:length(mean(cell2mat(valorCcg_mazeSide(FSI,:)))),...
    mean(dataFSI),...
    (ci(2,1) - ci(1,:))/2,...
    'k', 'LineStyle', 'none','lineW',1.5, 'HandleVisibility', 'off');
scatter(1:length(mean(cell2mat(valorCcg_mazeSide(FSI,:)))),mean(cell2mat(valorCcg_mazeSide(FSI,:))',2), 'MarkerEdgeColor',colorScatterFSI,'MarkerFaceColor',colorCirclScatterFSI)
xlim([.8 2.2])
xticks(1:2)
xticklabels({'Left','Right'})
xlabel('Spatial choice')
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13);
ylim([0, .6])
yticks(0:.1:.6)
yticklabels([])
text(-.15, 1.1, 'C', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

subplot(linhas,colunas,6)
hold on;
dataFSI = cell2mat(valorCcg_Rew(FSI,:));
plot(1:2, dataFSI,'color',[.988 .263 .4 0.3])
plot([0 6.5], [mean(ccgLagFrLocomotionCellsData{1}(:,FSI)) mean(ccgLagFrLocomotionCellsData{1}(:,FSI))], 'color',[.988 .263 .4],'lineS','--', 'lineW',2)

[h p ci stats] = ttest(dataFSI, mean(ccgLagFrLocomotionCellsData{1}(:,FSI)),...
    'alpha',pValueThreshold/2);
p<pValueThreshold/2
errorbar(1:length(mean(cell2mat(valorCcg_Rew(FSI,:)))),...
    mean(dataFSI),...
    (ci(2,1) - ci(1,:))/2,...
    'k', 'LineStyle', 'none','lineW',1.5, 'HandleVisibility', 'off');

scatter(1:length(mean(cell2mat(valorCcg_Rew(FSI,:)))),mean(cell2mat(valorCcg_Rew(FSI,:))',2), 'MarkerEdgeColor',colorScatterFSI,'MarkerFaceColor',colorCirclScatterFSI)

xlim([.8 2.2])
xticks(1:2)
xticklabels({'NoReward','Reward'})
xlabel('Trial outcome')
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13);
ylim([0, .6])
yticks(0:.1:.6)
yticklabels([])
text(-.15, 1.1, 'D', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

subplot(linhas, colunas,[7 9])
hold on;
dataMSN = cell2mat(valorCcg_portion(MSN,:));
plot(1:5, dataMSN,'color',[.682 .698 1 0.3], 'HandleVisibility', 'off')
plot([0 6.5], [mean(ccgLagFrLocomotionCellsData{1}(:,MSN)) mean(ccgLagFrLocomotionCellsData{1}(:,MSN))], 'color',[.682 .698 1],'lineS','--', 'lineW',2, 'HandleVisibility', 'off')
[h p ci stats] = ttest(dataMSN, mean(ccgLagFrLocomotionCellsData{1}(:,MSN)),...
    'alpha',pValueThreshold/5);
p<pValueThreshold/5
errorbar(1:length(mean(cell2mat(valorCcg_portion(MSN,:)))),...
    mean(dataMSN),...
    (ci(2,1) - ci(1,:))/2,...
    'k', 'LineStyle', 'none','lineW',1.5, 'HandleVisibility', 'off');

scatter(1:length(mean(cell2mat(valorCcg_portion(MSN,:)))),mean(cell2mat(valorCcg_portion(MSN,:))',2), 'MarkerEdgeColor',colorScatterMSN,'MarkerFaceColor',colorCirclScatterMSN)

xlim([.8 5.2])
xticks(1:5)
xticklabels({'1st','2nd','3rd','4th','5th'})
ylabel('Speed score (r)')
xlabel('Trial block')
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13);
ylim([0, .7])
yticks(0:.1:.7)

text(1.98, .46,'*','FontSize', 21)
text(3.98, .46,'*','FontSize', 21)
legend('MSN','FontWeight', 'bold','FontSize', 13,'NumColumns',2);
text(-.07, 1.1, 'E', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);


subplot(linhas, colunas,10)
hold on;

dataMSN = cell2mat(valorCcg_taskTrials(MSN,:));
plot(1:2, dataMSN,'color',[.682 .698 1 0.3])

plot([0 6.5], [mean(ccgLagFrLocomotionCellsData{1}(:,MSN)) mean(ccgLagFrLocomotionCellsData{1}(:,MSN))], 'color',[.682 .698 1],'lineS','--', 'lineW',2)
[h p ci stats] = ttest(dataMSN, mean(ccgLagFrLocomotionCellsData{1}(:,MSN)),...
    'alpha',pValueThreshold/2);
p<pValueThreshold/2
errorbar(1:length(mean(cell2mat(valorCcg_taskTrials(MSN,:)))),...
    mean(dataMSN),...
    (ci(2,1) - ci(1,:))/2,...
    'k', 'LineStyle', 'none','lineW',1.5, 'HandleVisibility', 'off');

scatter(1:length(mean(cell2mat(valorCcg_taskTrials(MSN,:)))),mean(cell2mat(valorCcg_taskTrials(MSN,:))',2), 'MarkerEdgeColor',colorScatterMSN,'MarkerFaceColor',colorCirclScatterMSN)

% text(1.98, .46,'*','FontSize', 21)
xlim([.8 2.2])
xticks(1:2)
xticklabels({'Visual','Spatial'})
ylabel('Speed score (r)')
xlabel('Cue')
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13);
ylim([0, .6])
yticks(0:.1:.6)
text(-.26, 1.1, 'F', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);


subplot(linhas, colunas,11)
hold on;

dataMSN = cell2mat(valorCcg_mazeSide(MSN,:));
plot(1:2, dataMSN,'color',[.682 .698 1 0.3])

plot([0 6.5], [mean(ccgLagFrLocomotionCellsData{1}(:,MSN)) mean(ccgLagFrLocomotionCellsData{1}(:,MSN))], 'color',[.682 .698 1],'lineS','--', 'lineW',2)
[h p ci stats] = ttest(dataMSN, mean(ccgLagFrLocomotionCellsData{1}(:,MSN)),...
    'alpha',pValueThreshold/2);
p<pValueThreshold/2
errorbar(1:length(mean(cell2mat(valorCcg_mazeSide(MSN,:)))),...
    mean(dataMSN),...
    (ci(2,1) - ci(1,:))/2,...
    'k', 'LineStyle', 'none','lineW',1.5, 'HandleVisibility', 'off');
scatter(1:length(mean(cell2mat(valorCcg_mazeSide(MSN,:)))),mean(cell2mat(valorCcg_mazeSide(MSN,:))',2), 'MarkerEdgeColor',colorScatterMSN,'MarkerFaceColor',colorCirclScatterMSN)

xlim([.8 2.2])
xticks(1:2)
xticklabels({'Left','Right'})
xlabel('Spatial choice')
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13);
ylim([0, .6])
yticks(0:.1:.6)
yticklabels([])
text(-.15, 1.1, 'G', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

subplot(linhas,colunas,12)
hold on;
dataMSN = cell2mat(valorCcg_Rew(MSN,:));
plot(1:2, dataMSN,'color',[.682 .698 1 0.3])

plot([0 6.5], [mean(ccgLagFrLocomotionCellsData{1}(:,MSN)) mean(ccgLagFrLocomotionCellsData{1}(:,MSN))], 'color',[.682 .698 1],'lineS','--', 'lineW',2)
[h p ci stats] = ttest(dataMSN, mean(ccgLagFrLocomotionCellsData{1}(:,MSN)),...
    'alpha',pValueThreshold/2);
p<pValueThreshold/2
errorbar(1:length(mean(cell2mat(valorCcg_Rew(MSN,:)))),...
    mean(dataMSN),...
    (ci(2,1) - ci(1,:))/2,...
    'k', 'LineStyle', 'none','lineW',1.5, 'HandleVisibility', 'off');
scatter(1:length(mean(cell2mat(valorCcg_Rew(MSN,:)))),mean(cell2mat(valorCcg_Rew(MSN,:))',2), 'MarkerEdgeColor',colorScatterMSN,'MarkerFaceColor',colorCirclScatterMSN)

xlim([.8 2.2])
xticks(1:2)
xticklabels({'NoReward','Reward'})
xlabel('Trial outcome')
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 13);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 13);
ylim([0, .6])
yticks(0:.1:.6)
yticklabels([])
text(-.15, 1.1, 'H', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

f = gcf;

%% SAVE IMAGE
clc
name = strcat('figure4_TemporalAndCueStability_OC_str_',timeInfoSession);
saveName = strcat(name, '.tiff');
filename = fullfile('./codes/paperFigures_v1.0/figures/', saveName);
exportgraphics(f, filename, 'Resolution',300)

'done'
%%
%%
%%
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