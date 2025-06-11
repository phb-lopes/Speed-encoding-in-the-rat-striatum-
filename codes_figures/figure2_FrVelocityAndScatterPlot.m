%%
clc, clf, clear
brainRegion = 'striatum\';
folder = 'G:\OneDrive\DocData\LFPDataMatlabFormat\spikes\spikes\msnFsiAndLocomotionCellsInfo\allCells\obertosClassification\1sInterval\onlyInTrials\';

dataFsiAndLocomotionCells_all = readData(folder);

cd G:\OneDrive\DocData

%%
clc
cells = {};
for i=1:length(dataFsiAndLocomotionCells_all)
    cells = [cells; dataFsiAndLocomotionCells_all{1,i}.data];
end

neuron = 114;%20;%114; % {'R268 20141016'}    {[7]}    {[10]}
data = cells(neuron,:)
dataFR =  data{4}(:,1);
dataVel = data{5};
%%
clc, clf
fig1 = gcf;
set(gcf,'color','w')
set(fig1,'position',[0, 0, 1400, 1050]);
pos1 = [0.13,0.489517819706499,0.775,0.435482180293501];
linhas = 5;
colunas = 3;
std_gaussiana = .5;
colorScatter = '#7575a3';
colorCirclScatter = '#7575a3';
tempo_firing_rate = 1:1:length(dataFR);
ini = 350;%min(tempo_firing_rate);
fim = 550;%length(tempo_firing_rate);

% smooth_fr = suavizarDisparos(tempo_firing_rate, dataFR, std_gaussiana);
% smooth_vel = suavizarDisparos(tempo_firing_rate, dataVel, std_gaussiana);
smooth_fr = dataFR;
smooth_vel = dataVel;
subplot(linhas,colunas,[1:3])
plot(tempo_firing_rate, smooth_vel, 'color','#e6ac00','lineW',3);
xticks([])
xlim([ini fim])
ylabel({'Locomotion','speed (cm/s)'});
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
box off
text(-.07, 1.2, 'A', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

subplot(linhas,colunas,[4:6])
plot(tempo_firing_rate, smooth_fr, 'color','#8300b4','lineW',3);
xticks([])
xlim([ini fim])
ylabel({'Firing rate (Hz)'});
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
box off
text(-.07, 1.2, 'B', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

subplot(linhas,colunas,[7:9])
hold on
plot(tempo_firing_rate, zscore(smooth_fr), 'color','#8300b4','lineW',3);
plot(tempo_firing_rate, zscore(smooth_vel), 'color','#e6ac00','lineW',3);

xlim([ini fim])
ylabel('Z-score');
xlabel('Time (s)');
xticks([ini:25:fim])
xticklabels([0:25:(fim-ini)])
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
text(-.07, 1.2, 'C', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

subplot(linhas,colunas,[[10 11],[13 14]])
hold on
scatter(smooth_vel, smooth_fr, 'MarkerEdgeColor',colorScatter,'MarkerFaceColor',colorCirclScatter)
scatter(smooth_vel, smooth_fr, 'MarkerEdgeColor',colorScatter,'MarkerFaceColor',colorCirclScatter)
[y_fit, R2_1] = getFit_polyfit(smooth_vel, smooth_fr', 1);
[r, p] = coeficiente_pearson(smooth_vel, smooth_fr);
fprintf('Valor de p: %.3e\n', p);  % Exibição em notação científica

text(50,50,strcat(['r = ' num2str(round(r,3))]), 'color','k','FontWeight', 'bold','FontSize', 16)
plot(smooth_vel, y_fit, 'k-', 'LineWidth', 3);



ylabel('Firing rate (Hz)');
xlabel({'Locomotion speed (cm/s)'});
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
yy = ylim;
yTicks = yticks;
text(-.11, 1.05, 'D', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

% subplot(linhas,colunas,[[11],[14]])
% hold on
% boxplot(smooth_fr, 'Symbol', '.', 'Colors',[0.459,0.459,0.639],'OutlierSize', 20)
% [values, bins] = hist(smooth_fr,200);
% plot(normalize(values,'range')*.3+1.2, bins,'color',colorCirclScatter)
% set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
% set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
% xlim([.7 1.5])
% yticklabels([])
% xticks([])
% ylim(yy)
% yticks(yTicks)
% box off
% text(-.16, .98, 'E', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

subplot(linhas,colunas,[[12],[15]])
surNumber = 10000;
maxLag = 10;
for surIdx=1:surNumber
    shiftAmount = randi([20,50],1)/100;
    circularShift = round(length(smooth_vel)*shiftAmount,0);
    [ccg, ~] = xcorr(zscore(circshift(smooth_fr',circularShift)), zscore(smooth_vel), maxLag, 'normalized');
    allSurCCG(:,surIdx) = ccg;
end
hold on;
%     [ccg, lags] = xcorr(smooth_fr, smooth_vel, maxLag, 'normalized');
[ccg, lags] = xcorr(zscore(smooth_fr), zscore(smooth_vel), maxLag, 'coeff');

[m i] = max(ccg);
Threshold = max(mean(allSurCCG, 2)+2*std(allSurCCG,1, 2));
m > Threshold;
plot(lags, ccg,'k','linew',3)
plot(lags, mean(allSurCCG, 2),'r-','linew',3)
% plot(lags, mean(allSurCCG, 2)+2*std(allSurCCG,1, 2),'color',colorScatter,'linew',3)
% plot(lags, mean(allSurCCG, 2)-2*std(allSurCCG,1, 2),'color',colorScatter,'linew',3)

A = mean(allSurCCG, 2)+2*std(allSurCCG,1, 2);
B = mean(allSurCCG, 2)-2*std(allSurCCG,1, 2);
fill([1:21,fliplr(1:21)]-11,[A' B(21:-1:1)'], 'r','facealpha',0.2,'EdgeColor','none', 'HandleVisibility', 'off')
% plot(lags, mean(allSurCCG, 2)+2*std(allSurCCG,1, 2),'color',colorScatter,'linew',3)
% plot(lags, mean(allSurCCG, 2)-2*std(allSurCCG,1, 2),'color',colorScatter,'linew',3)

plot([0 0], [-.2 .6], 'color','#9F9E9E','lineS','--','linew',2, 'HandleVisibility', 'off')
xlabel('Lag (s)')
ylabel('CCG (r)')
set(get(gca, 'YAxis'), 'FontWeight', 'bold','FontSize', 12);
set(get(gca, 'XAxis'), 'FontWeight', 'bold','FontSize', 12);
yticks([-.2:.1:.6])
legend('Real','Sur (± 2 SD)', 'FontWeight', 'bold', 'FontSize', 11)
set(legend,'Position',[0.81 0.38 0.1 0.01])

box off
text(-.2, 1.05, 'E', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

f = gcf;
%% SAVE IMAGE
clc
name = 'figure2_FRVelocityAndScatterPlot';
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

function [y_fit, R2] = getFit_polyfit(data_vel, data_fr, degree)
    % Verifica se as entradas são vetores coluna
    if isrow(data_vel)
        data_vel = data_vel'; % Transforma em coluna
    end
    if isrow(data_fr)
        data_fr = data_fr'; % Transforma em coluna
    end
    
    % Ajuste polinomial
    p = polyfit(data_vel, data_fr, degree);
    y_fit = polyval(p, data_vel);
    
    % Cálculo de SSR e SST
    SSR = sum((data_fr - y_fit).^2);
    SST = sum((data_fr - mean(data_fr)).^2);
    
    % Cálculo do R2
    R2 = 1 - (SSR / SST);
end

function [r, p] = coeficiente_pearson(data_vel, data_fr)
    % Verifica se as entradas possuem o mesmo tamanho
    if length(data_vel) ~= length(data_fr)
        error('As duas variáveis de entrada devem ter o mesmo tamanho');
    end
    
    % Calcula o coeficiente de correlação de Pearson e o valor de p
    [r, p] = corr(data_vel(:), data_fr(:), 'Type', 'Pearson');
end