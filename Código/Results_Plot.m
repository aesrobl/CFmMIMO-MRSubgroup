%% Clear Data

%close all;
clear;

%% Plot simulation results
%Initial Data
nbrOfSetups = 500;
Scenario = '1x100';
Antennas = 8;
Subgroup_Cases = [1, 10, 25, 50, 75, 100];
% Cluster_Cases = [100, 50, 10, 4, 1];
% Num_K_Cluster_Cases = [1, 2, 10, 25, 100];
Color_Vector = ['k','r','g','b','y','c'];
K = 100;
%Visualize the data
figure(1);
hold on; box on;
set(gca,'fontsize',16);

% For MR Cases
for i = 1:length(Subgroup_Cases)
    load(sprintf('100x%d-%s-MR-multi%d.mat', Antennas, Scenario, Subgroup_Cases(i)));
    plot(sort(SE_MR_multi(:)),linspace(0,1,nbrOfSetups*K), [Color_Vector(i) '-'],'LineWidth',2);
end

hold on;

% For Normalized MR Cases
for i = 1:length(Subgroup_Cases)
    load(sprintf('100x%d-%s-MR-normalized-multi%d.mat', Antennas, Scenario, Subgroup_Cases(i)));
    plot(sort(SE_MR_normalized_multi(:)),linspace(0,1,nbrOfSetups*K), [Color_Vector(i) '--'],'LineWidth',2);
end

hold on;

% For Enhanced MR Cases
for i = 1:length(Subgroup_Cases)
    load(sprintf('100x%d-%s-MR-enhanced-multi%d.mat', Antennas, Scenario, Subgroup_Cases(i)));
    plot(sort(SE_MR_enhanced_multi(:)),linspace(0,1,nbrOfSetups*K), [Color_Vector(i) ':'],'LineWidth',2);
end

xlabel('Eficiencia Espectral Acumulada [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'G = 1 (CB)', 'G = 10 (CB)', 'G = 25 (CB)', 'G = 50 (CB)', 'G = 75 (CB)', 'G = 100 (CB)','G = 1 (NCB)', 'G = 10 (NCB)', 'G = 25 (NCB)', 'G = 50 (NCB)', 'G = 75 (NCB)', 'G = 100 (NCB)','G = 1 (ECB)', 'G = 10 (ECB)', 'G = 25 (ECB)', 'G = 50 (ECB)', 'G = 75 (ECB)', 'G = 100 (ECB)'},'Interpreter','Latex','Location','SouthEast');
%title("CDF ASE for Clustered Cell-Free mMIMO (100x1) CB vs. NCB vs. ECB - 8 antennas");

% legend({'G = 1', 'G = 10', 'G = 25', 'G = 50', 'G = 75', 'G = 100','G = 1', 'G = 10', 'G = 25', 'G = 50', 'G = 75', 'G = 100','G = 1', 'G = 10', 'G = 25', 'G = 50', 'G = 75', 'G = 100'},'Interpreter','Latex','Location','SouthEast');
