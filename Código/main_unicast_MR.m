%This Matlab script can be used to reproduce Figures 6.3(b), 6.4, 6.5(b), and 6.6 in the monograph:
%
%Ozlem Tugfe Demir, Emil Bjornson and Luca Sanguinetti (2021),
%"Foundations of User-Centric Cell-Free Massive MIMO", 
%Foundations and Trends in Signal Processing: Vol. 14: No. 3-4,
%pp 162-472. DOI: 10.1561/2000000109
%
%This is version 1.0 (Last edited: 2021-01-31)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.

%Empty workspace and close figures
close all;
clear;



%% Define simulation setup

%Number of Monte-Carlo setups
nbrOfSetups = 250;

%Number of channel realizations per setup
nbrOfRealizations = 500;


%Number of APs 
L = 100;

%Number of antennas per AP
N = 4;

%Number of UEs in the network per cluster
Kc = 50;
Cluster = 10;
K = Kc*ones(1,Cluster);

%Length of coherence block
tau_c = 200;

%Length of pilot sequences
tau_p = 20;

%Angular standard deviation in the local scattering model (in radians)
ASD_varphi = deg2rad(15);  %azimuth angle
ASD_theta = deg2rad(15);   %elevation angle

%% Propagation parameters

%Total uplink transmit power per UE (mW)
p = 100;

%Total downlink transmit power per AP (mW)
rho_tot = 200;

%Prepare to save simulation results


%% Go through all setups
for n = 1:nbrOfSetups
     %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]); 
    [gainOverNoisedB,R,pilotIndex,D] = generateSetupUnicast(L,K,N,tau_p,1,0,ASD_varphi,ASD_theta);               
    %Generate channel realizations with estimates and estimation error correlation matrices
    [Hhat,H,B,C] = functionChannelEstimates(R,nbrOfRealizations,L,sum(K),N,tau_p,pilotIndex,p);
    %Compute SEs for DCC case
    [SE_MR_uni(:,n)] = functionComputeSE_downlink_MR(D,B,tau_c,tau_p,nbrOfRealizations,N,sum(K),L,R,pilotIndex,gainOverNoisedB,rho_tot);
    ASE_MR_uni(tau_p,n) = sum(SE_MR_uni(:,n));
end
  