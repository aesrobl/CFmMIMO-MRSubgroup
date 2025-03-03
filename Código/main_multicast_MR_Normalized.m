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

G_Values = [1, 10, 25, 50, 75, 100];
for g = 1:length(G_Values)
   
    %Empty workspace and close figures
    clearvars -except g G_Values
    close all;
    
    tic;

    %% Define simulation setup
    
    %Number of Monte-Carlo setups
    nbrOfSetups = 500;
    
    %Number of channel realizations per setup
    nbrOfRealizations = 100;
    
    %Number of APs 
    L = 100;
    
    %Number of antennas per AP
    N = 16;
    
    %Number of UEs in the network per cluster
    Kc = 10;
    Cluster = 10;
    K = Kc*ones(1,Cluster);
    
    %Number of subgroups
    G = G_Values(g);
    
    %Length of coherence block
    tau_c = 200;
    
    %Length of pilot sequences (depending on the number of subgroups)
    if G < 10
        tau_p = G;
    else
        tau_p = 10;
    end
    
    % %Angular standard deviation in the local scattering model (in radians)
    % ASD_varphi = 0;  %azimuth angle
    % ASD_theta = 0;   %elevation angle
    
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
        [gainOverNoisedB,R,pilotIndexMulti,spatial_sub,DMulti] = generateSetupMulticast(L,N,K,G,tau_p,1,n);          
        %Generate channel realizations with estimates and estimation error correlation matrices and compute SEs for subgrouping case    
        [HhatG,HhatMulti,HMulti,AMulti,BMulti,CMulti] = functionChannelEstimatesMulticast(R,nbrOfRealizations,L,sum(K),N,tau_p,pilotIndexMulti,p,spatial_sub);      
        [SE_MR_normalized_multi(:,n)] = functionComputeSE_downlink_multicast_MR_Normalized(HhatG,HMulti,DMulti,AMulti,tau_c,tau_p,nbrOfRealizations,N,sum(K),G,L,rho_tot,spatial_sub);
        [ASE_MR_normalized_multi(n)] = sum(SE_MR_normalized_multi(:,n));
    end
    
    results_filename = [num2str(L) 'x' num2str(N) '-' num2str(Cluster) 'x' num2str(Kc) '-MR-normalized-multi' num2str(G) '.mat'];
    save(results_filename,'SE_MR_normalized_multi','ASE_MR_normalized_multi');
    time = toc;
    disp(['El código se ha ejecutado en ' datestr(seconds(time), 'HH:MM:SS') '.']);
end
  