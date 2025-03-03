function [SE_P_MMSE] = functionComputeSE_downlink_P_MMSE(Hhat,H,D,R,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,gainOverNoisedB,rho_tot)
%Compute downlink SE for different transmit precoding schemes using the capacity
%bound in Theorem 6.1 for the centralized schemes and the capacity bound
%in Corollary 6.3 for the distributed schemes. Compute the genie-aided
%downlink SE from Corollary 6.6 for the centralized and the distributed operations. 
%
%INPUT:
%Hhat              = Matrix with dimension L*N  x nbrOfRealizations x K
%                    where (:,n,k) is the estimated collective channel to
%                    UE k in channel realization n.
%H                 = Matrix with dimension L*N  x nbrOfRealizations x K
%                    with the true channel realizations. The matrix is
%                    organized in the same way as Hhat.
%D                 = DCC matrix for cell-free setup with dimension L x K 
%                    where (l,k) is one if AP l serves UE k and zero otherwise
%B                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix of the channel estimate 
%                    between AP l and UE k, normalized by noise variance
%C                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix of the channel
%                    estimation error between AP l and UE k,
%                    normalized by noise variance
%tau_c             = Length of coherence block
%tau_p             = Length of pilot sequences
%nbrOfRealizations = Number of channel realizations
%N                 = Number of antennas per AP
%K                 = Number of UEs 
%L                 = Number of APs
%p                 = Uplink transmit power per UE (same for everyone)
%R                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix between AP l and UE k,
%                    normalized by noise
%pilotIndex        = Vector containing the pilot assigned to each UE
%rho_dist          = Matrix with dimension L x K where (l,k) is the power
%                    allocated to UE k by AP l in the distributed downlink
%                    operation
%gainOverNoisedB   = Matrix with dimension L x K where (l,k) is the channel
%                    gain (normalized by the noise variance) between AP l
%                    and UE k
%rho_tot           = Maximum allowed transmit power for each AP 
%
%OUTPUT:
%SE_MMSE           = SEs achieved with MMSE precoding in (6.16)
%SE_P_MMSE         = SEs achieved with P-MMSE precoding in (6.17)
%SE_P_RZF          = SEs achieved with P-RZF precoding in (6.18)
%SE_L_MMSE         = SEs achieved with L-MMSE precoding in (6.25)
%SE_LP_MMSE        = SEs achieved with LP-MMSE precoding in (6.33)
%SE_MR             = SEs achieved with MR precoding in (6.26)
%Gen_SE_P_MMSE     = Genie-aided SEs achieved with P-MMSE precoding in (6.17)
%Gen_SE_P_RZF      = Genie-aided SEs achieved with P-RZF precoding in (6.18)
%Gen_SE_LP_MMSE    = Genie-aided SEs achieved with LP-MMSE precoding in (6.33)
%Gen_SE_MR         = Genie-aided SEs achieved with MR precoding in (6.26)
%
%
%This Matlab function was developed to generate simulation results to:
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

%Store the N x N identity matrix
eyeN = eye(N);

%Compute the prelog factor assuming only downlink data transmission
prelogFactor = (1-tau_p/tau_c);

%Prepare to store the terms that appear in SEs
signal_P_MMSE = zeros(K,1);
interf_P_MMSE = zeros(K,1);
scaling_P_MMSE = zeros(K,1);
portionScaling_PMMSE = zeros(L,K);
interUserGains_P_MMSE = zeros(K,K,nbrOfRealizations);

%% Compute scaling factors for precoding
%Go through all channel realizations
for n = 1:nbrOfRealizations  
    %Go through all UEs
    for k = 1:K        
        %Determine the set of serving APs
        servingAPs = find(D(:,k)==1);
        La = length(servingAPs);        
        %Determine which UEs that are served by partially the same set
        %of APs as UE k, i.e., the set in (5.15)
        servedUEs = sum(D(servingAPs,:),1)>=1;                
        %Extract channel realizations and estimation error correlation
        %matrices for the APs that involved in the service of UE k
        Hhatallj_active = zeros(N*La,K);
        C_tot_blk_partial = zeros(N*La,N*La);        
        R_tot_blk_partial = zeros(N*La,N*La);        
        for l = 1:La
            Hhatallj_active((l-1)*N+1:l*N,:) = reshape(Hhat((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            C_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),servedUEs),4);
            R_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(R(:,:,servingAPs(l),~servedUEs),4);            
        end        
        %Compute MMSE, P-MMSE, and P-RZF precoding
        V_P_MMSE = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+p*(C_tot_blk_partial+R_tot_blk_partial)+eye(La*N))\Hhatallj_active(:,k));        
        %Compute scaling factor by Monte Carlo methods
        scaling_P_MMSE(k) = scaling_P_MMSE(k) + sum(abs(V_P_MMSE).^2,1)/nbrOfRealizations;        
        %Go through all the serving APs
        for l=1:La
            %Extract the portions of the centralized precoding vectors 
            V_P_MMSE2 = V_P_MMSE((l-1)*N+1:l*N,:);
            portionScaling_PMMSE(servingAPs(l),k) = portionScaling_PMMSE(servingAPs(l),k) ...
                + sum(abs(V_P_MMSE2).^2,1)/nbrOfRealizations;
        end
    end    
end
%Normalize the norm squares of the portions for the normalized centralized precoders
portionScaling_PMMSE = portionScaling_PMMSE./repmat(scaling_P_MMSE.',[L 1]);
%The parameters for the scalable centralized downlink power allocation in (7.43)
upsilon = -0.5;
kappa = 0.5;
%Compute the power allocation coefficients for centralized precoding according to (7.43)
rho_PMMSE = functionCentralizedPowerAllocation(K,gainOverNoisedB,D,rho_tot,portionScaling_PMMSE,upsilon,kappa);
%rho_PMMSE = rho_tot/K*ones(K,1);
%scaling_P_MMSE
%% Go through all channel realizations
%nbrOfRealizations = 100;
for n = 1:nbrOfRealizations       
    %Matrix to store Monte-Carlo results for this realization
    interf_P_MMSE_n = zeros(K,K);  
    %Consider the centralized schemes        
    %Go through all UEs
    for k = 1:K                
        %Determine the set of serving APs
        servingAPs = find(D(:,k)==1);       
        La = length(servingAPs);       
        %Determine which UEs that are served by partially the same set
        %of APs as UE k, i.e., the set in (5.15)
        servedUEs = sum(D(servingAPs,:),1)>=1;        
        %Extract channel realizations and estimation error correlation
        %matrices for the APs that involved in the service of UE k
        Hallj_active = zeros(N*La,K);        
        Hhatallj_active = zeros(N*La,K);
        C_tot_blk_partial = zeros(N*La,N*La);        
        R_tot_blk_partial = zeros(N*La,N*La);              
        for l = 1:La
            Hallj_active((l-1)*N+1:l*N,:) = reshape(H((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            Hhatallj_active((l-1)*N+1:l*N,:) = reshape(Hhat((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            C_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),servedUEs),4);
            R_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(R(:,:,servingAPs(l),~servedUEs),4);            
        end
        %Compute P-MMSE precoding
        w = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+p*(C_tot_blk_partial+R_tot_blk_partial)+eye(La*N))\Hhatallj_active(:,k));
        %Apply power allocation
        w = w*sqrt(rho_PMMSE(k)/scaling_P_MMSE(k));      
       % h1wuni = abs(H(:,n,1)'*w)
       % h2wuni = abs(H(:,n,2)'*w)
        %Compute realizations of the terms inside the expectations
        %of the signal and interference terms in Theorem 6.1
        signal_P_MMSE(k) = signal_P_MMSE(k) + (Hallj_active(:,k)'*w)/nbrOfRealizations;
        interf_P_MMSE_n(:,k) = interf_P_MMSE_n(:,k) + Hallj_active'*w;        
        %Compute gain of the signal from UE that arrives at other UEs
        interUserGains_P_MMSE(:,k,n) = interUserGains_P_MMSE(:,k,n) + Hallj_active'*w; 
    end       
    %Compute interference power in one realization
   interf_P_MMSE = interf_P_MMSE + sum(abs(interf_P_MMSE_n).^2,2)/nbrOfRealizations;   
end

%% Compute the SEs
%Compute SE in Theorem 6.1 with P-MMSE
SE_P_MMSE = prelogFactor*real(log2(1+(abs(signal_P_MMSE).^2) ./ (interf_P_MMSE - abs(signal_P_MMSE).^2 + 1)));
SE_unicast = sum(SE_P_MMSE)
%Remove unused large matrices
clear interUserGains_MR interUserGains_MMSE interUserGains_P_MMSE interUserGains_P_RZF interUserGains_L_MMSE interUserGains_LP_MMSE;
