function [SE_MR] = functionComputeSE_downlink_MR(D,B,tau_c,tau_p,nbrOfRealizations,N,K,L,R,pilotIndex,gainOverNoisedB,rho_tot)
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
signal_MR = zeros(K,1);
interf_MR = zeros(K,1);
cont_MR  = zeros(K,K);
scaling_MR = zeros(L,K);
interUserGains_MR = zeros(K,K,nbrOfRealizations);

%Compute the power allocation in (6.36) for distributed precoding
%Compute the power allocation in (6.36) for distributed precoding
rho_dist = zeros(L,K);
gainOverNoise = db2pow(gainOverNoisedB);

for l = 1:L
    
    %Extract which Subgroups are served by AP l
    servedUEs = find(D(l,:)==1);
    
    %Compute denominator in (6.36)
    normalizationAPl = sum(sqrt(gainOverNoise(l,servedUEs)));
    
    for ind = 1:length(servedUEs)
        
        rho_dist(l,servedUEs(ind)) = rho_tot*sqrt(gainOverNoise(l,servedUEs(ind)))/normalizationAPl;
        

    end
    
end

%% Compute scaling factors for precoding
%Computation for MR precoding
for l = 1:L
    
    %Extract which UEs are served by AP l
    servedUEs = find(D(l,:)==1);
    
    for ind = 1:length(servedUEs)
        
        %Compute scaling factor using the spatial correlation matrix of the
        %channel estimate
        scaling_MR(l,servedUEs(ind)) = trace(B(:,:,l,servedUEs(ind)));
        
    end
    
end

%% Compute MR closed-form expectations

%Go through all APs
for l = 1:L
    
    %Extract which UEs are served by the AP
    servedUEs = find(D(l,:)==1);
    
    %Go through all UEs served by the AP
    for ind = 1:length(servedUEs)
        
        %Extract UE index
        k = servedUEs(ind);
        
        %Desired signal term in (6.27)
        signal_MR(k) = signal_MR(k) + sqrt(rho_dist(l,k)*real(trace(B(:,:,l,k))));
        
        
        for i = 1:K
            
            
            %Non-coherent interference from UE k to UE i (the first term of
            %(6.28))
            interf_MR(i) = interf_MR(i) + rho_dist(l,k)*real(trace(B(:,:,l,k)*R(:,:,l,i)))/real(trace(B(:,:,l,k)));
            
            if pilotIndex(k) == pilotIndex(i)
                
                %Coherent interference from UE k to UE i (the second term
                %of (6.28))
                cont_MR(i,k) = cont_MR(i,k) + sqrt(rho_dist(l,k))*real(trace((B(:,:,l,k)/R(:,:,l,k))*R(:,:,l,i)))/sqrt(real(trace(B(:,:,l,k))));
                
            end
            
        end
        
    end
    
end

%% Compute the SEs
%Compute SE in Corollary 6.3 with MR  using the closed-form expressions in Corollary 6.4
SE_MR = prelogFactor*real(log2(1+(abs(signal_MR).^2) ./ (interf_MR + sum(abs(cont_MR).^2,2) - abs(signal_MR).^2 + 1)));

%Remove unused large matrices
clear interUserGains_MR interUserGains_MMSE interUserGains_P_MMSE interUserGains_P_RZF interUserGains_L_MMSE interUserGains_LP_MMSE;
