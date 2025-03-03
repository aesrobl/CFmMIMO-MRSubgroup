function [SE_MR] = functionComputeSE_downlink_multicast_MR_Normalized(HhatG,H,D,A,tau_c,tau_p,nbrOfRealizations,N,K,G,L,rho_tot,spatial_sub)
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
%prelogFactor = (1-min(tau_p,G)/tau_c);
prelogFactor = (1-tau_p/tau_c);
%Prepare to store the terms that appear in SEs
signal_MR = zeros(K,1);
interf_MR = zeros(K,1);
scaling_MR_normalized = zeros(L,G);
nu = 0.5;

for l = 1:L
    
    %Extract which Subgroups are served by AP l
    servedSubgroups = find(D(l,:)==1);
    
    %Compute denominator in (6.36)
    for ind = 1:length(servedSubgroups)
        aux_normalizationAPl(ind) = trace(A(:,:,l,servedSubgroups(ind)))^nu;
    end
    normalizationAPl = sum(aux_normalizationAPl);
    for ind = 1:length(servedSubgroups)
        rho_dist(l,servedSubgroups(ind)) = rho_tot*(trace(A(:,:,l,servedSubgroups(ind)))^nu)/normalizationAPl;
        %Compute scaling factor using the spatial correlation matrix of the
        %channel estimate       
        scaling_MR_normalized(l,servedSubgroups(ind)) = 1; 
        %Para normalized el escalado se puede simplificar como 
        %una matriz unitaria en comparación con el CB donde usamos 
        %la traza de B     
    end    
end

%% Go through all channel realizations
for n = 1:nbrOfRealizations   
    %Matrix to store Monte-Carlo results for this realization
    interf_MR_n = zeros(K,G);  
    %Consider the centralized schemes 
    %Go through all subgroups to obtain the precoders
    for g = 1:G        
        %Determine the set of serving APs that give connection to the g
        %subgroup
        servingAPs = find(D(:,g)==1);
        La = length(servingAPs);                     
        %Variables to store after extract channel realizations and estimation error correlation
        %matrices for the APs that involved in the service of UE k
        Hallj_active = zeros(N*La,K);        
        HhatGallj_active = zeros(N*La,G); 
        for l = 1:La
            Hallj_active((l-1)*N+1:l*N,:) = reshape(H((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            HhatGallj_active((l-1)*N+1:l*N,:) = reshape(HhatG((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N G]);
            %Compute MR precoding for Normalized Case
            w((l-1)*N+1:l*N) = HhatGallj_active((l-1)*N+1:l*N,g)/norm(HhatGallj_active((l-1)*N+1:l*N,g));
            %Apply power allocation
            w((l-1)*N+1:l*N) = w((l-1)*N+1:l*N)*sqrt(rho_dist(servingAPs(l),g)/(scaling_MR_normalized(servingAPs(l),g)));     
        end
        w = w.'; %Conjugate Transposed Operation on 
        %h1wmulti = abs(H(:,n,1)'*w)
        %h1wmulti = abs(H(:,n,2)'*w)
        %Compute realizations of the terms inside the expectations
        %of the signal and interference terms in Theorem 6.1
        %Go through all UEs in the subgroup   
        for k = find(g==spatial_sub)'                     
            %Extract channel realizations and estimation error correlation
            %matrices for the APs that involved in the service of UE k
            signal_MR(k) = signal_MR(k) + (Hallj_active(:,k)'*w)/nbrOfRealizations;
        end
        interf_MR_n(:,g) = interf_MR_n(:,g) + Hallj_active'*w;  
        clear w;
    end
    %Compute interference power in one realization
    interf_MR = interf_MR + sum(abs(interf_MR_n).^2,2)/nbrOfRealizations;   
end

%% Compute the SEs
%Compute SE in Theorem 6.1 with MR
SE_MR = prelogFactor*real(log2(1+(abs(signal_MR).^2) ./ (interf_MR - (abs(signal_MR).^2) + 1)));
for g = 1:G
    SE_MR(spatial_sub==g) = min(SE_MR(spatial_sub==g));
end
display(sum(SE_MR))