function [SE_P_MMSE] = functionComputeSE_downlink_multicast_IPMMSE(HhatG,Hhat,H,D,A,C,tau_c,tau_p,nbrOfRealizations,N,K,G,L,p,rho_tot,spatial_sub)
%Compute downlink SE using Improved Partial MMSE transmit precoding scheme.
%
%INPUT:
%HhatG
%Hhat              = Matrix with dimension L*N  x nbrOfRealizations x K
%                    where (:,n,k) is the estimated collective channel to
%                    UE k in channel realization n.
%H                 = Matrix with dimension L*N  x nbrOfRealizations x K
%                    with the true channel realizations. The matrix is
%                    organized in the same way as Hhat.
%D                 = DCC matrix for cell-free setup with dimension L x K 
%                    where (l,k) is one if AP l serves UE k and zero otherwise
%A
%C                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix of the channel
%                    estimation error between AP l and UE k,
%                    normalized by noise variance
%tau_c             = Length of coherence block
%tau_p             = Length of pilot sequences
%nbrOfRealizations = Number of channel realizations
%N                 = Number of antennas per AP
%K                 = Number of UEs 
%G
%L                 = Number of APs
%p                 = Uplink transmit power per UE (same for everyone)
%rho_tot           = Maximum allowed transmit power for each AP 
%spatial_sub
%OUTPUT:

%SE_P_MMSE         = SEs achieved with IP-MMSE precoding
%
%
%This Matlab function was developed to generate simulation results to:
%
%Alejandro de la Fuente (2023)
%
%Subgrouping in User-Centric Scalable Cell-Free Massive MIMO Multicasting,
%
%
%This is version 1.0 (Last edited: 2023-10-25)
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
scaling_P_MMSE = zeros(G,1);
portionScaling_PMMSE = zeros(L,G);
Kg = zeros(1,G);
for g = 1:G   
    %Determine number of users in each subgroup        
    Kg(g) = sum(spatial_sub==g);
end
%% Compute scaling factors for precoding
%Go through all channel realizations
for n = 1:nbrOfRealizations  
    %Go through all subgroups
    for g = 1:G   
        %Determine the set of serving APs
        servingAPs = find(D(:,g)==1);
        La = length(servingAPs);        
        %Determine which subgroups that are served by partially the same set
        %of APs as sugroup g, i.e., the set in (5.15)
        servedSubgroups = sum(D(servingAPs,:),1)>=1;  
        servedSubgroups_aux = find(servedSubgroups == 1);
        servedSubgroups_aux2 = find(servedSubgroups == 0);
        %Extract channel realizations and estimation error correlation
        %matrices for the APs that involved in the service of UE k
        HhatGallj_active_aux = zeros(N*La,G);
        HhatGallj_active = zeros(N*La,length(servedSubgroups_aux));
        C_tot_blk_partial = zeros(N*La,N*La);   
        R_tot_blk_partial = zeros(N*La,N*La);   
        C_tot_blk_partial_aux = zeros(N,N,length(servedSubgroups_aux)); 
        R_tot_blk_partial_aux = zeros(N,N,length(servedSubgroups_aux2));               
        for l = 1:La
            HhatGallj_active_aux((l-1)*N+1:l*N,:) = reshape(HhatG((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N G]);
            for ind = 1:G
                HhatGallj_active((l-1)*N+1:l*N,ind) = Kg(ind).*HhatGallj_active_aux((l-1)*N+1:l*N,ind);
            end                
            for ind = 1:length(servedSubgroups_aux)
                C_tot_blk_partial_aux(:,:,ind) = (Kg(servedSubgroups_aux(ind)).^2).*C(:,:,servingAPs(l),servedSubgroups_aux(ind));
            end 
            for ind = 1:length(servedSubgroups_aux2)
                R_tot_blk_partial_aux(:,:,ind) = (Kg(servedSubgroups_aux2(ind)).^2).*A(:,:,servingAPs(l),servedSubgroups_aux2(ind));
            end 
            C_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C_tot_blk_partial_aux,3);
            R_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(R_tot_blk_partial_aux,3);            
        end        
        %Compute P-MMSE precoding      
        V_P_MMSE = Kg(g)*sqrt(p/(tau_p*p))*(((p/(tau_p*p))*(HhatGallj_active(:,servedSubgroups)*HhatGallj_active(:,servedSubgroups)')+(p/(tau_p*p))*(C_tot_blk_partial+R_tot_blk_partial)+eye(La*N))\HhatGallj_active(:,g));       
        %Compute scaling factor by Monte Carlo methods
        scaling_P_MMSE(g) = scaling_P_MMSE(g) + sum(abs(V_P_MMSE).^2,1)/nbrOfRealizations;   
        %Go through all the serving APs
        for l=1:La
            %Extract the portions of the centralized precoding vectors 
            V_P_MMSE2 = V_P_MMSE((l-1)*N+1:l*N,:);
            portionScaling_PMMSE(servingAPs(l),g) = portionScaling_PMMSE(servingAPs(l),g) ...
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
for l = 1:L
    servedSubgroups = find(D(l,:)==1);
    for ind = 1:length(servedSubgroups)
        traceR(l,servedSubgroups(ind)) = trace(A(:,:,l,servedSubgroups(ind)));
    end
end
rho_PMMSE = functionCentralizedPowerAllocationMulticast(G,traceR,D,rho_tot,portionScaling_PMMSE,upsilon,kappa);

%% Go through all channel realizations
%nbrOfRealizations = 1;
for n = 1:nbrOfRealizations   
    %Matrix to store Monte-Carlo results for this realization
    interf_P_MMSE_n = zeros(K,G);  
    %Consider the centralized schemes 
    %Go through all subgroups to obtain the precoders
    for g = 1:G     
        %Determine the set of serving APs
        servingAPs = find(D(:,g)==1);
        La = length(servingAPs);        
        %Determine which subgroups that are served by partially the same set
        %of APs as subgroup g
        servedSubgroups = sum(D(servingAPs,:),1)>=1;       
        servedSubgroups_aux = find(servedSubgroups == 1);
        servedSubgroups_aux2 = find(servedSubgroups == 0);
        %Extract channel realizations and estimation error correlation
        %matrices for the APs that involved in the service of UE k
        Hallj_active = zeros(N*La,K);        
        Hhatallj_active = zeros(N*La,K);
        HhatGallj_active_aux = zeros(N*La,G);
        HhatGallj_active = zeros(N*La,length(servedSubgroups_aux));
        C_tot_blk_partial = zeros(N*La,N*La);   
        R_tot_blk_partial = zeros(N*La,N*La);   
        C_tot_blk_partial_aux = zeros(N,N,length(servedSubgroups_aux)); 
        R_tot_blk_partial_aux = zeros(N,N,length(servedSubgroups_aux2));               
        for l = 1:La
            Hallj_active((l-1)*N+1:l*N,:) = reshape(H((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            Hhatallj_active((l-1)*N+1:l*N,:) = reshape(Hhat((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);             
            HhatGallj_active_aux((l-1)*N+1:l*N,:) = reshape(HhatG((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N G]);
            for ind = 1:G
                HhatGallj_active((l-1)*N+1:l*N,ind) = Kg(ind).*HhatGallj_active_aux((l-1)*N+1:l*N,ind);
            end                          
            for ind = 1:length(servedSubgroups_aux)
                C_tot_blk_partial_aux(:,:,ind) = (Kg(servedSubgroups_aux(ind)).^2).*C(:,:,servingAPs(l),servedSubgroups_aux(ind));
            end 
            for ind = 1:length(servedSubgroups_aux2)
                R_tot_blk_partial_aux(:,:,ind) = (Kg(servedSubgroups_aux2(ind)).^2).*A(:,:,servingAPs(l),servedSubgroups_aux2(ind));
            end 
            C_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C_tot_blk_partial_aux,3);
            R_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(R_tot_blk_partial_aux,3);            
        end                   
        %Compute P-MMSE precoding
        w = Kg(g)*sqrt(p/(tau_p*p))*(((p/(tau_p*p))*(HhatGallj_active(:,servedSubgroups)*HhatGallj_active(:,servedSubgroups)')+(p/(tau_p*p))*(C_tot_blk_partial+R_tot_blk_partial)+eye(La*N))\HhatGallj_active(:,g));         
        %Apply power allocation
        w = w*sqrt(rho_PMMSE(g)/(scaling_P_MMSE(g)));     
   
        %Compute realizations of the terms inside the expectations
        %of the signal and interference terms
        %Go through all UEs in the subgroup   
        for k = find(g==spatial_sub)'                     
            %Extract channel realizations and estimation error correlation
            %matrices for the APs that involved in the service of UE k
            signal_P_MMSE(k) = signal_P_MMSE(k) + (Hallj_active(:,k)'*w)/nbrOfRealizations;
        end
        interf_P_MMSE_n(:,g) = interf_P_MMSE_n(:,g) + Hallj_active'*w;  
    end
    %Compute interference power in one realization
    interf_P_MMSE = interf_P_MMSE + sum(abs(interf_P_MMSE_n).^2,2)/nbrOfRealizations;   
end

%% Compute the SEs
%Compute SE with P-MMSE
SE_P_MMSE = prelogFactor*real(log2(1+(abs(signal_P_MMSE).^2) ./ (interf_P_MMSE - (abs(signal_P_MMSE).^2) + 1)));
for g = 1:G
    SE_P_MMSE(spatial_sub==g) = min(SE_P_MMSE(spatial_sub==g));
end
display(sum(SE_P_MMSE))

