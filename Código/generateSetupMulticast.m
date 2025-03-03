function [gainOverNoisedB,R,pilotIndexMulti,spatial_sub,DMulti,APpositions,UEpositions,distances] = generateSetupMulticast(L,N,K,G,tau_p,nbrOfSetups,seed,ASD_varphi,ASD_theta)
%This function generates realizations of the simulation setup described in
%
%Section 5.3.
%
%INPUT:
%L               = Number of APs per setup
%K               = Number of UEs in the network
%N               = Number of antennas per AP
%tau_p           = Number of orthogonal pilots
%nbrOfSetups     = Number of setups with random UE and AP locations
%seed            = Seed number of pseudorandom number generator
%ASD_varphi      = Angular standard deviation in the local scattering model 
%                  for the azimuth angle (in radians)
%ASD_theta       = Angular standard deviation in the local scattering model
%                  for the elevation angle (in radians)
%
%OUTPUT:
%gainOverNoisedB = Matrix with dimension L x K x nbrOfSetups where
%                  element (l,k,n) is the channel gain (normalized by the
%                  noise variance) between AP l and UE k in setup n
%R               = Matrix with dimension N x N x L x K x nbrOfSetups
%                  where (:,:,l,k,n) is the spatial correlation matrix
%                  between AP l and UE k in setup n, normalized by noise
%pilotIndex      = Matrix with dimension K x nbrOfSetups containing the
%                  pilots assigned to the UEs
%D               = DCC matrix with dimension L x K x nbrOfSetups where (l,k,n)
%                  is one if AP l serves UE k in setup n and zero otherwise
%                  for cell-free setup
%D_small         = DCC matrix with dimension L x K x nbrOfSetups where (l,k,n)
%                  is one if AP l serves UE k in setup n and zero otherwise
%                  for small-cell setup
%APpositions     = Vector of length L with the AP locations, where the real
%                  part is the horizontal position and the imaginary part
%                  is the vertical position
%UEpositions     = Vector of length K with UE positions, measured in the
%                  same way as APpositions
%distances       = Matrix with same dimension as gainOverNoisedB containing
%                  the distances in meter between APs and UEs
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


%% Define simulation setup

%Set the seed number if it is specified other than zero
if (nargin>6)&&(seed>0)
    rng(seed)
end

%Size of the coverage area (as a square with wrap-around)
squareLength = 1000; %meter

%Size of the clusters of users
C = length(K); %number of clusters of users
squareLengthCluster = 10*ones(1,C); %meter

%Communication bandwidth (Hz)
B = 20e6;

%Noise figure (in dB)
noiseFigure = 7;

%Compute noise power (in dBm)
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Pathloss parameters for the model in (5.42)
alpha = 36.7;
constantTerm = -30.5;

%Standard deviation of the shadow fading in (5.43)
sigma_sf = 4;

%Decorrelation distance of the shadow fading in (5.43)
decorr = 9;

%Height difference between an AP and a UE (in meters)
distanceVertical = 10;

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Prepare to save results
gainOverNoisedB = zeros(L,sum(K),nbrOfSetups);
R = zeros(N,N,L,sum(K),nbrOfSetups);
R_all = zeros(N*L,N*L,sum(K),nbrOfSetups);
distances = zeros(L,sum(K),nbrOfSetups);
pilotIndexUni = zeros(sum(K),nbrOfSetups);
pilotIndexMulti = zeros(G,nbrOfSetups);
masterAPs = zeros(sum(K),1); %the indices of master AP of each UE k 
masterAPsmulti = zeros(G,4); %the indices of master AP of each UE k 

%% Go through all setups
for n = 1:nbrOfSetups   
    %Random AP locations with uniform distribution
    APpositions = (rand(L,1) + 1i*rand(L,1)) * squareLength;   
    %Prepare to compute UE locations
    %UEpositions = zeros(sum(K),1);   
    %Compute alternative AP locations by using wrap around
    wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
    wrapVertical = wrapHorizontal';
    wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
    APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);   
    %Prepare to store shadowing correlation matrix
    shadowCorrMatrix = sigma_sf^2*ones(sum(K),sum(K));
    shadowAPrealizations = zeros(sum(K),L); 
    %Add UEs in clusters
    %Generate center of users' cluster location
    Clusterpositions = (rand(C,1) + 1i*rand(C,1)) * squareLength;
    for c = 1:C
        for k = 1:K(c)
            %Generate a random UE location in the area
            UEpositions(c,k) = Clusterpositions(c) + (rand(1,1) + 1i*rand(1,1)) * squareLengthCluster(c);
        end
    end
    UEpositions = UEpositions.';
    UEpositions = UEpositions(:);
%     plot(UEposition,'*b')
%     hold on
%     plot(APpositionsWrapped,'xr')

    for k = 1:sum(K)
        %Compute distances assuming that the APs are 10 m above the UEs
        [distanceAPstoUE,whichpos] = min(abs(APpositionsWrapped - repmat(UEpositions(k),size(APpositionsWrapped))),[],2);
     %   distanceAPstoUE
        distances(:,k,n) = sqrt(distanceVertical^2+distanceAPstoUE.^2);           
        %If this is not the first UE
        if k-1>0                
            %Compute distances from the new prospective UE to all other UEs
            shortestDistances = zeros(k-1,1);                
            for i = 1:k-1
                shortestDistances(i) = min(abs(UEpositions(k) - UEpositions(i) + wrapLocations));
            end                
            %Compute conditional mean and standard deviation necessary to
            %obtain the new shadow fading realizations, when the previous
            %UEs' shadow fading realization have already been generated.
            %This computation is based on Theorem 10.2 in "Fundamentals of
            %Statistical Signal Processing: Estimation Theory" by S. Kay
            newcolumn = sigma_sf^2*2.^(-shortestDistances/decorr);
            term1 = newcolumn'/shadowCorrMatrix(1:k-1,1:k-1);
            meanvalues = term1*shadowAPrealizations(1:k-1,:);
            stdvalue = sqrt(sigma_sf^2 - term1*newcolumn);                
        else %If this is the first UE                
            %Add the UE and begin to store shadow fading correlation values
            meanvalues = 0;
            stdvalue = sigma_sf;
            newcolumn = [];                
        end            
        %Generate the shadow fading realizations
        shadowing = meanvalues + stdvalue*randn(1,L);
        %shadowing = zeros(1,L);   
        %Compute the channel gain divided by noise power
        gainOverNoisedB(:,k,n) = constantTerm - alpha*log10(abs(distances(:,k,n))) + abs(shadowing)' - noiseVariancedBm;                                  
        %Update shadowing correlation matrix and store realizations
        shadowCorrMatrix(1:k-1,k) = newcolumn;
        shadowCorrMatrix(k,1:k-1) = newcolumn';
        shadowAPrealizations(k,:) = shadowing;  
        %Go through all APs
        for l = 1:L            
            %Compute nominal angle between UE k and AP l
            angletoUE_varphi = angle(UEpositions(k)-APpositionsWrapped(l,whichpos(l))); %azimuth angle
            angletoUE_theta = asin(distanceVertical/distances(l,k,n));  %elevation angle
            %Generate spatial correlation matrix using the local
            %scattering model in (2.18) and Gaussian angular distribution
            %by scaling the normalized matrices with the channel gain
            if nargin>7
                R(:,:,l,k,n) = db2pow(gainOverNoisedB(l,k,n))*functionRlocalscattering(N,angletoUE_varphi,angletoUE_theta,ASD_varphi,ASD_theta,antennaSpacing);
            else
                R(:,:,l,k,n) = db2pow(gainOverNoisedB(l,k,n))*eye(N);  %If angular standard deviations are not specified, set i.i.d. fading
            end
          %  R_all(N*(l-1)+1:N*(l-1)+N,N*(l-1)+1:N*(l-1)+N,k,n) = R(:,:,l,k,n);   
            Beta(k,l,n)=real(trace(R(:,:,l,k,n)))/N;
        end
    end
    
    %Subgroups are created automatically by clusters
%    spa_aux = 1;
%    spa_aux2 = 0;
%    for g = 1:G
%        spa_aux2 = spa_aux2 + K(g);
%        spatial_sub(1,spa_aux:spa_aux2) = g*ones(1,K(g));
%        spa_aux = spa_aux + K(g);
%    end
%    for k=1:sum(K)
%        r_aux(k,:,:)=R_all(:,:,k);
%        r(k,:)=r_aux(k,:);
%         for kk=1:sum(K)
%             Norm(k,kk)=real((trace(R_all(:,:,k)*trace(R_all(:,:,kk))))/trace(R_all(:,:,k)*R_all(:,:,kk)));
%             Nonorm(k,kk)=real(1/trace(R_all(:,:,k)*R_all(:,:,kk)));
%             Frob(k,kk)=norm(R_all(:,:,k)-R_all(:,:,kk),"fro");
%             Cos(k,kk)=1-sum(B(k,1:100).*B(kk,1:100))/(norm(B(k,1:100))*norm(B(kk,1:100)));
%         end
%    end 
%    r(:,end+1) = L*N;
%    spatial_sub(2,:) = kmedoids(r,G,'Distance',@distfun1);
%    Beta
%    G
    spatial_sub = kmedoids(Beta,G,'Distance','cosine');
    spatial_sub = spatial_sub';
    %Determine the master AP for subgroup 1 by looking for AP with best
    %channel condition
    for g = 1:G
        [~,master]=max(mean(gainOverNoisedB(:,spatial_sub==g,n),2));
        DMulti(master,g,n) = 1;
        masterAPsmulti(g) = master; 
        %Assign orthogonal pilots to the first tau_p subgroups according to
        %Algorithm 4.1
        if g <= tau_p
            pilotIndexMulti(g,n) = g;
        else
        %Assign pilot for remaining subgroups           
        %Compute received power to the master AP from each pilot
        %according to Algorithm 4.1
            pilotinterference = zeros(tau_p,1);            
            for t = 1:tau_p            
                aux = find(pilotIndexMulti(1:g-1,n)==t);
                for tt = 1:length(aux)
                    pilotinterference(t) = pilotinterference(t) + sum(db2pow(gainOverNoisedB(master,spatial_sub==aux(tt),n)));                                
                end
                clear aux;
            end
            %Find the pilot with the least receiver power according to
            %Algorithm 4.1
            [~,bestpilot] = min(pilotinterference);
            pilotIndexMulti(g,n) = bestpilot;            
        end
    end
    %Each AP serves the subgroup with the strongest channel condition on each of
    %the pilots in the cell-free setup            
    for l = 1:L        
        for t = 1:max(pilotIndexMulti(:,n))
            pilotsubgroups = find(t==pilotIndexMulti(:,n));
            for tt = 1:length(pilotsubgroups)
                aux(tt) = max(mean(gainOverNoisedB(l,spatial_sub==pilotsubgroups(tt),n),2));
            end
            [~,subgroupindex] = max(aux);
            clear aux;
            DMulti(l,pilotsubgroups(subgroupindex,n)) = 1;         
        end 
    end
end  

    
    


