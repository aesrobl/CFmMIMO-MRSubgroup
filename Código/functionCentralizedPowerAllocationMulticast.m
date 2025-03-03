function rho = functionCentralizedPowerAllocationMulticast(G,traceR,D,rho_tot,wg,upsilon,kappa)
%Compute the power coefficients for the scalable centralized power
%allocation scheme in (7.43)
%
%INPUT:
%K                  = Number of UEs              
%gainOverNoisedB    = Matrix with dimension L x K where (l,k) is the channel
%                     gain (normalized by the noise variance) between AP l
%                     and UE k
%D                  = DCC matrix for cell-free setup with dimension L x K 
%                     where (l,k) is one if AP l serves UE k and zero otherwise
%rho_tot            = Maximum allowed transmit power for each AP 
%wk                 = Matrix with dimension L x K where (l,k) is the
%                     expected value of the norm square of the portion
%                     of the normalized centralized transmit precoder of UE k
%                     corresponding to AP l in (7.37)
%upsilon            = The parameter \upsilon in (7.43)
%kappa              = The parameter \kappa in (7.43)
%
%OUTPUT:
%rho                = Vector of length K with centralized downlink power
%                     allocation coefficients, \rho_k in the monograph
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

%Compute the channel gains in linear domain, i.e., \beta_{kl} in (7.43)
poww = zeros(G,1);
maxPow = zeros(G,1);

for g = 1:G    
    %Extract which UEs are served by AP l
    servingAPs = find(D(:,g)==1);
    %Compute the numerator of (7.43)
    poww(g) =  sum(traceR(servingAPs,g)).^upsilon;
    poww(g) = poww(g)/(max(wg(servingAPs,g))^kappa);
    %Compute \omega_k in (7.41)
    maxPow(g) = max(wg(servingAPs,g));
end

%Prepare to store the normalization factor in (7.43)
normalizationFactor = zeros(G,1);

%Go through all subgroups
for g = 1:G
    %Extract which APs serve subgroup g
    servingAPs = find(D(:,g)==1);
    
    %Gp through all APs that serve subgroup g
    for ell = servingAPs.'
        %Extract which UEs are served by AP ell
        servedSubgroups = find(D(ell,:)==1);
        %Compute the normalization factor in (7.43)
        temporScalar = maxPow(servedSubgroups)'*poww(servedSubgroups)/rho_tot;
        normalizationFactor(g) = max(normalizationFactor(g),temporScalar);
    end
end

%Normalize the numerator terms in (7.43) to obtain \rho_g
rho = poww./normalizationFactor;


