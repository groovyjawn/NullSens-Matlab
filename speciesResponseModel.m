function [CDM,X,cvid,Y,YN,Mag,CDMB,B] = speciesResponseModel(n,p,q,N,M,CD,Type,Noise,Magnitude)

% This function generates a community data matrix (CDM) with n sites and 
% p species to q environmental factors (X), using a linear model. A 
% number, CD of species pairs will covary within the community. The type of 
% covariation depends on Type (0 negative, 1 positive, 2 mixed). The 
% magnitude of covariation depends on Magnitude(M), while the noise depends
% on Noise(N). Both the noise and the covariation magnitude are added to 
% the environmental species responses.

% INPUTS:
% n: number of sites
% p: number of species
% q: number of environmental factors (gradients)
% N: index of noise parameter to select
% M: index of covariation magnitude parameter to select
% CD: number of covarying pairs
% Type: type of community covariation
% Noise: vector of noise parameters
% Magnitude: vector of covariation magnitude parameters

% OUTPUTS:
% CDM:  Community Data Matrix (Sites x Species)
% X:    Environmental Data (Explanatory Variables) (Includes 1's intercept)
% Y:    Environmental Species Responses
% YN:   Noise (Additive)
% Mag:  Magnitude of Species Covariation
% CDMB: Community Data Matrix prior to Zero Censoring
% B:    Species Response Parameters
% cvid: Type of species covariation for each pair (mixed dynamics only)

% Preallocation
B = zeros(q+1,p);
Y = zeros(n,p);
YN = zeros(n,q);
CDM = zeros(n,p);
cvid = zeros(1,p);

% Generate Environmental Factors
X = -1 + (1-(-1)).*rand(n,q);
X = [ones(n,1),X];

% Generate Species Response to Environmental Factors: Covarying Pairs First
for i = 1:2:CD*2
    % Ensure that each covarying pair has at least 7 or 5% of total sites with non-zero abundance
    while length(find(CDM(:,i)>0 & CDM(:,i+1)>0)) <= 0.05*n || length(find(CDM(:,i)>0 & CDM(:,i+1)>0)) <= 7
        B(:,i:i+1) = -1 + (1-(-1)).*rand(q+1,2);
        Y(:,i:i+1) = X*B(:,i:i+1);
        YN(:,i:i+1) = randn(n,2)*Noise(N);
        CDM(:,i:i+1) = Y(:,i:i+1)+YN(:,i:i+1);
        Mag = Magnitude(M)*(-1 + (1-(-1)).*rand(1,n)'); % Covariation vector added/subtracted to species responses
        CDM(:,i) = CDM(:,i)+Mag; % Covariation vector added to first species in pair
        
        if Type == 1 % Positive Dynamics
            CDM(:,i+1) = CDM(:,i+1)+Mag; % If positive covariation add vector to second species in pair
        elseif Type == 0 % Negative Dynamics
            CDM(:,i+1) = CDM(:,i+1)-Mag; % If negative covariation subtract vector to second species in pair
        elseif Type == 2 % Mixed Dynamics
            flip = randi([0,1],1);
            cvid(i:i+1) = flip; % If mixed dynamics then randomly choose to add/subtract vector to second species pair
            if flip == 1 % Positive Dynamics
                CDM(:,i+1) = CDM(:,i+1)+Mag;
            else % Negative Dynamics
                CDM(:,i+1) = CDM(:,i+1)-Mag;
            end
        end
    end
end

% Generate Species Responses for Non-Covarying Pairs
for i = CD*2+1:p
    % Ensure that each species has at least 7 or 5% of total sites with non-zero abundance
    while length(find(CDM(:,i)>0)) <= 0.05*n || length(find(CDM(:,i)>0)) <= 7
        B(:,i) = -1 + (1-(-1)).*rand(q+1,1);
        Y(:,i) = X*B(:,i);
        YN(:,i) = randn(n,1)*Noise(N);
        CDM(:,i) = Y(:,i)+YN(:,i);
    end
end

CDMB = CDM; % CDM before Zero Censoring
rmneg = CDM <= 0;
CDM(rmneg) = 0; % Community Data Matrix Censored at Zero