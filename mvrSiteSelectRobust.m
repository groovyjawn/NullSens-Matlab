function [Yhat, Yres, BP, regress, R2, ComR2, AdjR2, ComAdjR2, C, T] = mvrSiteSelectRobust(CDM,X)

% This function first determines which sites to include in subsequent
% analysis, per species. It then uses Robust multivariate regression to
% estimate species response parameters to the environmental factors in X.
% Fitted and residual response are then obtained. The coefficient of
% determination is then computed to provide a measure of model fit.

% INPUTS:
% CDM: Community Data Matrix
% X: Environmental Factors

% OUTPUTS:
% Yhat: Fitted Species Response Matrix
% Yres: Residual Species Response Matrix
% BP: Estimate Species Response Parameters
% regress: Sites selected for analysis for each species
% R2: Coeff. of Det. for each species response
% ComR2: Community Averaged R2
% AdjR2: Adjusted Coeff. of Det. account for multiple environmental factors
% ComAdjR2: Community Average Adjusted R2
% C: Array of centroids used in site selection
% T: Array of thresholds used in site selection

[n,p] = size(CDM); % Get the number of sites, n and species, p
[~,q] = size(X); % Get the number of gradients

% Preallocation
BP = zeros(q,p); % Estimated Species response Parameters
regress = cell(1,p); % Array to hold sites selected for each species
Yhat = zeros(size(CDM)); % Fitted Species Responses
Yres = zeros(size(CDM)); % Residual Species Responses
R2 = zeros(1,p); % Coefficient of Determination
AdjR2 = zeros(1,p); % Adjusted Coeff. of Det. 
C = cell(1,p); % Centroids for species' site selection
T = cell(1,p); % Thresholds for species' site selection

% For each species in CDM:
for i = 1:p
    % SITE SELECTION PROCEDURE
    % Determine sites to include in regression analysis
    posAb = find(CDM(:,i)>0); % All sites with non-zero abundance
    remove = []; % Initialize vector of sites for removal from analysis
    centroid = zeros(1,q); % Preallocate
    thresh = zeros(1,q); % Preallocate
    for grad = 2:q % q = 1 refers to column vector of ones in X
        % Compute median value of environmental factor, grad using sites
        % with non-zero abundance
        centroid(grad) = prctile(X(posAb,grad),50);
        % Compute the distance of each non-zero abundance from centroid
        dist = squareform(pdist([centroid(grad); X(posAb,grad)]));
        % Set threshold at the mean distance plus one standard deviation
        thresh(grad) = mean(dist(1,2:end)) + std(dist(1,2:end));
        
        % Compute the distance of every site from the centroid. This
        % includes sites with both zero and non-zero abundance.
        dist = squareform(pdist([centroid(grad);X(:,grad)]));
        
        % Remove those sites that are outside of the threshold
        remove_temp = find(dist(1,2:end) > thresh(grad));
        remove = [remove,remove_temp];    %#ok<AGROW>
    end
    C{i} = centroid(2:q); % Record the centroid
    T{i} = thresh(2:q); % Record the threshold
   
    % Determine which sites are kept and store them in regress as True
    keep = setxor(1:n,remove);
    regress{i} = false(n,1);
    regress{i}(keep) = true; 
    
    % MULTIVARIATE REGRESSION
    % Try Robust Regression first, then Standard Regression if error
    try
        BP(:,i) = robustfit(X(regress{i},:),CDM(regress{i},i),[],[],'off');
    catch  %#ok<CTCH>
        ['Species ', num2str(i), ': Robustfit Exception - Using Standard MVR'] %#ok<NOPRT>
        BP(:,i) = (X(regress{i},:)'*X(regress{i},:))\(X(regress{i},:)'*CDM(regress{i},i));
    end
    
    % If NaN, set to zero
    rm = isnan(BP(:,i));
    BP(rm,i) = 0;
    
    % FITTED AND RESIDUAL RESPONSES
    Yhat(:,i) = X*BP(:,i); % Compute the fitted responses
    rmneg = Yhat(:,i) <= 0; % Find fitted responses with negative abundance
    Yhat(rmneg,i) = 0;  % Censor fitted responses to zero
    Yres(:,i) = CDM(:,i)-Yhat(:,i); % Compute residual responses
    
    
    % COEFFICIENT OF DETERMINATION (R2)
    ssres = sum((CDM(:,i)-Yhat(:,i)).^2);
    sstot = sum((CDM(:,i) - mean(CDM(:,i))).^2);
    R2(i) = 1-ssres/sstot;
    NN(i) = n; % Number of sites included in R2 calculation
end

% Community Average R2
ComR2 = mean(R2);

% Adjusted R2
AdjR2 = 1-((1-R2).*(NN-1)./(NN-q-1));

% Weighted Community Averaged Adjusted R2
ComAdjR2 = mean(AdjR2);