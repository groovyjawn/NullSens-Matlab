function RESULTS = NullSens(n,p,q,N,M,CD,Type,numIter,numReps,Name)

% This is the main function for the NullSens algorithm. It is formatted for
% experimental evaluation using simulated species responses. An external 
% script may be used to implement a parameter sweep by making multiple 
% calls this this function.

% INPUTS:
% n: number of sites
% p: number of species
% q: number of gradients
% N: noise parameter
% M: covariation parameter
% CD: number of covarying pairs
% Type: type of covariation (0 negative, 1 positive, 2 mixed)
% numIter: number of independent experiments to run
% numReps: number of random matrices generated for null distribution
% Name: File name to save 

% OUTPUT:
% RESULTS: Array containing output data from all experiments
% RESULTS{1}: input parameters in same order as above
% RESULTS{2}: type of covariation assigned to each covarying pair (Used for mixed dynamics only (2): Negative(0), Positive(1))
% RESULTS{3}: coeff. of det. for each species response
% RESULTS{4}: community averaged R2
% RESULTS{5}: adjusted coeff. of det. account for multiple environmental factors
% RESULTS{6}: community average adjusted R2
% RESULTS{7}: test statistic scores for each CDM under test
% RESULTS{8}: P-values for each CDM under test
% RESULTS{9}: array of pairwise covariation for each CDM under test
% RESULTS{10}: matrix of pairwise covariation for each CDM under test
% RESULTS{11}: number pairwise species set to covary using the species
% response, but don't have enough mutual sites to be a valid experiment
% RESULTS{12}: total runtime of all experiments (seconds)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tStart = tic; % Start the clock to measure runtime

% Hard Coded Parameters
Noise = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2]; % Noise Parameters
Magnitude = Noise(N).*[1,3,5,10,20]; % Covariation Magnitude Parameters

% Preallocation for reporting performance
cvid = zeros(numIter,p); % Type of covariation assigned to each covarying pair (Positive or Negative)
R2 = zeros(numIter,p); %Coeff. of Det. for each species response
ComR2 = zeros(numIter,1); % Community Averaged R2
AdjR2 = zeros(numIter,p); % Adjusted Coeff. of Det. account for multiple environmental factors
ComAdjR2 = zeros(numIter,1); % Community Average Adjusted R2
AvgRSq = zeros(1,numIter); % Test Statistic Scores for each CDM under Test
PRSqYres = zeros(1,numIter); % P-values for all experiemnts
CR = cell(1,numIter); % Pairwise Correlation
CV = cell(1,numIter); % Pairwise Covariation
invalid = zeros(1,numIter); % Count number of covarying pairs removed from analysis

for iter = 1:numIter
    % CLEAN START - Clear variables for after experiemnt 
    clearvars -except tStart iter numIter numReps n p q N M CD Type Name Noise Magnitude cvid R2 ComR2 AdjR2 ComAdjR2 AvgRSq AvgRSqYres PRSqYres CR CV invalid RESULTS
    
    % SPECIES RESPONSE MODEL
    % Generate Community Data Matrix (CDM) and Environmetal Factors (X)
    [CDM, X, cvid] = speciesResponseModel(n,p,q,N,M,CD,Type,Noise,Magnitude);
    
    % MULTIVARIATE REGRESSION
    % Sites Selected (regress), Residuals (Yres),
    [~, Yres, ~, regress, R2(iter,:), ComR2(iter), AdjR2(iter,:), ComAdjR2(iter)] = mvrSiteSelectRobust(CDM,X);

    % GENERATE NULL DISTRIBUTION
    AvgRSqYres = zeros(1,numReps); % Preallocate index score (Test Statistic) for each random matrix
    for i = 1:numReps-1
        RANDYres = nullmodel(Yres,regress); % Nullmodel - generates random matrix (RANDYres) for null distribution
        AvgRSqYres(i) = testStatistic(RANDYres,regress,CD); % Index computed on random matrix 
    end
    
    % Index computed on residuals of CDM under test
    [AvgRSqYres(i+1), invalid(iter), CR{iter}, CV{iter}] = testStatistic(Yres,regress,CD); % Index on redsiduals of CDM
    AvgRSq(iter) = AvgRSqYres(numReps);

    % SIGNIFICANCE TEST
    % Obtain P-value - Determine if residuals of CDM exhibit significant index score
    indRYres = find(AvgRSqYres >= AvgRSqYres(numReps));
    PRSqYres(iter) = length(indRYres)/numReps; % Store the P-value for current iteration
end
tStop = toc(tStart); % Total runtime for all experiments (seconds)

RESULTS =  {[n,p,q,N,M,CD,Type,numIter,numReps],cvid,R2,ComR2,AdjR2,ComAdjR2,AvgRSq,PRSqYres,CR,CV,invalid,tStop};
save(Name,'RESULTS')