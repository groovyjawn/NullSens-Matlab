function [AvgRSq, invalid, CR, CV] = testStatistic(matrix,regress,CD)

% This function computes the test statistic for matrix. The test statistic
% is a single number the quantifies the total covariation among all species
% pairs in matrix. Each species pair is considered only at mutual selected
% sites.

% INPUTS:
% matrix: matrix to compute test statistic
% regress: sites selected for each species
% CD: number of covarying pairs

% OUTPUTS:
% AvgRSq: test statistic for matrix
% invalid: number pairwise species set to covary using the species response
% model that do not have a sufficient number of mutual sites selected
% CR: matrix of pairwise correlation 
% CV: matrix of pairwise covariation

[n,p] = size(matrix); % n sites, p species

% Preallocate and Initialize
CR = zeros(p-1,p);
CV = zeros(p-1,p);
w = 0;
AvgRSq = 0;
invalid = 0;

% For each species pair (i,j):
for i = 1:p-1
    for j = i+1:p        
        % Count the number of mutual sites selected per species pair (i,j)
        temp = and(regress{i},regress{j}); 
        if sum(temp) <= n*0.05 || sum(temp) <= 7
            cvar = 0; % Set covariance between species pairs with
            % insufficent number of mutual to zero
            
            % Count number of covarying species pairs with insufficient
            % number of mutual sites - these are species pairs that have
            % been set to covary using the species response model
            if i == 1 && j == 2 && CD > 0
                invalid = invalid + 1;
            end
            if i == 3 && j == 4 && CD > 1
                invalid = invalid + 1;
            end
            if i == 5 && j == 6 && CD > 2
                invalid = invalid + 1;
            end
            if i == 7 && j == 8 && CD > 3
                invalid = invalid + 1;
            end
            if i == 9 && j == 10 && CD > 4
                invalid = invalid + 1;
            end
            if i == 11 && j == 12 && CD > 5
                invalid = invalid + 1;
            end
            if i == 13 && j == 14 && CD > 6
                invalid = invalid + 1;
            end
            if i == 15 && j == 16 && CD > 7
                invalid = invalid + 1;
            end
            if i == 17 && j == 18 && CD > 8
                invalid = invalid + 1;
            end
            if i == 19 && j == 20 && CD > 9
                invalid = invalid + 1;
            end
        else
            COV = cov(matrix(temp,[i,j])); % Covariance of species pair (i,j)
            CORR = corr(matrix(temp,[i,j])); % Correlation of species pair (i,j)
            if isnan(COV(1,2)) == 1 % Check for NaN, set to zero if True
                cvar = 0;
                CR(i,j) = 0;
                CV(i,j) = 0;
            else
                % COV and CORR are 2x2 matrices - the following extracts
                % the required value from the matrix
                cvar = COV(1,2);
                CR(i,j) = CORR(1,2);
                CV(i,j) = cvar;
            end 
        end
        
        % Test Statistic: 
        AvgRSq = AvgRSq + cvar.^2; % Sum of all pairwise covariance squared
        w = w + abs(cvar); % Sum of all pairwise absolute covariance  
    end
end


if w == 0
    AvgRSq = 0;
else
    AvgRSq = AvgRSq/w; % Test Statistic
end