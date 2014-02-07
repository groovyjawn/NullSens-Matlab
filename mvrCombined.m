function [Yhat, Yres, BP, regress, AR, ARC, AVGAR, AVGARC, C, T] = mvrCombined(CDM,X)

[n,p] = size(CDM);
[~,q] = size(X);

BP = zeros(q,p);
regress = cell(1,p);
Yhat = zeros(size(CDM));
Yres = zeros(size(CDM));
AR = zeros(1,p);
ARC = zeros(1,p);
C = cell(1,p);
T = cell(1,p);

for i = 1:p
    % Determine Points to Include in Regression Analysis
    posAb = find(CDM(:,i)>0);
    remove = [];
    centroid = zeros(1,q);
    thresh = zeros(1,q);
    for grad = 2:q
        centroid(grad) = prctile(X(posAb,grad),50);
        dist = squareform(pdist([centroid(grad); X(posAb,grad)]));
        thresh(grad) = mean(dist(1,2:end)) + std(dist(1,2:end));
        
        dist = squareform(pdist([centroid(grad);X(:,grad)]));
        remove_temp = find(dist(1,2:end) > thresh(grad));
        remove = [remove,remove_temp];    %#ok<AGROW>
    end
    C{i} = centroid(2:q);
    T{i} = thresh(2:q);
    
    keep = setxor(1:n,remove);
    regress{i} = false(n,1);
    regress{i}(keep) = true;
    
    % Multivariate Regression: Try Robust Regression first, then Standard Regression if necessary
    try
        BP(:,i) = robustfit(X(regress{i},:),CDM(regress{i},i),[],[],'off');
    catch  %#ok<CTCH>
        ['Species ', num2str(i), ': Robustfit Exception - Using Standard MVR'] %#ok<NOPRT>
        BP(:,i) = (X(regress{i},:)'*X(regress{i},:))\(X(regress{i},:)'*CDM(regress{i},i));
    end
    
    rm = isnan(BP(:,i));
    BP(rm,i) = 0;
    
    Yhat(:,i) = X*BP(:,i);
    rmneg = Yhat(:,i) <= 0;
    Yhat(rmneg,i) = 0;
    Yres(:,i) = CDM(:,i)-Yhat(:,i);

    % Compute the R2 Value
    mutual = CDM(:,i)~=0 & Yhat(:,i)~=0;
    if sum(mutual) <= 3
        AR(i) = 0;
    else
        AR(i) = corr(Yhat(mutual,i),CDM(mutual,i)).^2;
    end
    
    if isnan(AR(i)) == 1
        AR(i) = 0;
    end
end
    
ind = AR>0;
AVGAR = mean(AR(ind));
if isnan(AVGAR)
    AVGAR = 0;
end

ARC = AR;
AVGARC = AVGAR;