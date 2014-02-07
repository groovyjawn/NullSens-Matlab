function RESULTS = NullSens(n,p,q,N,M,CD,Type,Name)

tStart = tic;
Noise = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2];
Magnitude = Noise(N).*[1,3,5,10,20];
cvid = zeros(100,p);
AR = zeros(100,p);
AVGAR = zeros(100,1);
ARC = zeros(100,p);
AVGARC = zeros(100,1);
AvgRSq = zeros(1,100);
PRSqYres = zeros(1,100);
CR = cell(1,100);
CV = cell(1,100);
invalid = zeros(1,100);

for iter = 1:100
    clearvars -except tStart iter n p q N M CD Type Name Noise Magnitude cvid AR AVGAR ARC AVGARC AvgRSq AvgRSqYres PRSqYres CR CV invalid RESULTS
    
    [CDM, X, cvid] = SR_Uniform_Agnostic(n,p,q,N,M,CD,Type,Noise,Magnitude);
    
    [~, Yres, ~, regress, AR(iter,:), ARC(iter,:), AVGAR(iter), AVGARC(iter)] = mvrCombined(CDM,X);

    reps = 200; % Number of RAND matrices generated in nullmodel
    AvgRSqYres = zeros(1,reps);
    for i = 1:reps-1
        RANDYres = nullmodel(Yres,regress); % Nullmodel - generates RAND matrix
        AvgRSqYres(i) = COOV(RANDYres,regress,Type,CD); % Index on residuals of CDM 
    end
    
    % Index computation on CDM under test
    [AvgRSqYres(i+1), invalid(iter), CR{iter}, CV{iter}] = COOV(Yres,regress,Type,CD); % Index on redsiduals of CDM
    AvgRSq(iter) = AvgRSqYres(reps);

    % Get P-values - Determine if residuals of CDM is significant 
    indRYres = find(AvgRSqYres >= AvgRSqYres(reps));
    PRSqYres(iter) = length(indRYres)/reps;   
end

tStop = toc(tStart);
RESULTS =  {[n,p,q,N,M,CD,Type],cvid,AR,ARC,[AVGAR,AVGARC],AvgRSq,PRSqYres,CR,CV,invalid,tStop};
save(Name,'RESULTS')