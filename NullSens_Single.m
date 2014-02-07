delta  = [10,20,30,40,50,60,70,80,90,100,110,120];
for repetition = 1:length(delta)
n = 100;
p = delta(repetition);
q = 2;
N = 4;
M = 1;
CD = 0;
Type = 0;

tStart = tic;
Noise = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2];
Magnitude = Noise(N).*[1,3,5,10,20];
cvid = zeros(50,p);
AR = zeros(50,p);
AVGAR = zeros(50,1);
ARC = zeros(50,p);
AVGARC = zeros(50,1);
AvgRSq = zeros(1,50);
PRSqYres = zeros(1,50);
CR = cell(1,50);
CV = cell(1,50);
invalid = zeros(1,50);

for iter = 1:1
    clearvars -except delta TS repetition tStart iter n p q N M CD Type Name Noise Magnitude cvid AR AVGAR ARC AVGARC AvgRSq AvgRSqYres PRSqYres CR CV invalid RESULTS
    [CDM, X, cvid] = SR_Uniform_Agnostic(n,p,q,N,M,CD,Type,Noise,Magnitude);
    %for i = 1:n
    %    s = sum(CDM(i,:));
    %    CDM(i,:) = sqrt(CDM(i,:)./s);
    %end
    [~, Yres, ~, regress, AR(iter,:), ARC(iter,:), AVGAR(iter), AVGARC(iter)] = mvrCombined(CDM,X); %mvrCombined(CDM,X);

    reps = 2; % Number of RAND matrices generated in nullmodel
    AvgRSqYres = zeros(1,reps);
    ts = tic;
    for i = 1:reps-1
        RANDYres = nullmodel(Yres,regress); % Nullmodel - generates RAND matrix
        AvgRSqYres(i) = COOV(RANDYres,regress,Type,CD); % Index on residuals of CDM
    end
    TS(repetition) = toc(ts);

    % Index computation on CDM under test
    [AvgRSqYres(i+1), invalid(iter), CR{iter}, CV{iter}] = COOV(Yres,regress,Type,CD); % Index on redsiduals of CDM
    AvgRSq(iter) = AvgRSqYres(reps);

    % Get P-values - Determine if residuals of CDM is significant 
    indRYres = find(AvgRSqYres >= AvgRSqYres(reps));
    PRSqYres(iter) = length(indRYres)/reps;
end

tStop = toc(tStart);
RESULTS =  {[n,p,q,N,M,CD,Type],cvid,AR,ARC,[AVGAR,AVGARC],AvgRSq,PRSqYres,CR,CV,invalid,tStop};
end