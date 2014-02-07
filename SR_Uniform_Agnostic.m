function [CDM,X,cvid,Y,YN,Mag,CDMB,B] = SR_Uniform_Agnostic(n,p,q,N,M,CD,Type,Noise,Magnitude)

% CDM:  Community Data Matrix (Sites x Species)
% X:    Environmental Data (Explanatory Variables) (Includes 1's intercept)
% Y:    Environmental Species Responses
% YN:   Noise (Additive)
% Mag:  Magnitude of Species Covariation
% CDMB: Community Data Matrix prior to Truncation
% B:    Species Response Parameters
% cvid: Type of species covariation for mixed dynamics model

X = -1 + (1-(-1)).*rand(n,q);
X = [ones(n,1),X];
B = zeros(q+1,p);
Y = zeros(n,p);
YN = zeros(n,q);
CDM = zeros(n,p);
cvid = zeros(1,p);

for i = 1:2:CD*2
    while length(find(CDM(:,i)>0 & CDM(:,i+1)>0)) <= 0.05*n || length(find(CDM(:,i)>0 & CDM(:,i+1)>0)) <= 7
        B(:,i:i+1) = -1 + (1-(-1)).*rand(q+1,2);
        Y(:,i:i+1) = X*B(:,i:i+1);
        YN(:,i:i+1) = randn(n,2)*Noise(N);
        CDM(:,i:i+1) = Y(:,i:i+1)+YN(:,i:i+1);
        Mag = Magnitude(M)*(-1 + (1-(-1)).*rand(1,n)');
        CDM(:,i) = CDM(:,i)+Mag;
        if Type == 1 % Positive Dynamics
            CDM(:,i+1) = CDM(:,i+1)+Mag; 
        elseif Type == 0 % Negative Dynamics
            CDM(:,i+1) = CDM(:,i+1)-Mag;
        elseif Type == 2 % Mixed Dynamics
            flip = randi([0,1],1);
            cvid(i:i+1) = flip;
            if flip == 1 % Positive Dynamics
                CDM(:,i+1) = CDM(:,i+1)+Mag;
            else % Negative Dynamics
                CDM(:,i+1) = CDM(:,i+1)-Mag;
            end
        end
    end
end

for i = CD*2+1:p
    while length(find(CDM(:,i)>0)) <= 0.05*n || length(find(CDM(:,i)>0)) <= 7
        B(:,i) = -1 + (1-(-1)).*rand(q+1,1);
        Y(:,i) = X*B(:,i);
        YN(:,i) = randn(n,1)*Noise(N);
        CDM(:,i) = Y(:,i)+YN(:,i);
    end
end

CDMB = CDM; % CDM before Zero Truncation
rmneg = CDM <= 0;
CDM(rmneg) = 0; % Community Data Matrix Truncated 