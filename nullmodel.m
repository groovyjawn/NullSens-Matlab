function RANDYres = nullmodel(Yres,regress)

% This function randomizes the residuals of CDM as required for generating 
% the null distributuion. Species abundances are randomized using the sites
% selected for each species respectively.

% INPUTS: 
% Yres: residual species responses obtained after fitting CDM
% regress: array of sites selected for analysis, per species

% OUTPUT:
% RANDYres: randomized matrix of species response residuals

[~,p]=size(Yres); % Get number of species p
RANDYres = Yres; % Initialize random matrix with residuals of CDM

% Per species, randomize the abundances at sites selected (regress)
for i = 1:p
    permNon = find(regress{i});
    neworderNon = randsample(permNon,length(permNon));
    RANDYres(permNon,i) = Yres(neworderNon,i);
end