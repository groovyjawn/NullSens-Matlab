function RANDYres = nullmodel(Yres,regress)

[~,p]=size(Yres);
RANDYres = Yres;

for i = 1:p
    permNon = find(regress{i});
    neworderNon = randsample(permNon,length(permNon));
    RANDYres(permNon,i) = Yres(neworderNon,i);
end