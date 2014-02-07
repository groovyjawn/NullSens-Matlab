function [AvgRSq, invalid, CR, CV] = COOV(matrix,regress,Type,CD)

[n,p] = size(matrix);

CR = zeros(p-1,p);
CV = zeros(p-1,p);
w = 0;
AvgRSq = 0;
invalid = 0;
for i = 1:p-1
    for j = i+1:p        
        temp = and(regress{i},regress{j});
        if sum(temp) <= n*0.05 || sum(temp) <= 7
            cvar = 0;
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
            COV = cov(matrix(temp,[i,j]));
            CORR = corr(matrix(temp,[i,j]));
            if isnan(COV(1,2)) == 1
                cvar = 0;
                CR(i,j) = 0;
                CV(i,j) = 0;
            else
                cvar = COV(1,2);
                CR(i,j) = CORR(1,2);
                CV(i,j) = cvar;
            end 
        end
        
        AvgRSq = AvgRSq + cvar.^2;
        w = w + abs(cvar);
    end
end

if w == 0
    AvgRSq = 0;
else
    AvgRSq = AvgRSq/w;
end