function [vector] = extractCorrelations(RSAmatrix, labels, trialtype1, trialtype2)
% extract correlations from a RSA matrix matching these two trial types

vector = [];

[rs, cs] = size(RSAmatrix);

for r = 1:rs
    for c = 1:cs
        if strcmp(labels{r}, trialtype1) && ...
                strcmp(labels{c}, trialtype2) && ...
                RSAmatrix(r,c) ~= 0
            vector = horzcat(vector, RSAmatrix(r,c));
        end
    end
end

end