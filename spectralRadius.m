function s = spectralRadius(A)

if sum(sum(isnan(A)))
    s = nan;
else
    tmp = eig(A);
    s = max(sqrt(diag(tmp*tmp')));
end