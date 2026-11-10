function k = rpois(lam)
%RPOIS  Poisson(lam) draw without toolbox dependencies.
if lam <= 0, k = 0; return; end
if lam > 30
    k = max(0, round(lam + sqrt(lam)*randn));  % normal approximation
else
    L = exp(-lam); k = 0; p = rand;
    while p > L, k = k + 1; p = p*rand; end
end
end
