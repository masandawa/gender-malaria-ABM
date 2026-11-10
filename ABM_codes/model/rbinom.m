function k = rbinom(n, p)
%RBINOM  Binomial(n,p) draw without toolbox dependencies.
n = round(n);
if n <= 0 || p <= 0, k = 0; return; end
if p >= 1, k = n; return; end
v = n*p*(1-p);
if v > 25
    k = round(n*p + sqrt(v)*randn);           % normal approximation
elseif n <= 5000
    k = sum(rand(n,1) < p);                    % exact
else
    k = rpois(n*p);                            % Poisson approximation (small p)
end
k = max(0, min(n, k));
end
