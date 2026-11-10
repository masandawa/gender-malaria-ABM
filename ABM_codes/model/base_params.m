function P = base_params(X, m)
%BASE_PARAMS  Reference parameters with run-size settings X and density m.
P = abm_params();
P.N = X.N;
if nargin > 1 && ~isempty(m), P.m = m; end
end
