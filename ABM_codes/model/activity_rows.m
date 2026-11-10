function [Pout, Bed] = activity_rows(P, S, idx)
%ACTIVITY_ROWS  Hourly outdoor probabilities and in-bed masks for agents idx.
%   logit p_out,i(h) = logit p_base(h; age) + s_i*delta_A/2 + re_i,
%   with s_i = +1 for females and -1 for males, so delta_A < 0 means that
%   males spend more of the biting period outdoors.
if nargin < 3, idx = (1:P.N)'; end
n   = numel(idx);
Pout = zeros(n, 14); Bed = false(n, 14);
young = S.age(idx) < 5;
lA = log(P.poutA./(1-P.poutA)); lY = log(P.poutY./(1-P.poutY));
base = repmat(lA, n, 1); base(young,:) = repmat(lY, sum(young), 1);
sgn  = 1 - 2*S.sex(idx);                         % +1 female, -1 male
lo   = base + repmat(sgn*P.deltaA/2 + S.re(idx), 1, 14);
Pout = 1./(1+exp(-lo));
Bed  = repmat(P.bedA, n, 1); Bed(young,:) = repmat(P.bedY, sum(young), 1);
end
