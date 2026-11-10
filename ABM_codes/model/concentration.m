function Cn = concentration(O)
%CONCENTRATION  Event concentration and clinical-transmission overlap.
%   Uses all persons present during the evaluation period (current occupants
%   and those replaced through demographic turnover).
Nclin = [O.Nclin; O.ret(:,3)]; Z = [O.Zhv; O.ret(:,4)]; Ninf = [O.Ninf; O.ret(:,2)];
n = numel(Z);
Cn.pClin = mean(Nclin > 0); Cn.pTrans = mean(Z > 0);
Cn.inf0 = mean(Ninf==0); Cn.inf1 = mean(Ninf==1); Cn.inf2 = mean(Ninf==2); Cn.inf3p = mean(Ninf>=3);
Cn.meanInf = mean(Ninf);
for x = [0.5 0.8]
    qB = qfrac(Nclin, x, n); qT = qfrac(Z, x, n);
    f = sprintf('%02d', round(100*x));
    Cn.(['qB' f]) = qB; Cn.(['qT' f]) = qT;
    Cn.(['qBt' f]) = qB/max(Cn.pClin, eps); Cn.(['qTt' f]) = qT/max(Cn.pTrans, eps);
end
B = Nclin > 0; T = Z > 0;
Cn.J = sum(B & T)/max(sum(B | T), 1);
end

function q = qfrac(v, x, n)
s = sort(v, 'descend'); c = cumsum(s)/max(sum(s), eps);
k = find(c >= x, 1, 'first'); if isempty(k), k = n; end
q = k/n;
end
