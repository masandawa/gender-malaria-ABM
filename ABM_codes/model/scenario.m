function [P, C] = scenario(P, bits)
%SCENARIO  Mechanism configuration used in all analyses.
%
%   bits = [B A G S] switches four mechanism blocks on (1) or off (0):
%     B  gendered behaviour : female-male activity contrast delta_A = -0.15,
%                             habitual ITN use 0.865 (F) vs 0.835 (M),
%                             female effective treatment 0.54 vs 0.60 (M)
%     A  programme access   : policy-aligned ITN targets (78% pregnant women and
%                             children <5, 68% others) with ANC top-up, and IPTp
%     G  pregnancy biology  : attractiveness x2.0, clinical x1.5, severe x3.0
%     S  sex-specific clearance of chronic infection (F:M hazard ratio 1.82)
%   bits = [1 1 1 1] is the reference configuration; [0 0 0 0] is equality.
%   With A = 0, 70% of each sex receive an ITN.
%
%   The string 'REF' is accepted as shorthand for [1 1 1 1].
if ischar(bits) && strcmpi(bits, 'REF'), bits = [1 1 1 1]; end
C.alloc = 'sex'; C.cov = 0.70; C.ancTopUp = 0; C.realloc = 365;
C.covPol = [0.78 0.78 0.68];
% equality baseline
P.deltaA = 0; P.piUf = P.piU; P.piUm = P.piU; P.ptxF = P.ptx; P.ptxM = P.ptx;
P.pregBio = 0; P.sexClear = 0; P.iptpOn = 0;
if bits(1)
    P.deltaA = -0.15; P.piUf = 0.865; P.piUm = 0.835; P.ptxF = 0.90*P.ptx;
end
if bits(2)
    C.alloc = 'policyProb'; C.ancTopUp = 0.31; P.iptpOn = 1;
end
if bits(3), P.pregBio = 1; end
if bits(4), P.sexClear = 1; end
end
