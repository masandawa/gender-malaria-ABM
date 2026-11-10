%RUN_ALL  Reproduce every number, table and figure of the manuscript.
%
%   1. Set MODE below:
%        'quick' - small runs (minutes) to check that everything works;
%        'paper' - the settings used for the published results.
%   2. Run this script from the repository folder:  >> run_all
%
%   Each analysis writes CSV files to results/. Tables and figures are then
%   produced from those CSV files only, so they can be re-drawn instantly with
%        >> make_tables; make_figures
%   without re-running any simulation.
%
%   Set any RUN.Ax flag to false to skip an analysis whose results already exist.
%   Works in MATLAB (R2019b or later, no toolboxes) and GNU Octave (>= 6).

clear; clc;
MODE = 'paper';                       % 'quick' or 'paper'

RUN.A0 = true;   % calibration and identifiability profile
RUN.A1 = true;   % gendered night-time activity, event concentration
RUN.A2 = true;   % 2^4 factorial decomposition of mechanisms
RUN.A3 = true;   % pregnancy-parameter sensitivity
RUN.A4 = true;   % habitual vs nightly ITN use
RUN.A5 = true;   % fixed-resource ITN recipient selection
RUN.A6 = true;   % outdoor:indoor biting ratio (fixed and re-calibrated)
RUN.A7 = true;   % joint parameter uncertainty (Latin hypercube)
RUN.A8 = true;   % verification tests
RUN.A9 = true;   % mechanism checks: immunity structure, ITN protection components

root = fileparts(mfilename('fullpath')); if isempty(root), root = pwd; end
cd(root);
addpath(fullfile(root, 'model'), fullfile(root, 'analyses'), ...
        fullfile(root, 'plotting'), fullfile(root, 'config'));
if ~exist('results', 'dir'), mkdir('results'); end
if ~exist('figures', 'dir'), mkdir('figures'); end
X = settings(MODE);
t0 = tic;

% ---------------------------------------------------------------- A0
if RUN.A0
    m = A0_calibration(X);
else
    S = csv2struct('results/A0_calibration_summary.csv'); m = S.m_star;
end
fprintf('Calibrated mosquito-to-human ratio m* = %.2f\n', m);

if RUN.A1, A1_behaviour(X, m); end
if RUN.A2, A2_factorial(X, m); end
if RUN.A3, A3_pregnancy_sensitivity(X, m); end
if RUN.A4, A4_itn_use(X, m); end
if RUN.A5, A5_allocation(X, m); end
if RUN.A6, A6_biting_ratio(X, m); end
if RUN.A7, A7_uncertainty(X, m); end
if RUN.A8, A8_verification(X, m); end
if RUN.A9, A9_mechanism_checks(X, m); end

make_tables();
make_figures();
fprintf('Finished in %.1f minutes. See results/ and figures/.\n', toc(t0)/60);
