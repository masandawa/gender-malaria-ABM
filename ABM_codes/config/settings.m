function X = settings(mode)
%SETTINGS  Run sizes for the two modes of run_all.m.
%
%   'quick' : small populations and few realisations; checks that the whole
%             pipeline works and produces every table and figure (minutes).
%   'paper' : the settings used for the published results.
if nargin < 1, mode = 'paper'; end
X.mode = mode;
switch lower(mode)
    case 'paper'
        X.N = 1000;
        X.burn = 730; X.eval = 1825; X.hist = 365;
        X.calGrid = 10:2:20;  X.calReps = 8;           % A0 calibration
        X.profReps = 4;                                 % A0 profile
        X.A1reps = 12;                                  % behaviour
        X.A2reps = 16;                                  % factorial
        X.A3reps = 10;                                  % pregnancy sensitivity
        X.A4reps = 12;                                  % ITN use mode
        X.A5reps = 40;                                  % allocation
        X.A6reps = 8;  X.Rgrid = [0.70 1.43 2.50 4.37]; % biting ratio
        X.A7draws = 40; X.A7calReps = 2;                % Latin hypercube (m re-calibrated per draw)
        X.A9reps = 12; X.A9calReps = 4;                 % mechanism checks (Experiment 6)
    case 'quick'
        X.N = 400;
        X.burn = 365; X.eval = 730; X.hist = 365;
        X.calGrid = [12 16 20]; X.calReps = 2;
        X.profReps = 1;
        X.A1reps = 2; X.A2reps = 2; X.A3reps = 2; X.A4reps = 2; X.A5reps = 3;
        X.A6reps = 2; X.Rgrid = [0.70 4.37];
        X.A7draws = 4; X.A7calReps = 1;
        X.A9reps = 2; X.A9calReps = 1;
    otherwise
        error('mode must be ''quick'' or ''paper''');
end
end
