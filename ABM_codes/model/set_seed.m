function set_seed(s)
%SET_SEED  Reproducible random stream in MATLAB and GNU Octave.
if exist('OCTAVE_VERSION', 'builtin')
    rand('twister', s); randn('state', s); %#ok<RAND>
else
    rng(s, 'twister');
end
end
