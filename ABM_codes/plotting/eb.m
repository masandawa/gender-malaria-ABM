function h = eb(x, y, lo, hi, fmt, varargin)
%EB  Error bars from interval limits (MATLAB and Octave compatible).
h = errorbar(x, y, y - lo, hi - y, fmt);
if ~isempty(varargin), set(h, varargin{:}); end
end
