function mlegend(ax, labels, styles, xy, dy)
%MLEGEND  Compact legend drawn inside axes at normalised position xy = [x y]
%   (top-left corner), with row spacing dy (normalised). styles is a cell array
%   of {colour, marker, linestyle}. Renders identically in MATLAB and Octave.
if nargin < 5, dy = 0.085; end
xl = xlim(ax); yl = ylim(ax);
lx = strcmp(get(ax, 'XScale'), 'log'); ly = strcmp(get(ax, 'YScale'), 'log');
tox = @(u) map(u, xl, lx); toy = @(v) map(v, yl, ly);
for i = 1:numel(labels)
    c = styles{i}{1}; mk = styles{i}{2}; lst = styles{i}{3};
    yv = toy(xy(2) - (i - 1)*dy);
    if ~strcmp(lst, 'none')
        plot(ax, [tox(xy(1)) tox(xy(1) + 0.06)], [yv yv], lst, 'Color', c, 'LineWidth', 1.3);
    end
    if strcmp(mk, 'bar')
        plot(ax, tox(xy(1) + 0.03), yv, 's', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 6);
    elseif ~strcmp(mk, 'none')
        plot(ax, tox(xy(1) + 0.03), yv, mk, 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 4.5);
    end
    text(tox(xy(1) + 0.08), yv, labels{i}, 'FontSize', 7.5, 'VerticalAlignment', 'middle', 'Parent', ax);
end
end

function v = map(u, lim, islog)
if islog, v = 10.^(log10(lim(1)) + u*(log10(lim(2)) - log10(lim(1))));
else, v = lim(1) + u*(lim(2) - lim(1)); end
end
