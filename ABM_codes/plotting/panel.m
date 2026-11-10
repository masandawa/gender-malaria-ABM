function ax = panel(pos, letter, ttl)
%PANEL  Axes at normalised position pos with a bold panel letter and title.
ax = axes('Position', pos); hold(ax, 'on');
set(ax, 'TickDir', 'out', 'Box', 'off', 'FontSize', 8.5, 'LineWidth', 0.6);
title(ax, sprintf('(%s) %s', letter, ttl), 'FontSize', 9, 'FontWeight', 'bold');
end
