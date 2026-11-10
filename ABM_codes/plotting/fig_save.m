function fig_save(f, name)
%FIG_SAVE  Save a figure as 600-dpi PNG and vector PDF in figures/.
png = fullfile('figures', [name '.png']); pdf = fullfile('figures', [name '.pdf']);
if exist('OCTAVE_VERSION', 'builtin')
    print(f, png, '-dpng', '-r600');
    print(f, pdf, '-dpdf');
else
    exportgraphics(f, png, 'Resolution', 600);
    exportgraphics(f, pdf, 'ContentType', 'vector');
end
close(f);
end
