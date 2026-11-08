function export_figure(fig, filename, varargin)
% EXPORT_FIGURE  Export figure as PDF and/or PNG with consistent settings

    p = inputParser;
    addParameter(p, 'Resolution', 600);
    addParameter(p, 'ExportPDF', false);
    addParameter(p, 'ExportPNG', true);
    addParameter(p, 'OutputFolder', fullfile(pwd, 'figures'));
    parse(p, varargin{:});
    s = p.Results;

    set(fig, 'Color', 'w');

    if ~exist(s.OutputFolder, 'dir')
        mkdir(s.OutputFolder);
    end

    png_path = fullfile(s.OutputFolder, [filename '.png']);
    pdf_path = fullfile(s.OutputFolder, [filename '.pdf']);

    if s.ExportPDF
        exportgraphics(fig, pdf_path, 'ContentType', 'vector');
    end

    if s.ExportPNG
        exportgraphics(fig, png_path, 'Resolution', s.Resolution);
    end
end