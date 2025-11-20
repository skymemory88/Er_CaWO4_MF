function [T_top, idxTop] = plot_S1S2_top_zefoz(pertResults, const, Options, varargin)
% plot_S1S2_top_zefoz_compact
%   Visualize S1 and S2 for the top ZEFOZ candidates with:
%     - Consistent 2×2 grid layout (4 subplots per figure), paging as needed,
%     - compact per-candidate panel: S1 (bars, left y) + S2 (markers, right y),
%     - optional separate tensor-heatmap figures (4 heatmaps per figure),
%     - robust tensor scaling to avoid "solid color" panels.
%
% Inputs
%   pertResults : struct containing outputs of perturb2nd (tensor and/or axis fields)
%   const       : constants with const.Gh2mV (meV per GHz)
%   Options     : struct for compute_all_S1S2 (supports Options.ndE or Options.pairs)
%
% Name-Value pairs
%   'TopN'         : number of top ZEFOZ candidates to visualize (default 8, optimized for 2×2 layouts)
%   'PerFigure'    : candidates per figure (fixed at 4 for consistent 2×2 layouts)
%   'ShowTensor'   : true/false to also plot S2 tensor heatmaps in separate figures (default false)
%   'HeatmapScale' : 'panel' (default) | 'global' — color limits per heatmap or shared
%   'NoiseCov'     : scalar sigma_B (T) or 3×3 covariance (T^2); if provided and tensor
%                    exists, predicted quasi-static linewidth σ_ν is annotated
%
% Outputs
%   T_top  : table of the plotted candidates (subset of compute_all_S1S2 output)
%   idxTop : indices into the full, sorted table returned by compute_all_S1S2
%
% Notes
%   * S2 shown is the TRUE curvature (Hessian), i.e., 2× coefficient already.
%   * Axis-wise S1/S2 always shown; full tensor only when available.

    % ---------- Parse options ----------
    p = inputParser;
    addParameter(p,'TopN',8);  % Default 8 for consistent 2×2 layouts
    addParameter(p,'PerFigure',4);
    addParameter(p,'ShowTensor',false);
    addParameter(p,'HeatmapScale','panel');  % 'panel' or 'global'
    addParameter(p,'NoiseCov',[]);
    parse(p,varargin{:});
    TopN       = p.Results.TopN;
    PerFigure  = 4; % Fixed at 4 for consistent 2×2 layouts
    ShowTensor = logical(p.Results.ShowTensor);
    HeatScale  = validatestring(p.Results.HeatmapScale, {'panel','global'});
    NoiseCov   = p.Results.NoiseCov;

    % ---------- Compute S1/S2 for requested transitions ----------
    % Optimize Options for memory efficiency when only plotting top candidates
    if TopN < 1000 && ~isfield(Options, 'maxTransitions')
        Options.maxTransitions = min(5000, TopN * 50); % Reasonable limit for plotting
    end
    
    % Pre-filter for ZEFOZ if plotting only top candidates (optimized for TopN=8)
    if TopN <= 32 && ~isfield(Options, 'onlyZefoz')  % Optimized for multiples of 4
        Options.onlyZefoz = true;
        % Use Options.zefoz(3) as the criteria, with fallback
        if isfield(Options, 'zefoz') && length(Options.zefoz) >= 3
            Options.zefozThreshold = abs(Options.zefoz(3)) * 1000; % Convert GHz/mT to GHz/T
        else
            Options.zefozThreshold = 1e-2; % Fallback threshold for initial filtering
        end
    end
    
    [T_all, Trans] = compute_all_S1S2(pertResults, const, Options);
    if isempty(T_all) || height(T_all)==0
        error('No transitions produced by compute_all_S1S2.');
    end
    nPlot = min(TopN, height(T_all));
    idxTop = (1:nPlot).';
    T_top  = T_all(idxTop,:);
    
    % ---------- Map table rows to Trans entries ----------
    mapIdx = zeros(nPlot,1);
    for r = 1:nPlot
        m = T_top.m(r); n = T_top.n(r); ord = T_top.order(r);
        mapIdx(r) = find([Trans.m].'==m & [Trans.n].'==n & [Trans.order].'==ord, 1, 'first');
    end
    
    % Memory cleanup: clear large arrays if not needed
    if nPlot < height(T_all) / 10
        clear T_all; % Free memory for large unused data
    end
    
    % Display S1, max(S2) pairs for identified ZEFOZ modes
    fprintf('\n         ZEFOZ candidates = %d\n', nPlot);
    fprintf('Mode  |m⟩→|n⟩   Freq(GHz)   |S1|(GHz/T)   max|S2|(GHz/T²)\n');
    fprintf('----  --------   ---------   -----------   ---------------\n');
    for r = 1:nPlot
        freqGHz = T_top.freq_GHz(r);
        s1norm = T_top.S1norm_GHzT(r);
        % Calculate max|S2| from individual components
        tr_idx = mapIdx(r);
        s2_vals = [Trans(tr_idx).s2X, Trans(tr_idx).s2Y, Trans(tr_idx).s2Z];
        max_s2 = max(abs(s2_vals));
        fprintf('%2d    |%d⟩→|%d⟩   %8.4f   %10.3e   %13.3e\n', ...
                r, Trans(tr_idx).n-1, Trans(tr_idx).m-1, freqGHz, s1norm, max_s2);
    end
    fprintf('========================================\n\n');

    % ---------- Chunk into pages of <=4 panels ----------
    nPages = ceil(nPlot / PerFigure);
    xlbls = {'X','Y','Z'};

    % Precompute noise-aware sigma_nu if requested & tensor available
    wantNoise = ~isempty(NoiseCov) && exist('predict_linewidth_from_S1S2','file')==2;

    for page = 1:nPages
        i1 = (page-1)*PerFigure + 1;
        i2 = min(page*PerFigure, nPlot);
        idxChunk   = idxTop(i1:i2);
        transChunk = mapIdx(i1:i2);
        Kchunk     = numel(idxChunk);

        % Determine y-limits per page for comparability
        S1vals = []; S2vals = [];
        for j = 1:Kchunk
            tr = Trans(transChunk(j));
            if ~isempty(tr.S1_vec_GHzT), s1v = tr.S1_vec_GHzT(:).';
            else, s1v = [tr.s1X, tr.s1Y, tr.s1Z];
            end
            s2v = [tr.s2X, tr.s2Y, tr.s2Z];
            S1vals = [S1vals, s1v];
            S2vals = [S2vals, s2v];
        end
        S1max = max(1e-12, max(abs(S1vals)));
        S2max = max(1e-12, max(abs(S2vals)));

        % ---- MAIN FIGURE: compact S1/S2 panels (consistent 2×2 layout) ----
        fh = figure('Name',sprintf('ZEFOZ: S1-S2 page %d/%d', page, nPages), 'Color','k');
        % Always use 2×2 layout for visual consistency
        rows = 2; cols = 2;
        tl = tiledlayout(rows, cols, 'TileSpacing','compact', 'Padding','compact');

        for j = 1:Kchunk
            r  = idxChunk(j);
            tr = Trans(transChunk(j));
            s1norm  = T_top.S1norm_GHzT(r);

            ax = nexttile(tl);
            % S1 (left y): bars
            if ~isempty(tr.S1_vec_GHzT), s1v = tr.S1_vec_GHzT(:).';
            else, s1v = [tr.s1X, tr.s1Y, tr.s1Z];
            end
            yyaxis(ax,'left');
            bar(1:3, s1v, 0.65, 'EdgeColor','none');
            ylim(ax, [-1.05*S1max, 1.05*S1max]);
            yline(ax,0,'k-');
            ylabel(ax,'S1 (GHz/T)');

            % S2 (right y): markers/line
            s2v = [tr.s2X, tr.s2Y, tr.s2Z];
            yyaxis(ax,'right');
            plot(1:3, s2v, 'o-', 'LineWidth',1);
            ylim(ax, [-1.05*S2max, 1.05*S2max]);
            ylabel(ax,'S2 (GHz/T^2)');

            set(ax, 'XTick',1:3, 'XTickLabel', xlbls);
            box(ax,'on');

            % Title line
            tstr = sprintf('|%d> --> |%d>, |S1|=%.3g', tr.n, tr.m, s1norm);

            % Optional noise-aware linewidth (tensor required)
            if wantNoise && ~isempty(tr.S1_vec_GHzT) && ~isempty(tr.S2_ten_GHzT2)
                [sigmaGHz, muShiftGHz] = predict_linewidth_from_S1S2(tr.S1_vec_GHzT, tr.S2_ten_GHzT2, NoiseCov);
                tstr = sprintf('%s, \\sigma_\\nu≈%.3g kHz (⟨Δν⟩≈%.3g kHz)', tstr, 1e6*sigmaGHz, 1e6*muShiftGHz);
            end
            title(ax, tstr, 'Interpreter','none');
        end

        sgtitle(tl, sprintf('Top ZEFOZ candidates (page %d/%d) — Options.ndE=%s', ...
                            page, nPages, optionsToString(Options)), 'FontWeight','bold');

        % ---- OPTIONAL: separate tensor heatmaps page (<=4 per figure) ----
        if ShowTensor
            % Gather tensors available on this chunk
            haveTen = arrayfun(@(idx) ~isempty(Trans(idx).S2_ten_GHzT2), transChunk);
            if any(haveTen)
                % Compute global heat scale for this page if requested
                if strcmpi(HeatScale,'global')
                    smax = 0;
                    for j = 1:Kchunk
                        if haveTen(j)
                            Kten = Trans(transChunk(j)).S2_ten_GHzT2;
                            smax = max(smax, max(abs(Kten),[],'all'));
                        end
                    end
                    smax = max(smax, 1e-12);
                end

                fh2 = figure('Name',sprintf('ZEFOZ: S2 tensors page %d/%d', page, nPages), 'Color','k');
                tl2 = tiledlayout(2, 2, 'TileSpacing','compact', 'Padding','compact');  % Always 2×2

                for j = 1:Kchunk
                    axh = nexttile(tl2);
                    if haveTen(j)
                        Kten = Trans(transChunk(j)).S2_ten_GHzT2;
                        imagesc(axh, Kten);
                        axis(axh,'image');
                        colorbar(axh);
                        set(axh,'XTick',1:3,'XTickLabel',xlbls,'YTick',1:3,'YTickLabel',xlbls);
                        % Robust scaling to avoid solid color panels
                        if strcmpi(HeatScale,'panel')
                            smaxj = max(1e-12, max(abs(Kten),[],'all'));
                            caxis(axh, [-smaxj, smaxj]);
                        else
                            caxis(axh, [-smax, smax]);
                        end
                        title(axh, sprintf('m=%d, n=%d (order %d)', ...
                              Trans(transChunk(j)).m, Trans(transChunk(j)).n, Trans(transChunk(j)).order));
                    else
                        text(0.5,0.5,'tensor N/A','Parent',axh,'HorizontalAlignment','center');
                        axis(axh,'off');
                    end
                end
                sgtitle(tl2, sprintf('S2 tensors (page %d/%d) — scaling: %s', ...
                                     page, nPages, HeatScale), 'FontWeight','bold');
            end
        end
    end
end

% helper function
function s = optionsToString(Options)
    if isfield(Options,'pairs') && ~isempty(Options.pairs)
        s = 'custom';
    elseif isfield(Options,'ndE')
        if isscalar(Options.ndE), s = sprintf('%d', Options.ndE);
        else, s = sprintf('[%s]', num2str(Options.ndE));
        end
    else
        s = '1';
    end
end
