function VisualizeNetwork(obj, Atoms, Bonds, Nvec)
% -------------------------------------------------------------------------
% VisualizeNetwork
% - Visualize generated network and optional distributions.
%
% Atoms column layout: [ID | molID | X | Y | Z | deg | nbrs...]
% X at col 3, Y at col 4.
% -------------------------------------------------------------------------

    if ~obj.flags.iplot
        obj.log.print('   Skipping visualization because obj.flags.iplot = false.\n');
        return;
    end

    obj.log.print('   Visualizing network...\n');

    N_atom = size(Atoms,1);
    Total_bond = size(Bonds,1);

    if N_atom == 0
        obj.log.print('   Visualization skipped because Atoms is empty.\n');
        return;
    end

    mode = lower(obj.architecture.strand_typology.mode);
    b = obj.domain.b;
    scale = obj.domain.scale;
    if isempty(scale)
        scale = 1;
    end

    % ---------------------------------------------------------------------
    % Main network plot
    % ---------------------------------------------------------------------
    figure; hold on;

    scatter(Atoms(:,3), Atoms(:,4), 8, 'k', 'filled');

    for k = 1:Total_bond
        if Bonds(k,1) == 0
            continue;
        end

        i1 = Bonds(k,2);
        i2 = Bonds(k,3);

        if i1 < 1 || i1 > N_atom || i2 < 1 || i2 > N_atom
            continue;
        end

        if strcmp(mode, 'bimodal') && size(Bonds,2) >= 5 && Bonds(k,5) ~= 1
            plot([Atoms(i1,3) Atoms(i2,3)], ...
                 [Atoms(i1,4) Atoms(i2,4)], ...
                 'r-', 'LineWidth', 1.5);
        else
            plot([Atoms(i1,3) Atoms(i2,3)], ...
                 [Atoms(i1,4) Atoms(i2,4)], ...
                 'k-');
        end
    end

    axis equal tight;
    title('Final bonds (post-cleanup)');
    set(gca, 'FontSize', 14, 'LineWidth', 1.5);
    hold off;

    % ---------------------------------------------------------------------
    % Optional histogram blocks
    % ---------------------------------------------------------------------
    if isempty(Bonds) || isempty(Nvec)
        obj.log.print('   Skipping histogram plots because Bonds or Nvec is empty.\n');
        return;
    end

    if any(strcmp(mode, {'poly', 'polydisperse'}))

        nbinsN = max(10, min(80, ceil(sqrt(numel(Nvec)))));
        figure; hold on;
        histogram(Nvec, nbinsN, ...
            'FaceColor', [0.2 0.2 0.2], ...
            'FaceAlpha', 0.8, ...
            'LineWidth', 0.0005);
        xlabel('Assigned N per bond');
        ylabel('Count');
        title('N distribution (polydisperse)');
        set(gca, 'FontSize', 14, 'LineWidth', 1.5);
        hold off;

        figure; hold on;
        histogram(Bonds(:,4), 60, ...
            'FaceColor', [0.1 0.1 0.9], ...
            'FaceAlpha', 0.8, ...
            'LineWidth', 0.0005);
        xlabel('Bond length L');
        ylabel('Count');
        title('Length distribution (polydisperse)');
        axis tight;
        set(gca, 'FontSize', 14, 'LineWidth', 1.5);
        hold off;

        lamvec = (Bonds(:,4)) ./ (Nvec * b);

        dlam = 0.01;
        lam_min = min(lamvec);
        lam_max = max(lamvec);
        nbinsLam = max(10, min(100, ceil((lam_max - lam_min) / dlam)));

        figure; hold on;
        histogram(lamvec, nbinsLam, ...
            'FaceColor', [0.0 0.6 0.0], ...
            'FaceAlpha', 0.8, ...
            'LineWidth', 0.0005);

        xlabel('Prestretch \lambda = L/(N b)');
        ylabel('Count');
        title('Prestretch distribution (polydisperse)');
        axis tight;
        set(gca, 'FontSize', 14, 'LineWidth', 1.5);

        yl = ylim;
        lam_med = median(lamvec);
        lam_ref = 1 / sqrt(mean(Nvec));

        plot([lam_med lam_med], yl, '-', ...
            'LineWidth', 1.5, 'Color', [0.3 0.3 0.3]);
        plot([lam_ref lam_ref], yl, '--', ...
            'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);

        legend('All bonds', 'median \lambda', '1/\surd(mean N)', ...
            'Location', 'best');
        hold off;

        figure; hold on;

        nbinsN3D = 40;
        nbinsLam3D = 40;

        edgesN = linspace(min(Nvec), max(Nvec), nbinsN3D + 1);
        edgesLam = linspace(min(lamvec), max(lamvec), nbinsLam3D + 1);

        [counts, edgesN, edgesLam] = histcounts2(Nvec, lamvec, edgesN, edgesLam);

        centersN = 0.5 * (edgesN(1:end-1) + edgesN(2:end));
        centersLam = 0.5 * (edgesLam(1:end-1) + edgesLam(2:end));

        surf(centersN, centersLam, counts', ...
            'EdgeColor', 'none', 'FaceColor', 'interp');

        colormap(parula);
        colorbar;
        xlabel('Kuhn segments N');
        ylabel('Prestretch \lambda = L/(N b)');
        zlabel('Count');
        title('Joint distribution of N and prestretch');
        set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'View', [45 30]);
        grid on;
        box on;
        hold off;

    elseif strcmp(mode, 'bimodal')

        if size(Bonds,2) >= 5
            type1 = (Bonds(:,5) == 1);
        else
            type1 = true(size(Bonds,1),1);
        end

        nbins = max(10, min(80, ceil(sqrt(numel(Nvec)))));
        figure; hold on;
        histogram(Nvec(type1), nbins, ...
            'FaceColor', [0.2 0.2 0.2], ...
            'FaceAlpha', 0.8, ...
            'LineWidth', 0.0005);
        histogram(Nvec(~type1), nbins, ...
            'FaceColor', [1 0 0], ...
            'FaceAlpha', 0.8, ...
            'LineWidth', 0.0005);
        xlabel('Assigned N per bond');
        ylabel('Count');
        title('N distribution (bimodal)');
        set(gca, 'FontSize', 14, 'LineWidth', 1.5);
        hold off;

        figure; hold on;
        histogram(Bonds(type1,4), 50, ...
            'FaceColor', [0.2 0.2 0.2], ...
            'FaceAlpha', 1.0, ...
            'LineWidth', 0.0005);
        histogram(Bonds(~type1,4), 50, ...
            'FaceColor', [1 0 0], ...
            'FaceAlpha', 1.0, ...
            'LineWidth', 0.0005);
        xlabel('Bond length L');
        ylabel('Count');
        title('Length distribution (bimodal)');
        axis tight;
        set(gca, 'FontSize', 14, 'LineWidth', 1.5);
        hold off;

        lamvec1 = (Bonds(type1,4)) ./ (Nvec(type1) * b);
        lamvec2 = (Bonds(~type1,4)) ./ (Nvec(~type1) * b);

        N1 = obj.perbond.kuhn.bimodal.mean_1;
        N2 = obj.perbond.kuhn.bimodal.mean_2;

        figure; hold on;

        if sum(type1) > 1
            dlam1 = max(1, ceil((max(lamvec1) - min(lamvec1)) / 0.01));
            histogram(lamvec1, dlam1, ...
                'FaceColor', [0.2 0.2 0.2], ...
                'FaceAlpha', 0.8, ...
                'LineWidth', 0.0005);
        end

        if sum(~type1) > 1
            dlam2 = max(1, ceil((max(lamvec2) - min(lamvec2)) / 0.01));
            histogram(lamvec2, dlam2, ...
                'FaceColor', [1 0 0], ...
                'FaceAlpha', 0.8, ...
                'LineWidth', 0.0005);
        end

        axis tight;
        yl = ylim;

        plot([1/sqrt(N1) 1/sqrt(N1)], yl, '--', ...
            'LineWidth', 2, 'Color', [0.2 0.2 0.2]);
        plot([1/sqrt(N2) 1/sqrt(N2)], yl, '--', ...
            'LineWidth', 2, 'Color', [1 0 0]);

        xlabel('Prestretch \lambda = L/(N b)');
        ylabel('Count');
        title('Prestretch distribution (bimodal)');
        set(gca, 'FontSize', 14, 'LineWidth', 1.5);
        hold off;

    else
        obj.log.print('   No additional histograms defined for mode "%s".\n', mode);
    end

end
