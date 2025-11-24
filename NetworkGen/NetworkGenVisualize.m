function NetworkGenVisualize(Domain,Atoms,Bonds,Nvec,scale,options)
% NetworkGenVisualize - Visualize the generated network
% INPUT: 
%   Domain  : struct with fields xlo, xhi, ylo, yhi, zlo, zhi
%   Atoms   : [ ID | X | Y | Z | num_bond | nbr1 | nbr2 | nbr3 | nbr4 | ... ]
%   Bonds   : [ bondID | id1 | id2 | L | bondType ]
%   options : structure with visualization options (not used currently)
% OUTPUT:
%   Figure showing the network structure
% Assumptions:
%   - 2D network (Z=0) 
%   - R2016a compatible, no implicit expansion

if options.iplot ~= true
    warning('Skipping visualization...');
    return; % skip visualization
end

fprintf('   Visualizing network...\n');

% --------------------- Unpack & Setup ---------------------
N_atom = length(Atoms);
Total_bond = size(Bonds, 1);

isPeriodic  =  strcmpi(options.boundary_box,'Periodic');

% --------------------- Plot final bonds ---------------------
if strcmp(options.dist_type,'bimodal')
    figure; hold on;

    scatter(Atoms(1:N_atom,2), Atoms(1:N_atom,3), 6, 'k', 'filled');
    for k = 1:Total_bond
        if Bonds(k,1) == 0, continue; end
        i1 = Bonds(k,2); i2 = Bonds(k,3);
        if Bonds(k,5) == 1
            plot([Atoms(i1,2) Atoms(i2,2)], [Atoms(i1,3) Atoms(i2,3)], 'Color',[150 150 150]/255,'LineWidth',0.25);
        else
            plot([Atoms(i1,2) Atoms(i2,2)], [Atoms(i1,3) Atoms(i2,3)], 'r-','LineWidth',1.75);
        end
    
    end
    axis equal tight; title('Final bonds (post-prune)');
else
    figure; hold on;
    scatter(Atoms(1:N_atom,2), Atoms(1:N_atom,3), 6, 'k', 'filled');
    for k = 1:Total_bond
        if Bonds(k,1) == 0, continue; end
        i1 = Bonds(k,2); i2 = Bonds(k,3);
        plot([Atoms(i1,2) Atoms(i2,2)], [Atoms(i1,3) Atoms(i2,3)], 'k-');
    end
    axis equal tight; title('Final bonds (post-prune)');
end
% --------------------- Network stats (Kuhn length) ---------------------

% ---------- Optional histogram ----------

if strcmp(options.dist_type,'polydisperse')
    % --- Kuhn-segment distribution ---
    nbinsN = max(10, min(80, ceil(sqrt(numel(Nvec)))));
    figure; hold on
    histogram(Nvec, nbinsN, ...
        'FaceColor', [0.2 0.2 0.2], 'FaceAlpha', 0.8, 'LineWidth', 0.0005);
    xlabel('Assigned N per bond');
    ylabel('Count');
    title('N distribution (polydisperse)');
    set(gca, 'FontSize', 16, 'LineWidth', 2);
    hold off

    % --- Bond-length distribution ---
    figure; hold on
    nbinsL = 60;
    histogram(Bonds(:,4), nbinsL, ...
        'FaceColor', [0.1 0.1 0.9], 'FaceAlpha', 0.8, 'LineWidth', 0.0005);
    axis tight
    xlabel('Bond length L');
    ylabel('Count');
    title('Length distribution (polydisperse)');
    set(gca, 'FontSize', 16, 'LineWidth', 2);
    hold off

    % --- Prestretch distribution ---
    b = options.b;
    lamvec = (Bonds(:,4)/scale) ./ (Nvec * b);   % ? = L/(N·b)
    dlam = 0.01;
    lam_min = min(lamvec);
    lam_max = max(lamvec);
    nbinsLam = max(10, min(100, ceil((lam_max - lam_min)/dlam)));

    figure; hold on
    histogram(lamvec, nbinsLam, ...
        'FaceColor', [0.0 0.6 0.0], 'FaceAlpha', 0.8, 'LineWidth', 0.0005);
    axis tight
    xlabel('Prestretch \lambda = L/(N*b)');
    ylabel('Count');
    title('Prestretch distribution (polydisperse)');
    set(gca, 'FontSize', 16, 'LineWidth', 2);

    % Optional reference lines (R2016-safe)
    yl = ylim;
    lam_med = median(lamvec);
    lam_ref = 1 / sqrt(mean(Nvec));  % approximate reference
    plot([lam_med lam_med], yl, '-',  'LineWidth', 1.5, 'Color', [0.3 0.3 0.3]);
    plot([lam_ref lam_ref], yl, '--', 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
    legend('All bonds','median \lambda','1/\surd(mean N)','Location','best');
    hold off
    
    % --- 3D histogram: joint distribution of N and ? ---
    figure; hold on

    % Define bins for N and ?
    nbinsN3D = 40;
    nbinsLam3D = 40;

    edgesN = linspace(min(Nvec), max(Nvec), nbinsN3D+1);
    edgesLam = linspace(min(lamvec), max(lamvec), nbinsLam3D+1);

    % Compute 2D histogram (count matrix)
    [counts, edgesN, edgesLam] = histcounts2(Nvec, lamvec, edgesN, edgesLam);

    % Convert to bin centers for plotting
    centersN = 0.5*(edgesN(1:end-1) + edgesN(2:end));
    centersLam = 0.5*(edgesLam(1:end-1) + edgesLam(2:end));

    % Create a surface plot of counts
    surf(centersN, centersLam, counts', ...
        'EdgeColor', 'none', 'FaceColor', 'interp');

    colormap(parula);
    colorbar;
    xlabel('Kuhn segments N');
    ylabel('Prestretch \lambda = L/(N b)');
    zlabel('Count');
    title('Joint distribution of N and prestretch (\lambda)');
    set(gca, 'FontSize', 16, 'LineWidth', 2, 'ZScale', 'linear', 'View', [45 30]);

    % Optional grid and lighting for clarity
    grid on;
    box on;
    hold off;

end

% Bimodal histograms
if strcmp(options.dist_type,'bimodal')
    type1 = Bonds(:,5) == 1;

    % Kuhn Segments
    nbins = max(10, min(80, ceil(sqrt(numel(Nvec)))));
    figure; hold on
    histogram(Nvec(type1),nbins,'FaceColor',[0.2 0.2 0.2],'FaceAlpha',0.8,'LineWidth',0.0005)
    histogram(Nvec(~type1),nbins,'FaceColor',[1 0 0],'FaceAlpha',0.8,'LineWidth',0.0005)
    xlabel('Assigned N per bond'); ylabel('Count');
    title('N distribution');

    % bond length
    figure; hold on; 
    histogram(Bonds(type1,4),50,'FaceColor',[0.2 0.2 0.2],'FaceAlpha',1,'LineWidth',0.0005);
    histogram(Bonds(~type1,4),50,'FaceColor',[1 0 0],'FaceAlpha',1,'LineWidth',0.0005);
    axis tight
    %axis([0.3 1.3 0 350]); 
    xlabel('Bond length L'); ylabel('Count'); title('Length distribution (final)')
    set(gca,'FontSize',16,'LineWidth',2)

    % prestretch
    b = options.b;
    lamvec1 = (Bonds(type1,4))./(Nvec(type1)*b);
    lamvec2 = (Bonds(~type1,4))./(Nvec(~type1)*b);
    
    N1 = options.bimodal.N1;
    N2 = options.bimodal.N2;

    figure; hold on;
    hold on
    yl = ylim;  % get current y-axis limits
    plot([1/sqrt(N1) 1/sqrt(N1)], yl, '--', 'LineWidth', 2, 'Color', [0.2 0.2 0.2]);
    plot([1/sqrt(N2) 1/sqrt(N2)], yl, '--', 'LineWidth', 2, 'Color', [1 0 0]);


    dlam = 0.01;
    if sum(type1) > 1
        histogram(lamvec1,ceil((max(lamvec1) - min(lamvec1))/dlam),'FaceColor',[0.2 0.2 0.2],'FaceAlpha',0.8,'LineWidth',0.0005);
    end
    if sum(~type1) > 1
        histogram(lamvec2,ceil((max(lamvec2) - min(lamvec2))/dlam),'FaceColor',[1 0 0],'FaceAlpha',0.8,'LineWidth',0.0005);
    end
    axis tight
    xlabel('Prestretch \lambda = L/(N*b)'); ylabel('Count'); title('Prestretch distribution (final)')
    set(gca,'FontSize',16,'LineWidth',2)
    axis([0 1 0 1000])

 
end


end