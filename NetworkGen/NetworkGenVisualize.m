function NetworkGenVisualize(Domain,Atoms,Bonds,Nvec,options)
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

% --------------------- Plot final bonds ---------------------
figure; hold on;

scatter(Atoms(1:N_atom,2), Atoms(1:N_atom,3), 6, 'k', 'filled');
for k = 1:Total_bond
    if Bonds(k,1) == 0, continue; end
    i1 = Bonds(k,2); i2 = Bonds(k,3);
    if Bonds(k,5) == 1
        plot([Atoms(i1,2) Atoms(i2,2)], [Atoms(i1,3) Atoms(i2,3)], 'k-');
    else
        plot([Atoms(i1,2) Atoms(i2,2)], [Atoms(i1,3) Atoms(i2,3)], 'r-');
    end
    
end
axis equal tight; title('Final bonds (post-prune)');

% --------------------- Network stats (Kuhn length) ---------------------

% ---------- Optional histogram ----------
nbins = max(10, min(80, ceil(sqrt(numel(Nvec)))));
figure; hist(Nvec, nbins);
xlabel('Assigned N per bond'); ylabel('Count');
title('N distribution');

% Bimodal histograms
if strcmp(options.dist_type,'bimodal')
    type1 = Bonds(:,5) == 1;
    figure; hold on; 
    histogram(Bonds(type1,4),50,'FaceColor',[0.2 0.2 0.2],'FaceAlpha',1,'LineWidth',0.0005);
    histogram(Bonds(~type1,4),50,'FaceColor',[1 0 0],'FaceAlpha',1,'LineWidth',0.0005);
    axis tight
    %axis([0.3 1.3 0 350]); 
    xlabel('Bond length L'); ylabel('Count'); title('Length distribution (final)')
    set(gca,'FontSize',16,'LineWidth',2)
end


end