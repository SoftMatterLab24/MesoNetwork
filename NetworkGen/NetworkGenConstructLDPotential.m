function [LDpot] = NetworkGenConstructLDPotential(Domain,Atoms,Bonds,Nvec,options)
%NetworkGenConstructLDPotential - Construct local density potential for LAMMPS
 % 
 % Files:
 %   NetworkGenConstructLDPotential.m
 %   mytab.localdensity.table

% ---------- Unpack domain & defaults ----------
Atom_count = size(Atoms,1);
Bond_count = size(Bonds,1);
Total_kuhn_segment = sum(Nvec);

xlo = Domain.xlo; xhi = Domain.xhi;
ylo = Domain.ylo; yhi = Domain.yhi;

b = options.b;              % Kuhn length

kLD = options.LDpot_strength; % strength factor
N_rho = options.LDpot_N_rho; % number of density points
rho_min = options.LDpot_rho_min; % minimum density
rho_max = options.LDpot_rho_max; % maximum density
drho = (rho_max - rho_min) / (N_rho - 1); % density step


% ---------- Construct local density potential ----------

sig_c = 0.5*b*(Total_kuhn_segment/Atom_count)^(1/2);   % desired equlibrium length
atom_density = Atom_count / ((xhi - xlo)*(yhi - ylo)); % number density

R2 = 4*sig_c;            % outer radius of density calculation
rho0 = 0.8*(R2/sig_c)^2; % desired equilibrium density for given R2

R1 = 0.8*sig_c;      % inner radius of density calculation
rc = 1.0*2*sig_c;    % cutoff radius for bpm/spring repulsion

% Construct potential density function
rho_vec = linspace(rho_min, rho_max, N_rho+1)'; % density vector

pot_density = kLD * (rho_vec - rho0).^2; % harmonic potential around rho0


% ----------- Pack LDpot struct ----------
LDpot.N_LD          = 1;        % single local density potential
LDpot.N_rho         = N_rho;    % number of density points
LDpot.R_lower       = R1;       % 10% overlap with other repulisve potentials 
LDpot.R_upper       = R2;       % outer radius of density calculation
LDpot.rc            = rc;       % cutoff radius for bpm/spring repulsion
LDpot.rho_min       = rho_min; 
LDpot.rho0          = rho0;
LDpot.rho_max       = rho_max;
LDpot.drho          = drho;
LDpot.pot_density   = pot_density;

% ----------- Print summary -----------
fprintf('   Constructed LD potential with parameters:\n');
fprintf('   Target equilibrium density rho0 = %.4f\n', rho0);
fprintf('   Lower Cutoff R1 = %.4f * b, Upper Cutoff R2 = %.4f * b', R1/b, R2/b);
fprintf('   Cutoff for BPM/Spring repulsion rc = %.4f * b \n', rc/b);

end
