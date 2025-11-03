function [LDpot] = NetworkGenConstructLDPotential(Atoms,Bonds,Nvec,options)
%NetworkGenConstructLDPotential - Construct local density potential for LAMMPS
 % 
 % Files:
 %   NetworkGenConstructLDPotential.m
 %   mytab.localdensity.table

% ---------- Unpack domain & defaults ----------
Atom_count = size(Atoms,1);
Bond_count = size(Bonds,1);
Total_kuhn_segment = sum(Nvec);

k = options.LDpot_strength; % strength factor

% ---------- Construct local density potential ----------

sig_c = b*(Total_kuhn_segment/Atom_count)^(1/2); % 


% k (rho-rho0)^2






% ----------- Pack LDpot struct ----------
LDpot.N_LD          = 1; % single local density potential
LDpot.N_rho         = 1000; % number of density points
LDpot.R_lower       = 0; % 10% overlap with other repulisve potentials 
LDpot.R_upper       = 0;
LDpot.rho_min       = 0.0;
LDpot.rho_max       = 100;
LDpot.drho          = 0.01;
LDpot.pot_density   = pot_density;

end
