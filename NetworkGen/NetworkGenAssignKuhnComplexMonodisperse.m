function Nvec = NetworkGenAssignKuhnComplexMonodisperse(Bonds, options)
% Monodisperse N per bond type (1/2/3). Type 4 reserved.
% Bonds layout: [bondID id1 id2 L bondType]

nb = size(Bonds,1);
Nvec = zeros(nb,1);
if nb == 0, return; end

btype = Bonds(:,5);

N1 = options.complex.N_type1;
N2 = options.complex.N_type2;
N3 = options.complex.N_type3;

Nvec(btype==1) = N1;
Nvec(btype==2) = N2;
Nvec(btype==3) = N3;

% Safety for any weird types
maskOther = ~(btype==1 | btype==2 | btype==3);
if any(maskOther)
    Nvec(maskOther) = N1;
end

fprintf('   Complex monodisperse N: type1=%g, type2=%g, type3=%g\n', N1, N2, N3);
end