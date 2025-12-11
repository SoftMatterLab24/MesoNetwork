n = network;

n.alpha = 2;
n.stdN2 = 11.072*log(n.alpha) + 2.1987;

[Domain,Atoms,Bonds,Nvec,order] = generateNetwork(n);