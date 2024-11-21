function h = mass_to_height(x, p)

A=p(1:end-1);                   %[cm2] Cross sectional area in all tanks
rho=p(end);                  % [g/cm^3] Density of water

h= x./ (rho * A');          % converting mass in [g] to height in [cm]