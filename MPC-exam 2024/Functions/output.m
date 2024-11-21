function [z] = output(X, At, rho)

z = mass_to_height(X(:,1:2), [At(1:2); rho]);