function [Y] = sensor_wo_noise(X, At, rho)
Y = mass_to_height(X(:,1:2)', At(1:2), rho);
