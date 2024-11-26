function [Y] = sensor_wo_noise(X, At, rho)

Y = mass_to_height(X, At, rho);
