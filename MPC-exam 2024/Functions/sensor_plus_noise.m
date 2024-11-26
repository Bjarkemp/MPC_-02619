function [Y] = sensor_plus_noise(X, At, rho, R)

Lr = chol(R,'lower');                                   % Cholesky-dekomposition. It just gives me the standard deviation instead of variance.
v = (Lr*randn(width(X),length(X)))';                    % Measurement noise. Follows normal distribution with mean=0 and has st.dev of Lr

Y = mass_to_height(X', At, rho)+v'; % Heights in all tanks