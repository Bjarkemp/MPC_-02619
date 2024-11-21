function z = FourTankSystemOutput(mass, p)
    % Extract parameters
    At = p(5:8); % Tank cross-sectional areas [cm2]
    rho = p(12); % Density of water [g/cm3]

    % Calculate heights from masses (mass = height * At * rho)
    z = mass ./ (At * rho);
end