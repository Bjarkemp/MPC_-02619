function xdot = QuadrupleTankProcess(t, x, u, p)
    % Parameters
    a1 = p(1); a2 = p(2); a3 = p(3); a4 = p(4);
    A1 = p(5); A2 = p(6); A3 = p(7); A4 = p(8);
    g = p(9); gamma1 = p(10); gamma2 = p(11); rho = p(12);

    % masses
    m1 = x(1); m2 = x(2); m3 = x(3); m4 = x(4);

    % Flows
    q1in = gamma1 * u(1);
    q2in = gamma2 * u(2);
    q3in = (1 - gamma2) * u(2);
    q4in = (1 - gamma1) * u(1);

    q1 = a1 * sqrt(2 * g * m1/(rho * A1));
    q2 = a2 * sqrt(2 * g * m2/(rho * A2));
    q3 = a3 * sqrt(2 * g * m3/(rho * A3));
    q4 = a4 * sqrt(2 * g * m4/(rho * A4));

    % Differential equations (mass)
    dm1dt = (q1in + q3 - q1) * rho;
    dm2dt = (q2in + q4 - q2) * rho;
    dm3dt = (q3in - q3) * rho;
    dm4dt = (q4in - q4) * rho;

    % Output xdot
    xdot = [dm1dt; dm2dt; dm3dt; dm4dt];
end
