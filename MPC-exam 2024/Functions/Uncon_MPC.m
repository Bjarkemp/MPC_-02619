function u_mpc = Uncon_MPC(MPC_sys, x, R)
    % Extract matrices for QP problem
    g = MPC_sys.M_x0 * x + MPC_sys.M_r * R + MPC_sys.M_u1 * u_prev;
    Hs = MPC_sys.Hs;
    H = MPC_sys.H;

    % Solve unconstrained QP
    u_mpc = qpsolver(H+Hs, g, [], [], [], [], [], []);
end
