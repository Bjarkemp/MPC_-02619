function u_mpc = Uncon_MPC(MPC_sys, x, R, u_prev)
    % Extract matrices for QP problem
    x
    MPC_sys.M_x0
    g = MPC_sys.M_x0 * x + MPC_sys.M_r * R + MPC_sys.M_u1 * u_prev
    H = MPC_sys.H;

    % Solve unconstrained QP
    u_mpc = qpsolver(H, g, [], [], [], [], [], []);
end
