function u_mpc = Uncon_MPC(MPC_sys, x, R)
    % Extract matrices for QP problem
    g = MPC_sys.M_x0 * x + MPC_sys.M_r * R;
    H = MPC_sys.H;

    % Solve unconstrained QP
    u_mpc = qpsolver(H, g, [], [], [], [], [], []);
end
