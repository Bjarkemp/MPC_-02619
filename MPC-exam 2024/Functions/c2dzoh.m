% function [Abar,Bbar]=c2dzoh(A,B,Ts)
%     [nx,nu]=size(B);
%     M = [A B; zeros(nu,nx) zeros(nu,nu)];
%     Phi = expm(M*Ts);
%     Abar = Phi(1:nx,1:nx);
%     Bbar = Phi(1:nx,nx+1:nx+nu);
% end

function [Abar,Bbar,Gbar]=c2dzoh(A,B,G,Ts)
    % A: Systemmatrix
    % B: Inputmatrix
    % G: Støjinputmatrix
    % Ts: Samplingstiden
    [nx,nu]=size(B);  % Dimensioner af B
    [nx,ng]=size(G);  % Dimensioner af G

    % Udvidet matrix til at inkludere G
    M = [A B G; 
         zeros(nu,nx+nu+ng); 
         zeros(ng,nx+nu+ng)];
    
    % Matrixeksponential
    Phi = expm(M*Ts);

    % Ekstrahere de diskretiserede matricer
    Abar = Phi(1:nx,1:nx);         % Diskretiseret A
    Bbar = Phi(1:nx,nx+1:nx+nu);   % Diskretiseret B
    Gbar = Phi(1:nx,nx+nu+1:end);  % Diskretiseret G
end