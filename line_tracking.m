clear all
close all
clc

% Initialization - STRAIGHT LINE
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           SOSPROGRAM  (Theorem 3.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.1;
k = 1;
%%%%% SOSPROGRAM FOR EACH POINT
pvar e1 e2 e4 e5
e = [e1;e2;e4;e5];              % Z(et) = et, "Application to the robot" section
ubar = [1.8,0];                 % equilibrium (v_r, w_r) for the line, eq. (4.3)
Ur = ubar;
umax = [2.5;2.5];                % actuator saturation limits, eq. (4.4)
f = [Ur(1)*e5+e2*ubar(2); Ur(1)*e4-e1*ubar(2); Ur(2)*e5; 0];
%g = [-1 e2; 0 -e1; 0 -e5-1; 0 e4];
% B(et) rescaled by umax - normalization to +-1, Section 2.4.2
g = [-1*umax(1) e2*umax(2); 0*umax(1) -e1*umax(2); 0*umax(1) (-e5-1)*umax(2); 0*umax(1) e4*umax(2)];
% Vector of monomials chosen
Z = e;
% Pseudo-linear model of the saturated error, eq. (3.19)-(3.20)
A =  [0 ubar(2) 0 Ur(1);
    -ubar(2) 0 Ur(1) 0;
    0 0 0 Ur(2);
    0 0 0 0];
B = g;
% Dimensions
n = length(e);
m = size(g,2);
nz = length(Z);
% definition of X = E(Sx), analysis region, eq. (3.22)
Sx = eye(nz)*13;                % S_X = 13 I, Results section - straight-line trajectory
epsi = 1e-6;                    % epsilon slack in condition eq. (3.31)
% -- SOS PROGRAM (decision variables of Theorem 3.5)
prog = sosprogram(e);
[prog,Q] = sospolymatrixvar(prog,monomials(e,0),[nz,nz],'symmetric');   % Q, degree 0
[prog,K] = sospolymatrixvar(prog,monomials(e,0:1),[m,nz]);              % K(et), degree 0 to 1
[prog,T] = sospolymatrixvar(prog,monomials(e,0:1),[m,nz]);              % T(et), degree 0 to 1
gx = 1-Z'*Sx*Z;                 % l(et), boundary of X, eq. (3.22)
geq = e4^2+e5^2-1;              % algebraic constraint G(et)=0, eq. (3.4)
for i=1:nz
    for j=1:n
    M(i,j) = diff(Z(i),e(j));   % matrix M(et) = dZ/de, used in Lemma 3.2
    end
end
% -- Theorem 3.5, condition eq. (3.31) (S-procedure over Vdot, Lemma 3.2)
% vertOmega, vertices of Omega=[0,1]^2, Proposition 2.12
th1 = [0;0];
th2 = [0;1];
th3 = [1;0];
th4 = [1;1];
vert_omega = [th1 th2 th3 th4];
n_vert = 2^m; % number of vertices of Omega, |vertOmega|=2^m
for i=1:n_vert
    [prog,Sr] = sospolymatrixvar(prog,monomials(e,1:2),[nz,nz],'symmetric');   % S_R^i, SOS multiplier
    [prog,Seq] = sospolymatrixvar(prog,monomials(e,0:1),[nz,nz],'symmetric'); % S_eq, free multiplier for the equality
    TH = eye(m)*vert_omega(i);
    % F1^i(et,theta) from condition eq. (3.31)
    F1 = -(M*A*Q+M*B*(TH*K+(eye(m)-TH)*T) + Q*A'*M' + (K'*TH'+T'*(eye(m)-TH'))*B'*M') - Sr*gx - Seq*geq;
    S1 = F1 - epsi*eye(nz);     % epsilon*I slack for the strict inequality
    prog = sosineq(prog,S1);
end
% Theorem 3.5, condition eq. (3.32) (Schur complement, Lemma 3.3)
for j=1:m
    [prog,Sh] = sospolymatrixvar(prog,monomials(e,1:2),[nz,nz],'symmetric');  % S_H^j
    [prog,Seq] = sospolymatrixvar(prog,monomials(e,0:1),[nz,nz],'symmetric');
    F2 = Q-Sh*gx-Seq*geq;       % F2^j(et)
    S2 = [1 T(j); T(j)' F2(j)]; % Schur block of condition eq. (3.32)
    prog = sosineq(prog,S2);
end
% -- Theorem 3.5, condition eq. (3.33) (Lemma 3.4, E(Q^-1) subset X)
[prog,Seq] = sospolymatrixvar(prog,monomials(e,0:1),[nz,nz],'symmetric');
S3 = inv(Sx)-Q-Seq*geq;
prog = sosineq(prog,S3);
% -- Optimization criterion: maximize tr(Q), enlarges the ellipsoid E(Q^-1)
obj = trace(Q);
prog = sossetobj(prog,-obj);
% SDP solution, eq. (3.35), via MOSEK
sol = sossolve(prog);
Q = double(sosgetsol(sol,Q));
K = sosgetsol(sol,K);
T = sosgetsol(sol,T);
% Lyapunov function and gains, eq. (3.34)
V = Z'*inv(Q)*Z;                % V(et) = Z' Q^-1 Z
F = K*inv(Q);                   % F(et) = K(et) Q^-1
H = T*inv(Q);                   % H(et) = T(et) Q^-1
u = F*Z;                        % control law u(et) = F(et) Z(et)
v_u = H*Z;                      % auxiliary control nu(et) = H(et) Z(et)

toc

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           TRAJECTORY GENERATION 2 - closed-loop simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initializations for reference trajectory - Euler integration, eq. (4.1)
dt = 0.1;          % sampling time
Tsim = 20;         % simulation time
h = dt;            % integrating time
xr = [0;-1;0];
Xr = xr;           %Euler

% Initializations for real trajectory
xo = [-2;-1.5;0];
x = xo;
X = x;
Ubar = [];

syms E1 E2 E4 E5
E = [E1;E2;E4;E5-1];
Z_mon = monomials(E,1:4); % Note: renamed to avoid clashing with Z from pvar

% Error history, used in the phase portrait (Figure 13)
E1_hist = zeros(1, round(Tsim/h));
E2_hist = zeros(1, round(Tsim/h));
E3_hist = zeros(1, round(Tsim/h));

figure('Name', '2D Trajectory');
for k = 1:round(Tsim/h)  
    %%%%%%% 
    Xr = [Xr xr+h*([Ur(1)*cos(xr(3));Ur(1)*sin(xr(3));Ur(2)])];  % eq. (4.1)
    xr = Xr(:,k+1);
    plot(Xr(1,:),Xr(2,:),'b.'); drawnow
    hold on
    
    a = Xr(:,k)-x;
    % postural error in the robot's own frame, eq. (2.14)
    E1 = cos(x(3))*a(1)+sin(x(3))*a(2);
    E2 = -sin(x(3))*a(1)+cos(x(3))*a(2);
    E3 = a(3);
    
    E1_hist(k) = E1;
    E2_hist(k) = E2;
    E3_hist(k) = E3;
    
    E4 = sin(E3);
    E5 = cos(E3);
    
    uk = matlabFunction(p2s(u));
    U = uk(E1,E2,E4,E5-1); 
    
    % de-normalization and saturation Phi(u), Section 2.4.2 / eq. (2.3)
    v(k) = U(1)*umax(1)+ubar(1);
    v(k) = min(umax(1),max(-umax(1),v(k)));
    w(k) = U(2)*umax(2)+ubar(2);
    w(k) = min(umax(2),max(-umax(2),w(k)));
    
    X = [X x+dt*[v(k)*cos(x(3));v(k)*sin(x(3));w(k)]];
    x = X(:,k+1);
    plot(X(1,:),X(2,:),'+'); drawnow
end

% 2D tracking, Figure 10
plot(Xr(1,:),Xr(2,:),'','LineWidth',2);
hold on
plot(Xr(1,1),Xr(2,1),'o','MarkerEdgeColor','b','MarkerFaceColor','b');
hold on
plot(X(1,:),X(2,:),'+','MarkerEdgeColor','r','LineWidth',1);
hold on
plot(X(1,1),X(2,1),'o','MarkerEdgeColor','r','MarkerFaceColor','r');
legend('Virtual Robot','q_{r}(0)','Real Robot','q_{c}(0)')
xlabel('x (m)')
ylabel('y (m)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Plotting v and w - Figure 11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name', 'Control Signals');
time = dt*(0:k-1);
subplot(2,1,1)
plot(time,v,'b','LineWidth',1.2)
xlabel('Time (s)')
ylabel('v (m/s)')
title('Linear Velocity (v)')

hold on
subplot(2,1,2)
plot (time,w,'r','LineWidth',1.2)
xlabel('Time (s)')
ylabel('\omega (rad/s)')
title('Angular Velocity (\omega)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Plotting the x, y and theta error - Figure 12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Error = Xr-X;
Error4 = sin(Error(3,:));
Error5 = cos(Error(3,:));
figure('Name', 'Error Evolution');
plot(time,Error(1,1:end-1),time,Error(2,1:end-1),time,Error(3,1:end-1),time,Error4(1:end-1),time,Error5(1:end-1),'LineWidth',1.2);
legend('e_1','e_2','e_3','e_4','e_5')
xlabel('Time (s)')
ylabel('Error (m)')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Phase portrait: E(Q^-1), region X, boundaries of L(H) and trajectory
%      Figure 13 - Lemma 3.3, Lemma 3.4, Theorem 3.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Define symbolic variables
syms e1 e2 e3 e4 e5

% 2. Translated state vector Z_sym, with e4=sin(e3), e5=cos(e3)-1
Z_sym = [e1; e2; sin(e3); cos(e3)-1];

% 3. Invariant ellipsoid E(Q^-1) = {V(et) <= 1}, red
V_Lyap = Z_sym' * inv(Q) * Z_sym;

% 3.5 Analysis region X = E(Sx), eq. (3.22) (coincides here with E(Q^-1),
%     since S_X = Q^-1 = 13 I at the optimum, condition eq. (3.33) active)
X_ellip = Z_sym' * Sx * Z_sym;

% 4. Boundaries of L(H): nu(et) = H(et) Z(et), region of validity of the
%    polytopic saturation model, Lemma 3.3, |h_j(et)Z(et)| <= 1
v_u_sym = p2s(v_u); 
v_u_sym = subs(v_u_sym, [e4, e5], [sin(e3), cos(e3)-1]);
f1_sym = v_u_sym(1); % nu_1, v_x actuator
f2_sym = v_u_sym(2); % nu_2, w actuator

% 5. Figure
figure('Name', '3D Stability Analysis with Trajectory and Region X', 'Color', 'w');
hold on; grid on; view(3);
xlabel('e_1 (m)', 'FontWeight', 'bold');
ylabel('e_2 (m)', 'FontWeight', 'bold');
zlabel('e_3 (rad)', 'FontWeight', 'bold');

% Viewing limits (zoom in)
lim = 1.0; 
plot_limits = [-lim lim -lim lim -lim lim];
res = 50; 

% Region X (cyan)
fimplicit3(X_ellip - 1, plot_limits, 'FaceColor', 'c', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'MeshDensity', res);

% Invariant ellipsoid E(Q^-1) (red)
fimplicit3(V_Lyap - 1, plot_limits, 'FaceColor', 'r', 'FaceAlpha', 0.8, 'EdgeColor', 'none', 'MeshDensity', res);

% Boundary of L(H) associated with nu_1 (v_x), |nu_1| <= 1
fimplicit3(f1_sym - 1, plot_limits, 'FaceColor', 'y', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'MeshDensity', res);
fimplicit3(f1_sym + 1, plot_limits, 'FaceColor', 'y', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'MeshDensity', res);

% Boundary of L(H) associated with nu_2 (w), |nu_2| <= 1
fimplicit3(f2_sym - 1, plot_limits, 'FaceColor', 'g', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'MeshDensity', res);
fimplicit3(f2_sym + 1, plot_limits, 'FaceColor', 'g', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'MeshDensity', res);

% ---> ERROR TRAJECTORY <---
% Postural error evolution over time (blue line)
plot3(E1_hist, E2_hist, E3_hist, 'b-', 'LineWidth', 2.5);

% Initial error condition (blue dot)
plot3(E1_hist(1), E2_hist(1), E3_hist(1), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);

% Origin, equilibrium e=0 (black star)
plot3(0, 0, 0, 'k*', 'MarkerSize', 12, 'LineWidth', 2);

% Legend
camlight; lighting gouraud;
legend_elements = [patch(NaN, NaN, 'c', 'FaceAlpha', 0.15), ...
                   patch(NaN, NaN, 'r', 'FaceAlpha', 0.8), ...
                   patch(NaN, NaN, 'y', 'FaceAlpha', 0.15), ...
                   patch(NaN, NaN, 'g', 'FaceAlpha', 0.15), ...
                   plot(NaN, NaN, 'b-', 'LineWidth', 2.5), ...
                   plot(NaN, NaN, 'k*', 'MarkerSize', 10)];
legend(legend_elements, {'Region \chi', 'Ellipsoid V(e)=1', 'Boundary \nu_1', 'Boundary \nu_2', 'Error Trajectory e(t)', 'Origin (e=0)'}, 'Location', 'best');