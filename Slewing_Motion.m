%% slewing_crane_14state.m
% 14-state slewing crane model with rope elasticity + boom bending + torsion
clear; close all; clc;

%% ----------------- PARAMETERS -----------------
g = 9.81;

% Mass / geometry
m_payload = 500;      % payload mass (kg)
L0 = 30;              % nominal rope length (m)
J_phi = 5e5;          % slewing inertia (kg*m^2)
m_boom = 1000;        % effective boom mass for bending mode (kg)

% Damping / stiffness (tune as needed)
d_phi_drive = 2e3;    % drive/shaft viscous damping (N*m*s/rad)
k_torsion = 1e6;      % torsional stiffness of shaft (N*m/rad)
d_torsion = 1e4;      % torsional damping (N*m*s/rad)

d_eta = 200;          % sway radial damping (N*s/m equivalent / scaled)
d_zeta = 200;         % sway tangential damping

k_l = 5e4;            % rope axial stiffness (N/m)
d_l = 1e3;            % rope axial damping (N*s/m)

% Boom bending (single-mode per plane approximated)
k_boom_y = 1e6;       % vertical bending stiffness (N/m)
d_boom_y = 5e3;       % vertical damping
k_boom_z = 1e6;       % horizontal bending stiffness (N/m)
d_boom_z = 5e3;       % horizontal damping
m_bend_eff = 500;     % effective modal mass for a single bending mode

% Other params
R = 20;               % radial offset of rope attach (m) - geometry
% note: many coupling coefficients are linearized approximations

%% ----------------- INPUT: Slew command -----------------
t_end = 5; dt = 0.001;
t = 0:dt:t_end;
t_data = t;

omega_max = 0.25; t_ramp = 2;
omega_cmd = zeros(size(t));
for i=1:length(t)
    if t(i) < t_ramp
        omega_cmd(i) = omega_max * (t(i)/t_ramp);
    elseif t(i) < (t_end - t_ramp)
        omega_cmd(i) = omega_max;
    else
        omega_cmd(i) = omega_max*(1 - (t(i)-(t_end-t_ramp))/t_ramp);
    end
end

%% ----------------- INITIAL CONDITIONS -----------------
phi0 = 0; phidot0 = 0;
eta0 = 0.05; etadot0 = 0;
zeta0 = 0.02; zetadot0 = 0;
l0 = L0; ldot0 = 0;
y0 = 0.01; ydot0 = 0;
z0 = 0.005; zdot0 = 0;
theta0 = 0.01; thetadot0 = 0;

x0 = [phi0; phidot0; eta0; etadot0; zeta0; zetadot0; l0; ldot0; y0; ydot0; z0; zdot0; theta0; thetadot0];

%% ----------------- SOLVE ODE -----------------
odefun = @(tt,x) crane14(tt, x, t_data, omega_cmd, ...
    m_payload, L0, J_phi, ...
    d_phi_drive, k_torsion, d_torsion, ...
    d_eta, d_zeta, ...
    k_l, d_l, ...
    k_boom_y, d_boom_y, k_boom_z, d_boom_z, m_bend_eff, g, R);

opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.02);
[t_sol, x_sol] = ode45(odefun, t, x0, opts);

% Interpolate onto uniform t (solver already used t, but keep consistent)
phi = interp1(t_sol, x_sol(:,1), t);
phidot = interp1(t_sol, x_sol(:,2), t);
eta = interp1(t_sol, x_sol(:,3), t);
etadot = interp1(t_sol, x_sol(:,4), t);
zeta = interp1(t_sol, x_sol(:,5), t);
zetadot = interp1(t_sol, x_sol(:,6), t);
l = interp1(t_sol, x_sol(:,7), t);
ldot = interp1(t_sol, x_sol(:,8), t);
y_b = interp1(t_sol, x_sol(:,9), t);
y_bdot = interp1(t_sol, x_sol(:,10), t);
z_b = interp1(t_sol, x_sol(:,11), t);
z_bdot = interp1(t_sol, x_sol(:,12), t);
theta_t = interp1(t_sol, x_sol(:,13), t);
theta_tdot = interp1(t_sol, x_sol(:,14), t);

%% ----------------- Derived quantities -----------------
% Bending moments (approx): M = k_boom * deflection
M_y = k_boom_y .* y_b;   % N-m approximated as stiffness * deflection
M_z = k_boom_z .* z_b;
M_torsion = k_torsion .* theta_t;

%% ----------------- PLOTS -----------------
figure('Position',[100 100 1400 900]);

subplot(4,3,1);
plot(t, omega_cmd,'k','LineWidth',1.5, Color='White'); grid on; title('Slew Command \omega (rad/s)');
xlabel('t (s)');

subplot(4,3,2);
plot(t, rad2deg(eta),'b','LineWidth',1.2); grid on; title('Radial Sway \eta (deg)');
xlabel('t (s)');

subplot(4,3,3);
plot(t, rad2deg(zeta),'r','LineWidth',1.2); grid on; title('Tangential Sway \zeta (deg)');
xlabel('t (s)');

subplot(4,3,4);
plot(t, rad2deg(phi),'g','LineWidth',1.2); grid on; title('\phi (deg)');
xlabel('t (s)');

%subplot(4,3,5);
%plot(t, l - L0,'m','LineWidth',1.2); grid on; title('Rope Elastic Extension (m)');
%xlabel('t (s)');

subplot(4,3,5);
plot(t, ldot,'c','LineWidth',1.2); grid on; title('Rope Extension Rate (m/s)');
xlabel('t (s)');

%subplot(4,3,7);
%plot(t, y_b,'LineWidth',1.2); grid on; title('Boom vertical deflection y_b (m)');
%xlabel('t (s)');

%subplot(4,3,8);
%plot(t, z_b,'LineWidth',1.2); grid on; title('Boom horiz deflection z_b (m)');
%xlabel('t (s)');

subplot(4,3,6);
plot(t, M_y/1e3,'LineWidth',1.2); grid on; title('Bending moment M_y (kN*m)');
xlabel('t (s)');

subplot(4,3,7);
plot(t, M_z/1e3,'LineWidth',1.2); grid on; title('Bending moment M_z (kN*m)');
xlabel('t (s)');

subplot(4,3,8);
plot(t, M_torsion/1e3,'LineWidth',1.2); grid on; title('Torsion moment (kN*m)');
xlabel('t (s)');

subplot(4,3,9);
plot(t, rad2deg(theta_t),'LineWidth',1.2); grid on; title('Torsion angle (deg)');
xlabel('t (s)');

sgtitle('14-state Slewing Crane Dynamics');

%% ----------------- SAVE RESULTS -----------------
Results = table(t', omega_cmd', phi', phidot', eta', etadot', zeta', zetadot', ...
    l', ldot', y_b', y_bdot', z_b', z_bdot', theta_t', theta_tdot', ...
    'VariableNames', {'Time_s','Omega_cmd','Phi','PhiDot','Eta','EtaDot','Zeta','ZetaDot', ...
    'RopeLength_m','RopeLenDot_mps','BoomY_m','BoomYdot_mps','BoomZ_m','BoomZdot_mps','Theta_t_rad','Theta_tdot'});
writetable(Results,'slewing_crane_results_14state.csv');
disp('Saved slewing_crane_results_14state.csv');

%% ----------------- END MAIN -----------------
%% ----------------- ODE function -----------------
function dx = crane14(t, x, t_data, omega_cmd, ...
        m_payload, L0, J_phi, ...
        d_phi_drive, k_torsion, d_torsion, ...
        d_eta, d_zeta, ...
        k_l, d_l, ...
        k_boom_y, d_boom_y, k_boom_z, d_boom_z, m_bend_eff, g, R)

% Unpack states
phi = x(1); phidot = x(2);
eta = x(3); etadot = x(4);
zeta = x(5); zetadot = x(6);
l = x(7); ldot = x(8);
y_b = x(9); y_bdot = x(10);
z_b = x(11); z_bdot = x(12);
theta_t = x(13); theta_tdot = x(14);

% Interpolate omega command
omega = interp1(t_data, omega_cmd, t, 'linear', 'extrap');

% ---------------- Slewing (torsion spring + damping drive) ----------------
% Torque from shaft/drive: -k_torsion*theta_t - d_torsion*theta_tdot - d_phi_drive*(phidot - omega)
% Here theta_t is torsion twist; assume platform rotation phi closely follows torsion twist
% Use torque balance for platform rotation (simplified)
T_drive = -k_torsion*(theta_t) - d_torsion*(theta_tdot) - d_phi_drive*(phidot - omega);
phidd = T_drive / J_phi;

% ---------------- Rope elasticity (axial) ----------------
% Simple mass-spring-damper for rope extension
% Stiffness force (pulling payload toward boom): F_spring = -k_l*(l - L0) - d_l*ldot
% Effective axial acceleration (approx linearized)
% Add a small geometric coupling from centripetal terms:
f_axial = -k_l*(l - L0) - d_l*ldot + m_payload*(phidot^2 * (R + l*eta)); % simple centrifugal approx
ldd = f_axial / m_payload;

% ---------------- Sway equations (coupled) ----------------
% Effective pendulum length is current l
% Note: linearized small-angle model with coupling to phidd and phidot
etadd = -(d_eta/m_payload)*etadot - (g / l)*eta ...
        - phidd * zeta - 2*phidot * zetadot ...
        + ( - (y_b / l) ); % couple small boom vertical deflection (geometric)

zetadd = -(d_zeta/m_payload)*zetadot - (g / l)*zeta ...
         + phidd * eta + 2*phidot * etadot ...
         + ( - (z_b / l) ); % couple small boom horizontal deflection

% ---------------- Boom bending modes ----------------
% Simple second-order modal equation per plane (mass-spring-damper)
% Coupled to sway through forcing proportional to payload lateral acceleration
% Forcing approximated as m_payload * lateral accel projection / modal mass
% Vertical boom mode (y)
F_boom_y = m_payload * ( etadd * l + 2*etadot*ldot ); % crude coupling
y_bdd = ( -k_boom_y * y_b - d_boom_y * y_bdot + F_boom_y ) / m_bend_eff;

% Horizontal boom mode (z)
F_boom_z = m_payload * ( zetadd * l + 2*zetadot*ldot );
z_bdd = ( -k_boom_z * z_b - d_boom_z * z_bdot + F_boom_z ) / m_bend_eff;

% ---------------- Torsion dynamics ----------------
% Torsion twist dynamics around shaft: torque balance
% Torsional damping & stiffness plus reaction from drive torque
% We model theta_t acceleration from net torque (drive torque transmitted minus internal)
% For numerical simplicity use:
theta_tdd = ( -k_torsion*theta_t - d_torsion*theta_tdot + d_phi_drive*(phidot - omega) ) / (J_phi*0.1);
% note: dividing by a fractional inertia to produce reasonable twist dynamics
% (tune to match actual shaft inertia if known)

% ---------------- Assemble dx ----------------
dx = zeros(14,1);
dx(1) = phidot;
dx(2) = phidd;
dx(3) = etadot;
dx(4) = etadd;
dx(5) = zetadot;
dx(6) = zetadd;
dx(7) = ldot;
dx(8) = ldd;
dx(9) = y_bdot;
dx(10) = y_bdd;
dx(11) = z_bdot;
dx(12) = z_bdd;
dx(13) = theta_tdot;
dx(14) = theta_tdd;

end
