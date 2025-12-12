%% Pendulum with Variable Cable Length, Stiffness, Damping, and Tension Calculation
% Parameters
g = 9.81;      % Gravity (m/s^2)
m = 1000;      % Payload mass (kg)
k = 1e5;       % Cable stiffness (N/m)
c = 10;        % Damping coefficient (N·s/m)
r0 = 30;       % Unstretched cable length (m)
theta0 = 5;    % Initial angle in degrees
t0 = 0;        % Start time (s)
tf = 30;       % End time (s)
dt = 0.02;     % Time step (s)

% Initial conditions (r, theta, r_dot, theta_dot)
IC = [r0; deg2rad(theta0); 0; 0];

% ODE function for variable-length pendulum with spring/damp
function dy = pendulum_ode(~, y)
    % y = [r; theta; r_dot; theta_dot]
    g = 9.81;
    m = 1000;
    k = 1e5;
    c = 10;
    r0 = 30;
    r = y(1);
    theta = y(2);
    r_dot = y(3);
    theta_dot = y(4);

    dy = zeros(4,1);
    dy(1) = r_dot;
    dy(2) = theta_dot;
    dy(3) = r*theta_dot^2 ...
          - (k/m)*(r - r0) ...
          - (c/m)*r_dot ...
          - g*cos(theta);
    dy(4) = -2*r_dot*theta_dot/r ...
          - (c/m)*theta_dot ...
          - (g/r)*sin(theta);
end

% Solve using ode45
tspan = t0:dt:tf;
[t, y] = ode45(@pendulum_ode, tspan, IC);

% Extract variables
r = y(:,1);
theta = y(:,2);
r_dot = y(:,3);
theta_dot = y(:,4);

% Cable tension
T = m*(r.*theta_dot.^2 - g*cos(theta)) + k.*(r - r0);

% Horizontal and vertical tension components
Rx = T .* sin(theta);
Ry = T .* cos(theta);

% Payload velocity
V = sqrt( (r_dot).^2 + (r.*theta_dot).^2 );

% Plot cable tension versus time and angle, like in thesis
figure;
subplot(2,1,1);
plot(t, T, 'b', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Cable Tension (N)');
grid on;
title('Cable Tension vs. Time');

subplot(2,1,2);
plot(t, theta*180/pi, 'r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Swing Angle (deg)');
grid on;
title('Swing Angle vs. Time');

% Optional: plot horizontal and vertical components
figure;
plot(t, Rx, 'g', t, Ry, 'm');
xlabel('Time (s)');
ylabel('Tension Components (N)');
legend('Horizontal Rx','Vertical Ry');
grid on;
title('Tension Components vs. Time');

%% ----- ADDITIONAL REQUIRED PLOTS -----

% 1) Trolley Position vs Time (fixed crane = 0)
x_trolley = zeros(size(t));  % change later if motion profile exists

figure;
plot(t, x_trolley, 'b', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Trolley Position (m)');
title('Trolley Position vs Time');
grid on;

% 2) Payload Velocity vs Time
figure;
plot(t, V, 'm', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Payload Velocity (m/s)');
title('Payload Velocity vs Time');
grid on;

% 3) Angular Velocity vs Swing Angle (Phase Diagram)
figure;
plot(theta * 180/pi, theta_dot, 'r', 'LineWidth', 2);
xlabel('Swing Angle (deg)'); ylabel('Angular Velocity (rad/s)');
title('Phase Plot: Angular Velocity vs Swing Angle');
grid on;
