% Define the transcendental equation as a function handle, f(x)
f = @(x) ((cos(x).*sinh(x) - sin(x).*cosh(x))./(cos(x).*cosh(x)+1) ...
         + 1./x) .* ...
         ((cos(x).*sinh(x) + sin(x).*cosh(x))./(cos(x).*cosh(x)+1) ...
         + (cos(x).*cosh(x)+1)./(cos(x).*sinh(x) - sin(x).*cosh(x))) + ...
         ((sin(x).*sinh(x))./(cos(x).*cosh(x)+1)).^2;

% Plot the function over a range to visualize zeros (roots)
figure;
fplot(f,[0.01,15],'LineWidth',2); % Avoid x=0 since 1/x diverges!
xlabel('\lambda l');
ylabel('f(\lambda l)');
title('Plot of transcendent equation f(\lambda l)');
ylim([-10 10]); % Adjust to better see zeros
grid on;
hold on;

% (Optional) Draw a zero line for reference
yline(0, 'k--');


roots = zeros(4,1);
guesses = [1, 2, 3, 4.05];   % Use values based on your own plot!
for k = 1:4
    roots(k) = fzero(f, guesses(k));
end
disp('First 4 roots:');
disp(roots);



% Physical and geometric parameters (example values, replace as needed)
E  = 2.1e11;      % Young's modulus (N/m^2)
I  = 8.5e-6;      % Area moment of inertia (m^4)
A  = 0.015;       % Cross-sectional area (m^2)
rho = 7850;       % Density (kg/m^3)
l   = 10;         % Length of beam (m)

% Roots (lambda_n*l), e.g. from your root finding step  % Example: replace with your calculated roots

% Calculate omega_n for each root
omega_n = ((roots).^2 ./ l^2) * sqrt(E*I/(A*rho));

% Display results
disp('Natural frequencies omega_n (rad/s):');
disp(omega_n);

% (Optional) Convert to Hz
freq_Hz = omega_n/(2*pi);
disp('Natural frequencies (Hz):');
disp(freq_Hz);
