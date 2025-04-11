% Define constants (equivalent to Python's constants module)
AU = 1.496e11;
G = 6.67430e-11;
M = 1.989e30;

% Main script execution (equivalent to Python's if __name__ == "__main__":)
initial_speed = 15 * 1000;  % 15 km/s
initial_angle_rad = 2 * pi * (30 / 360);  % 30 degrees
initial_x = 0;  % m
initial_y = -1 * AU;  % m

% Simulate Python's InitialConditions object creation by creating the vector directly
y0 = [
    initial_x;
    initial_y;
    initial_speed * cos(initial_angle_rad);
    initial_speed * sin(initial_angle_rad)
];

solved = solve_ode_and_get_values(y0, G, M); % Pass G, M explicitly
plot_trajectory(solved, G, M, 'trajectory.png'); % Pass G, M for consistency if needed by plot_trajectory

% --- Function Definitions ---

% Note: In MATLAB, helper functions used by top-level script code
%       are typically defined below the main script or in separate files.
%       Local functions (like ode_func here) must be at the end of the file.


function solved_output = solve_ode_and_get_values(y0, G, M) % Changed from initial_conditions object to y0 vector
    % Using RK45 method equivalent (ode45) to solve the ODE
    % d2r/dt2 = -G * M / r^2, r = (x^2 + y^2) ** 0.5 (in 2d)
    % Initial conditions passed in y0
    t_points = linspace(0, 2e7, 200000);

    t_span = [t_points(1), t_points(end)]; % Define t_span for ode45

    % solve_ivp -> ode45 This function was swaped from the following original python import (from scipy.integrate import solve_ivp)
    [t, y_matrix] = ode45(@(t_arg, y_arg) ode_func(t_arg, y_arg, G, M), t_span, y0, odeset('OutputFcn', @odeprog, 'Refine',1), t_points); % Pass G, M to ode_func, evaluate at t_points

    % Replicate the output structure of solve_ivp roughly
    solved_output.t = t.'; % Transpose t to row vector like solve_ivp
    solved_output.y = y_matrix.'; % Transpose y_matrix so rows are variables, columns are time steps
    % solved_output.sol = {t, y_matrix}; % Store raw ode45 output if needed later? No, match solve_ivp fields.
    solved_output.status = 0; % Assume success, ode45 doesn't return status like solve_ivp easily
    % Can add more fields like 'message' if needed based on try/catch around ode45
end


function plot_trajectory(solved, G, M, filename) % Added G, M for consistency
    % Extract the solution
    y_values = solved.y; % y_values is already [vars; time]
    x = y_values(1, :);
    y = y_values(2, :);
    t = solved.t; % Get time points

    % Plotting the trajectory

    r = (x.^2 + y.^2).^0.5; % Use element-wise operators .^
    % Find indexs where the sign of the derivative changes
    dr_dt = gradient(r(:), t(:)); % Ensure column vectors for gradient time spacing
    dr_dt_sign = sign(dr_dt);
    sign_change_indices = find(diff(dr_dt_sign) ~= 0); % Indices k where sign(k+1) != sign(k)


    % Limit data to the second periapsis
    try
        limit_index = sign_change_indices(3); % Python's sign_change[2] corresponds to 3rd change (index 3)
        x = x(1:limit_index);
        y = y(1:limit_index);
        r = r(1:limit_index); % Also limit r
    catch ME % MATLAB equivalent of except
        fprintf('Could not find second periapsis, using all data\n'); % Use fprintf for print
        % Let x, y, r use all data if fewer than 3 changes found
    end

    % Find the periapsis
    [min_r, periapsis_index] = min(r); % min returns value and index
    periapsis_x = x(periapsis_index);
    periapsis_y = y(periapsis_index);

    figure; % Equivalent to creating a plot figure implicitly with plt.plot
    plot(x, y);
    hold on; % Need hold on before scatter
    xlabel('x position (m)');
    ylabel('y position (m)');
    scatter(periapsis_x, periapsis_y, 'red', 'filled', 'DisplayName', sprintf('Periapsis at r = %.2e m', min_r)); % Use filled scatter and DisplayName
    legend('show'); % Explicitly show legend
    title('Trajectory of the object');
    grid on;
    axis equal;
    hold off; % Release hold
    % Save to windows clipboard (Not directly possible - save to file)
    % Replicate exact path string
    fullFilePath = ['/mnt/c/Users/hegra/Documents/', filename];
    try
        saveas(gcf, fullFilePath); % gcf gets current figure handle
    catch saveME
        fprintf('Failed to save plot to %s: %s\n', fullFilePath, saveME.message);
    end
end


function dy = ode_func(t, y, G, M)
    % Function to define the system of ODEs.
    % y(1) = x position
    % y(2) = y position
    % y(3) = x velocity
    % y(4) = y velocity

    r = (y(1)^2 + y(2)^2)^0.5;
    % Avoid division by zero if r is extremely small
    if r < eps
        ax = 0;
        ay = 0;
    else
        inv_r3 = 1 / r^3; % Calculate 1/r^3 once
        ax = -G * M * y(1) * inv_r3;
        ay = -G * M * y(2) * inv_r3;
    end

    % Return derivatives as a column vector
    dy = [
        y(3);
        y(4);
        ax;
        ay
    ];
end
