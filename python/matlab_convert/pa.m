% ============== FILE: pa.m ==============

% Constants (defined globally for script context)
AU = 1.496e11;  % m
G = 6.67430e-11;  % m^3 kg^-1 s^-2
M = 1.989e30;  % kg (mass of the sun)


% Main script execution (equivalent to if __name__ == "__main__":)
initial_speed = 15 * 1000;  % 15 km/s
initial_angle_rad = 2 * pi * (30 / 360);  % 30 degrees
initial_x = 0;  % m
initial_y = -1 * AU;  % m
% Simulating InitialConditions object creation by defining the initial state vector y0
% y0 = [x, y, vx, vy] - Column vector for ode45
y0 = [
    initial_x;
    initial_y;
    initial_speed * cos(initial_angle_rad);
    initial_speed * sin(initial_angle_rad)
];

[t_solved, y_solved_matrix] = solve_ode_and_get_values(y0, G, M); % Pass G, M explicitly
plot_trajectory(t_solved, y_solved_matrix, G, M, 'trajectory.png'); % Pass G, M explicitly


% --- Function Definitions ---

% Note: In MATLAB, functions usually reside in separate files or as local functions
% at the end of a script file. Placing them before the main script logic works too.

function dy = ode_func(t, y, G, M) % Pass G, M explicitly
    % Function to define the system of ODEs.
    % y(1) = x position
    % y(2) = y position
    % y(3) = x velocity
    % y(4) = y velocity
    % No docstring equivalent needed per instructions

    r = (y(1)^2 + y(2)^2)^0.5;
		inv_r3 = 1 / (r^3); % Calculate 1/r^3 once
		ax = -G * M * y(1) * inv_r3;
		ay = -G * M * y(2) * inv_r3;

    dy = [y(3); y(4); ax; ay]; % Return column vector
end


function [t_out, y_values_matrix] = solve_ode_and_get_values(y0, G, M) % Pass G, M explicitly
    % Using RK45 method to solve the ODE
    % d2r/dt2 = -G * M / r^2, r = (x^2 + y^2) ** 0.5 (in 2d)
    % Initial conditions are passed in y0
    t_points = linspace(0, 2e7, 200000);
    t_span = [0, 2e7]; % Define t_span explicitly for ode45

    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9); % Set tolerances similar to default RK45

    % --- MATLAB function swap ---
    [t_out, y_values_matrix] = ode45(@(t, y) ode_func(t, y, G, M), t_points, y0, options); % solve_ivp -> ode45 % This function was swaped from the following original python import (from scipy.integrate import solve_ivp)
    % --- End function swap ---

    % ode45 returns t and y matrix directly, no 'solved' object like solve_ivp
    % t_out corresponds to solved.t
    % y_values_matrix corresponds to solved.y.', where each *row* is a time step
end


function plot_trajectory(t, y_values_matrix, G, M, filename) % Take t and y_matrix, G, M as input
    % Extract the solution
    % y_values = solved.y -> y_values_matrix.' (transpose)
    x = y_values_matrix(:, 1); % First column is x
    y = y_values_matrix(:, 2); % Second column is y
    % Plotting the trajectory

    r = (x.^2 + y.^2).^0.5; % Use element-wise power .^
    % Find indexs where the sign of the derivative changes
    dr_dt = gradient(r); % gradient is MATLAB equivalent of np.gradient
    dr_dt_sign = sign(dr_dt); % sign is MATLAB equivalent of np.sign
    sign_change_indices = find(diff(dr_dt_sign)); % find(diff(...)) similar to np.where(np.diff(...))[0]

    % Limit data to the second periapsis
    try
        limit_index = sign_change_indices(3); % MATLAB uses 1-based indexing
        x = x(1:limit_index);
        y = y(1:limit_index);
        r = r(1:limit_index); % Also limit r for finding min later
    catch ME % Use catch ME for error handling
        fprintf('Could not find second periapsis, using all data\n'); % Use fprintf for print
    end

    % Find the periapsis
    [min_r, periapsis_index] = min(r); % min returns value and index
    periapsis_x = x(periapsis_index);
    periapsis_y = y(periapsis_index);

    % Plotting commands - MATLAB equivalents
    figure; % Create a new figure window explicitly
    plot(x, y);
    hold on; % Keep plot active for scatter
    xlabel('x position (m)');
    ylabel('y position (m)');
    scatter(periapsis_x, periapsis_y, 'red', 'filled', 'DisplayName', sprintf('Periapsis at r = %.2e m', min_r)); % Use scatter, 'filled', DisplayName for legend label
    legend('show'); % Show legend
    title('Trajectory of the object');
    grid on;
    axis equal;
    hold off; % Release plot hold
    % Save to file
    try
        saveas(gcf, filename); % Use saveas, gcf gets current figure handle
    catch saveME
        fprintf('Error saving file to %s: %s\n', filename, saveME.message);
        % Attempt saving to current directory as fallback
        [~, name, ext] = fileparts(filename);
        fallback_filename = [name, ext];
        try
            saveas(gcf, fallback_filename);
            fprintf('Saved plot to current directory as: %s\n', fallback_filename);
        catch finalSaveME
            fprintf('Failed to save plot even in current directory: %s\n', finalSaveME.message);
        end
    end

end
