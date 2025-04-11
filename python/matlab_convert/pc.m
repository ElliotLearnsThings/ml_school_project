% MATLAB equivalent of the third Python script (Newton-Raphson Solver for 2 variables)

% Ensure the functions from the first conversion are available:
% - solve_ode_and_get_values(y0, G, M, t_span_seconds) -> [t, y_values_matrix]
%   (Note: Added t_span_seconds argument to ensure simulation runs long enough)
% - plot_trajectory(t, y_values_matrix, G, M, filename)
% - ode_func(t, y, G, M) -> dy (as a local function or separate .m file)
% Make sure G, M, AU constants are defined.

% Define constants (can be defined here or loaded/passed)
AU = 1.496e11; % Astronomical Unit in meters
G = 6.67430e-11; % Gravitational constant in m^3 kg^-1 s^-2
M = 1.989e30; % Mass of the central body (e.g., Sun) in kg
MAX_ITER = 1000; % Max iterations for the main solver loop (Python had 1,000,000?)
                 % Let's use a more reasonable number like 1000 first.

% --- Main script execution starts here ---
target_x = 0;         % Target x position (m)
target_y = 2.7e10;    % Target y position (m)
target_time_days = 210; % Target time in days (matches default in Python handle_solved)

guess_diff = 0.1;     % Perturbation size for Jacobian recovery (degrees/m/s)
max_attempts = 10;    % Max recovery attempts for singular Jacobian

initial_angle = 320;  % Initial guess for launch angle in degrees
initial_speed = 26000;% Initial guess for launch speed in m/s (26 km/s)

fprintf('Starting solver...\n');
fprintf('Target X: %.3e m, Target Y: %.3e m at %d days\n', target_x, target_y, target_time_days);
fprintf('Initial Guess Angle: %.4f deg, Speed: %.2f m/s\n', initial_angle, initial_speed);

% Call the solver function
[final_angle, final_speed, t_final, y_final] = solver( ...
    target_x, target_y, target_time_days, ...
    initial_angle, initial_speed, ...
    guess_diff, max_attempts, ...
    G, M, AU, MAX_ITER ...
);

% --- Optional: Display final results if solver succeeded ---
if ~isempty(final_angle)
    fprintf('\n--- Solver Finished ---\n');
    % Verify final position
    [final_calc_x, final_calc_y] = handle_solved(t_final, y_final, target_time_days);
    final_residual_vec = [target_x - final_calc_x; target_y - final_calc_y];
    final_residual_mag = norm(final_residual_vec);

    fprintf('Target Position:      (%.4e, %.4e) m\n', target_x, target_y);
    fprintf('Calculated Position:  (%.4e, %.4e) m\n', final_calc_x, final_calc_y);
    fprintf('Final Residual Mag:   %.4e m\n', final_residual_mag);
    fprintf('Final Angle:          %.6f degrees\n', final_angle);
    fprintf('Final Speed:          %.4f m/s\n', final_speed);

    % Plot the final trajectory
    % Assumes plot_trajectory is available
    plot_trajectory(t_final, y_final, G, M, 'trajectory_p3.png');
else
    fprintf('\n--- Solver Failed ---\n');
end


% --- Function Definitions ---

function [target_x, target_y] = handle_solved(t, y_values_matrix, target_time_days)
    % Finds the (x, y) position at a specific target time using interpolation.

    % Extract the solution components
    x = y_values_matrix(:, 1);
    y = y_values_matrix(:, 2);

    % Calculate target time in seconds
    target_time_sec = target_time_days * 24 * 3600;

    % Check if target time is within the simulation time span
    min_t = min(t);
    max_t = max(t);
    if target_time_sec < min_t || target_time_sec > max_t
        warning('Target time (%.2f s) is outside the simulation time span [%.2f s, %.2f s]. Result is extrapolated.', target_time_sec, min_t, max_t);
        % Consider throwing an error or adjusting simulation span if extrapolation is not desired
        % error('Target time is outside the simulation time span. Increase t_span in solve_ode_and_get_values.');
    end

    % Interpolate using cubic splines (requires Curve Fitting Toolbox for spline())
    % If toolbox is unavailable, consider 'pchip' or 'interp1' with 'spline' method
    try
        ppx = spline(t, x); % Piecewise polynomial structure for x(t)
        ppy = spline(t, y); % Piecewise polynomial structure for y(t)

        % Evaluate the splines at the target time
        target_x = ppval(ppx, target_time_sec);
        target_y = ppval(ppy, target_time_sec);
    catch ME
         if strcmp(ME.identifier, 'MATLAB:UndefinedFunction') && contains(ME.message, 'spline')
             fprintf('Warning: `spline` function (Curve Fitting Toolbox) not found. Using `interp1`.\n');
             target_x = interp1(t, x, target_time_sec, 'spline', 'extrap'); % Use interp1 as fallback
             target_y = interp1(t, y, target_time_sec, 'spline', 'extrap');
         else
             rethrow(ME); % Rethrow other errors
         end
    end
end


function [calc_x, calc_y, t, y_values_matrix] = solve_ode(initial_angle, initial_speed, G, M, AU, target_time_days)
    % Solves the ODE for given initial conditions and returns the calculated
    % position (x, y) at target_time_days, along with the full solution.

    initial_angle_rad = 2 * pi * (initial_angle / 360); % Convert angle to radians
    initial_x = 0;  % m
    initial_y = -1 * AU;  % m

    % Initial conditions vector [x; y; vx; vy] (column vector)
    y0 = [
        initial_x;
        initial_y;
        initial_speed * cos(initial_angle_rad);
        initial_speed * sin(initial_angle_rad)
    ];

    % Ensure the simulation runs long enough
    % Add a buffer to the required time
    t_end_seconds = target_time_days * 24 * 3600 * 1.1; % Simulate 10% longer

    % Call the external/previous ODE solver function
    % Assumes solve_ode_and_get_values is modified to accept t_end_seconds
    % or that its default span is sufficiently long.
    % If it takes t_span = [t_start, t_end]:
    [t, y_values_matrix] = solve_ode_and_get_values(y0, G, M, [0, t_end_seconds]);
    % If it takes only t_end:
    % [t, y_values_matrix] = solve_ode_and_get_values(y0, G, M, t_end_seconds);

    % Get the position at the target time using interpolation
    [calc_x, calc_y] = handle_solved(t, y_values_matrix, target_time_days);
end


function [final_angle, final_speed, t_final, y_final] = solver( ...
    target_x, target_y, target_time_days, ...
    initial_angle, initial_speed, ...
    guess_diff, max_attempts, ...
    G, M, AU, MAX_ITER ...
)
    % Uses a Newton-like method with finite differences for the Jacobian
    % to find the (initial_angle, initial_speed) that results in reaching
    % (target_x, target_y) at target_time_days.

    new_angle = initial_angle;
    new_speed = initial_speed;
    attempts = 0;       % Counter for recovery attempts
    min_residual = Inf; % Track minimum residual achieved
    found = false;      % Convergence flag
    t_final = [];       % To store the last valid trajectory time vector
    y_final = [];       % To store the last valid trajectory solution matrix

    fprintf('Starting Newton-Raphson iterations...\n');

    for i = 1:MAX_ITER
        if attempts > max_attempts
            fprintf('Maximum number of recovery attempts (%d) reached. Stopping.\n', max_attempts);
            final_angle = []; final_speed = []; % Indicate failure
            return;
        end

        % --- Step 1: Evaluate function F(angle, speed) ---
        try
            [x, y, t_curr, y_curr] = solve_ode(new_angle, new_speed, G, M, AU, target_time_days);
            t_final = t_curr; % Store latest valid solution data
            y_final = y_curr;
        catch ME
            fprintf('Error during ODE solve in iteration %d: %s\n', i, ME.message);
            fprintf('Attempting recovery by perturbing current guess.\n');
            % Perturb current guess slightly and try again
            new_angle = mod(new_angle + (rand*2 - 1) * guess_diff, 360);
            new_speed = new_speed + (rand*2 - 1) * 10; % Perturb speed by +/- 10 m/s
            attempts = attempts + 1;
            continue; % Skip to next iteration
        end

        % --- Step 2: Calculate Residual Vector R = Target - F ---
        residual_1 = target_x - x;
        residual_2 = target_y - y;
        residual_vector = [residual_1; residual_2]; % Column vector
        residual = norm(residual_vector); % Magnitude of residual

        % --- Print Progress ---
        if residual < min_residual
            min_residual = residual;
        end
        fprintf('Iter %d: MinRes=%.3e, CurRes=%.3e, Angle=%.5f, Speed=%.2f\n', ...
                 i, min_residual, residual, new_angle, new_speed);
        % Debugging: print residuals
        % fprintf('  ResX: %.3e, ResY: %.3e\n', residual_1, residual_2);


        % --- Step 3: Check Convergence ---
        tolerance = 1e4; % Convergence tolerance in meters (adjust as needed)
                         % Python used 'residual < 1', which is extremely tight.
                         % Let's use 10km = 1e4 m first.
        if residual < tolerance
            fprintf('Converged in %d iterations with residual %.4e m.\n', i, residual);
            found = true;
            break;
        end

        % --- Step 4: Calculate Jacobian J = dF/d(vars) using Finite Differences ---
        % Note: F = Target - Current, so dF/d(vars) = - d(Current)/d(vars)
        % Python code calculates -d(Current)/d(vars) directly.

        h_angle = 1e-3; % Step size for angle (degrees) - small is better! Python used 10.
        h_speed = 1.0;  % Step size for speed (m/s) - Python used 10.

        try
            % Perturb angle
            [x_ap, y_ap, ~, ~] = solve_ode(new_angle + h_angle, new_speed, G, M, AU, target_time_days);
            [x_am, y_am, ~, ~] = solve_ode(new_angle - h_angle, new_speed, G, M, AU, target_time_days);
            % Perturb speed
            [x_sp, y_sp, ~, ~] = solve_ode(new_angle, new_speed + h_speed, G, M, AU, target_time_days);
            [x_sm, y_sm, ~, ~] = solve_ode(new_angle, new_speed - h_speed, G, M, AU, target_time_days);

             % Central difference derivatives for d(Current)/d(vars)
            current_dx_dangle = (x_ap - x_am) / (2 * h_angle);
            current_dy_dangle = (y_ap - y_am) / (2 * h_angle);
            current_dx_dspeed = (x_sp - x_sm) / (2 * h_speed);
            current_dy_dspeed = (y_sp - y_sm) / (2 * h_speed);

            % Jacobian J = dF/d(vars) = - d(Current)/d(vars)
            J = -[current_dx_dangle, current_dx_dspeed;
                  current_dy_dangle, current_dy_dspeed];

            % Debugging: Check Jacobian values
            % disp(J);

            % Check if Jacobian calculation resulted in NaNs (e.g., if ODE failed)
            if any(isnan(J(:)))
                error('Jacobian contains NaN values. ODE solve likely failed during finite difference.');
            end

        catch ME % Catch errors during the 4 extra ODE solves
             fprintf('Error during Jacobian calculation in iteration %d: %s\n', i, ME.message);
             fprintf('Attempting recovery by perturbing current guess.\n');
             new_angle = mod(new_angle + (rand*2 - 1) * guess_diff, 360);
             new_speed = new_speed + (rand*2 - 1) * 10; % Perturb speed by +/- 10 m/s
             attempts = attempts + 1;
             continue; % Skip to next iteration
        end


        % --- Step 5: Solve Linear System J * delta = R ---
        try
            % Check condition number before solving (optional but good practice)
            if cond(J) > 1e12 % Check if matrix is ill-conditioned
                 fprintf('Warning: Jacobian is ill-conditioned (cond=%.2e). Results may be inaccurate.\n', cond(J));
            end

            % Use backslash operator to solve J * delta = residual_vector
            delta = J \ residual_vector;

            % Check if solution contains NaNs
             if any(isnan(delta))
                 error('Solver returned NaN delta. Jacobian might be effectively singular.');
             end

        catch ME
            % Handle cases where J is singular or nearly singular
             if strcmp(ME.identifier, 'MATLAB:singularMatrix') || strcmp(ME.identifier, 'MATLAB:nearlySingularMatrix') || contains(ME.message, 'NaN')
                 fprintf('Jacobian is singular or solver failed in iteration %d. Perturbing guess.\n', i);
                 % Perturb current guess slightly as a recovery strategy
                 new_angle = mod(new_angle + (rand*2 - 1) * guess_diff, 360);
                 new_speed = new_speed + (rand*2 - 1) * 10; % Perturb speed by +/- 10 m/s
                 attempts = attempts + 1;
                 continue; % Skip to next iteration
             else
                 rethrow(ME); % Rethrow unexpected errors
             end
        end

        % --- Step 6: Update Variables: vars_new = vars_old - delta ---
        % Newton's method: x_k+1 = x_k - J(x_k)^-1 * F(x_k)
        % We solved J * delta = F, so delta = J^-1 * F
        % Therefore, x_k+1 = x_k - delta
        delta_angle = delta(1);
        delta_speed = delta(2);

        % Apply update (maybe with damping factor alpha if needed, alpha=1 for standard Newton)
        alpha = 1.0;
        new_angle = new_angle - alpha * delta_angle;
        new_speed = new_speed - alpha * delta_speed;

        % Keep angle in [0, 360) range
        new_angle = mod(new_angle, 360);

        % Optional: Add constraint on speed if necessary (e.g., prevent negative speed)
        % if new_speed < 1000 % Example lower bound
        %     new_speed = 1000;
        % end

    end % End of iteration loop

    % --- Post-Loop Actions ---
    if ~found
        fprintf('Solver did not converge within %d iterations.\n', MAX_ITER);
        fprintf('Returning the best result found (residual: %.4e m).\n', min_residual);
        % Keep the last calculated angle/speed/trajectory
        final_angle = new_angle;
        final_speed = new_speed;
        % t_final and y_final hold the last successful trajectory data
        if isempty(t_final) % If solver failed on first iteration
             final_angle = []; final_speed = []; % Indicate complete failure
        end
    else
        % If found, the loop broke after the successful check.
        % new_angle and new_speed hold the converged values.
        % t_final and y_final hold the trajectory data for the converged solution.
        final_angle = new_angle;
        final_speed = new_speed;
    end

end


% --- Placeholder/Reminder for required functions ---
% Ensure these functions from the previous steps are available AND UPDATED.
% Specifically, solve_ode_and_get_values needs to handle the time span.

% --- Example Updated solve_ode_and_get_values ---
function [t, y_values_matrix] = solve_ode_and_get_values(y0, G, M, t_span_seconds)
    % Solves the ODE over the specified time span [t_start, t_end].
    % y0: initial conditions [x; y; vx; vy]
    % G, M: constants
    % t_span_seconds: [t_start, t_end] for the simulation in seconds

    if nargin < 4 || isempty(t_span_seconds)
        t_span_seconds = [0, 2e7]; % Default span if not provided
        warning('Using default time span for ODE solver.');
    end

    % Define time points for output (optional, ode45 adapts step size)
    % Using just the span allows ode45 to choose optimal steps.
    % t_points = linspace(t_span_seconds(1), t_span_seconds(2), 200000); % Optional fixed points

    options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9); % Adjust tolerances if needed

    % Solve the ODE using ode45
    [t, y_values_matrix] = ode45(@(t_arg, y_arg) ode_func(t_arg, y_arg, G, M), t_span_seconds, y0, options);
    % If using specific t_points:
    % [t, y_values_matrix] = ode45(@(t_arg, y_arg) ode_func(t_arg, y_arg, G, M), t_points, y0, options);
end

% --- plot_trajectory (Assumed unchanged from first conversion) ---
function plot_trajectory(t, y_values_matrix, G, M, filename)
    % Plots the trajectory and saves it.
    if isempty(t) || isempty(y_values_matrix)
         fprintf('No trajectory data to plot.\n');
         return;
    end
    x = y_values_matrix(:, 1);
    y = y_values_matrix(:, 2);
    figure; % New figure
    plot(x, y);
    hold on;
    % Add start/end points, origin, etc. if desired
    plot(x(1), y(1), 'go', 'MarkerFaceColor','green', 'DisplayName','Start');
    plot(x(end), y(end), 'rs', 'MarkerFaceColor','red', 'DisplayName','End');
    plot(0,0, 'bo', 'MarkerSize', 8, 'DisplayName','Origin (Star)');
    xlabel('x position (m)');
    ylabel('y position (m)');
    title(sprintf('Trajectory (End Time: %.1f days)', t(end)/(24*3600)));
    grid on;
    axis equal;
    legend('show', 'Location','best');
    hold off;
    try
        saveas(gcf, filename);
        fprintf('Plot saved to %s\n', filename);
    catch saveME
         fprintf('Could not save plot: %s\n', saveME.message);
    end
end

% --- ode_func (Assumed unchanged from first conversion) ---
function dy = ode_func(t, y, G, M)
    % Defines the system of ODEs for the two-body problem.
    r_sq = y(1)^2 + y(2)^2;
    r = sqrt(r_sq);
    if r < eps % Avoid division by zero / instability near origin
        ax = 0;
        ay = 0;
    else
        inv_r3 = 1 / (r_sq * r); % More efficient than r^3
        ax = -G * M * y(1) * inv_r3;
        ay = -G * M * y(2) * inv_r3;
    end
    dy = [y(3); y(4); ax; ay]; % Return column vector [vx; vy; ax; ay]
end
