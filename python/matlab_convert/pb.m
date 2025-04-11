% Define constants if not already defined/imported from pa.m
% Assuming pa.m defines these globally or they are loaded
if ~exist('AU', 'var'); AU = 1.496e11; end
if ~exist('G', 'var'); G = 6.67430e-11; end
if ~exist('M', 'var'); M = 1.989e30; end

% Define MAX_ITER (equivalent to Python global)
MAX_ITER = 1000;

% Main script execution (equivalent to Python's if __name__ == "__main__":)
target_periapsis = 2.7e10;
initial_angle = 30;
guess_diff = 0.1;

% Call solver function, passing constants explicitly
solver(target_periapsis, initial_angle, guess_diff, G, M, AU, MAX_ITER); % Pass MAX_ITER


% --- Function Definitions ---

function periapsis_value = handle_solved(solved)
    % finds the periapsis of first orbit given a solved trajectory

    % Extract the solution
    y_values = solved.y;
    x = y_values(1, :);
    y = y_values(2, :);
    % Extract the time points
    t_points = solved.t;

    % Calculate the distance from the origin
    r = (x.^2 + y.^2).^0.5;

    % Interpolate the distance over time with splines
    % from scipy.interpolate import CubicSpline -> spline/ppval
    pp = spline(t_points, r); % Create spline structure

    % Find derivate
    % r_interpolated.derivative() -> fnder (Curve Fitting Toolbox) or gradient
    % Using fnder requires toolbox, gradient is built-in.
    % Let's use fnder to match Python's use of spline derivative, requires toolbox.
    % If toolbox unavailable, gradient(r, t_points) could approximate but isn't identical.
    try
        pp_deriv = fnder(pp); % Get derivative spline structure
        dr_dt_vals = ppval(pp_deriv, t_points); % Evaluate derivative at t_points
    catch ME
        if contains(ME.identifier, 'UndefinedFunction') && contains(ME.message, 'fnder')
             fprintf('Warning: fnder (Curve Fitting Toolbox) not found. Using gradient for derivative.\n');
             dr_dt_vals = gradient(r(:), t_points(:)); % Use gradient as fallback
         else
             rethrow(ME);
         end
    end


    % Find the first 0 crossing of the derivative
    sign_diff = diff(sign(dr_dt_vals));
    zero_crossing_indices = find(sign_diff ~= 0);

    if isempty(zero_crossing_indices)
        error('Could not find periapsis time, try increasing the time span'); % Use error for exceptions
    end

    first_zero_crossing_idx_before = zero_crossing_indices(1);
    periapsis_time = t_points(first_zero_crossing_idx_before); % Exact Python logic: time *before* crossing
    % Get the value at the zero crossing (evaluating original spline at that time)
    periapsis_value = ppval(pp, periapsis_time); % Evaluate original spline pp
end


function [periapsis_value, solved] = solve_ode(new_angle, G, M, AU) % Pass constants
    initial_speed = 23 * 1000;  % 15 km/s -> Corrected to 23 km/s as per Python
    initial_angle_rad = 2 * pi * (new_angle / 360);  % 30 degrees
    initial_x = 0;  % m
    initial_y = -1 * AU;  % m

    % Simulate InitialConditions object creation
    y0 = [
        initial_x; % 1 AU -> No, initial_x is 0
        initial_y;
        initial_speed * cos(initial_angle_rad); % 15 km/s -> No, uses calculated speed
        initial_speed * sin(initial_angle_rad)  % 29.78 km/s (orbital speed of Earth) -> No, uses calculated speed
        % type=InitialConditionsType.FIXED_VELOCITY -> Not needed in MATLAB version
    ];

    % Assume solve_ode_and_get_values is available (from pa.m or path)
    solved = solve_ode_and_get_values(y0, G, M); % Pass G, M
    periapsis_value = handle_solved(solved);
    % Return tuple (handled by MATLAB returning multiple output arguments)
end


function [final_angle, final_solved] = solver(target_periapsis, initial_angle, guess_diff, G, M, AU, MAX_ITER, attempts) % Pass constants and MAX_ITER

    % Handle optional argument 'attempts'
    if nargin < 8
        attempts = 0;
    end

    if attempts > 10
        fprintf('Too many attempts, giving up\n'); % Use fprintf for print
        final_angle = NaN; % Return NaN or empty to indicate failure
        final_solved = [];
        return;
    end

    % Get initial points for secant method
    [p1, ~] = solve_ode(initial_angle, G, M, AU);
    residual_1 = p1 - target_periapsis;

    initial_angle_2 = initial_angle + guess_diff;
    [p2, solved] = solve_ode(initial_angle_2, G, M, AU); % Keep track of last 'solved'
    residual_2 = p2 - target_periapsis;

    % Calculate first new angle using secant method
    delta_residual = residual_2 - residual_1;
    if abs(delta_residual) < eps % Check for near-zero division (eps is machine epsilon)
        fprintf('Initial residuals identical or too close. Restarting with random angle.\n');
        % If residuals are identical, try a different guess recursively
        % Use rand for np.random.uniform(0, 360) -> 360*rand
        [final_angle, final_solved] = solver(target_periapsis, 360*rand, guess_diff, G, M, AU, MAX_ITER, attempts + 1);
        return; % Exit current call
    else
        new_angle = initial_angle_2 - residual_2 * (initial_angle_2 - initial_angle) / delta_residual;
    end

    % Iterative solver
    found = false;
    final_solved = solved; % Initialize final_solved with the one from p2 calculation
    i = 0; % Initialize loop counter correctly

    for i = 1:MAX_ITER % MATLAB loop syntax
        fprintf('Iteration %d, angle: %.6f, residual: %.6e\n', i, new_angle, residual_2); % Use fprintf

        % Save current values before updating
        angle_prev = initial_angle_2; % Not used in Python logic, but good practice? No, remove if python didnt use. Python *did* use it implicitly via initial_angle = initial_angle_2 below.
        residual_prev = residual_2; % Not used in Python logic. Remove.

        % Store values for next secant step (as Python did implicitly)
        angle_k_minus_1 = initial_angle;   % Keep track of x_{k-1}
        angle_k = initial_angle_2; % Keep track of x_{k}
        residual_k_minus_1 = residual_1;   % Keep track of f(x_{k-1})
        residual_k = residual_2;         % Keep track of f(x_{k})

        % Solve with new angle
        [p_new, solved] = solve_ode(new_angle, G, M, AU);
        residual_new = p_new - target_periapsis;

        % Update for next iteration (Secant method state)
        initial_angle = angle_k;            % x_{k-1} becomes x_k
        initial_angle_2 = new_angle;        % x_k becomes x_{k+1} (the latest calculated angle)
        residual_1 = residual_k;            % f(x_{k-1}) becomes f(x_k)
        residual_2 = residual_new;          % f(x_k) becomes f(x_{k+1}) (the latest calculated residual)
        final_solved = solved; % Store the solution corresponding to the latest angle

        % Check for convergence
        if abs(residual_new) < 1e-4
            found = true;
            break; % Exit loop
        end

        % Calculate next angle using secant method
        delta_residual = residual_2 - residual_1;
        if abs(delta_residual) < eps
             fprintf('Secant method division by zero or near-zero. Perturbing angle.\n');
             % If secant method fails, try with a small perturbation
             % Use (rand*2-1) for np.random.uniform(-1, 1)
             new_angle = initial_angle_2 + (rand*2-1) * guess_diff;
        else
             new_angle = initial_angle_2 - residual_2 * (initial_angle_2 - initial_angle) / delta_residual;
        end


        % Ensure angle stays in reasonable range
        % Use mod() for angle wrapping, handles negatives correctly unlike rem() sometimes.
        new_angle = mod(new_angle, 360);
        % Python's while loop is equivalent to a single mod operation.

    end % End of loop

    % Print results
    if found
        fprintf('Found solution in %d iterations\n', i); % Use iteration counter i
    else
        fprintf('Could not find solution within tolerance\n');
    end

    fprintf('Initial angle: %.6f degrees\n', initial_angle); % This is angle_k (previous)
    fprintf('Final angle: %.6f degrees\n', initial_angle_2); % This is angle_{k+1} (the final one)
    fprintf('Final Residual: %.6e\n', residual_2); % The residual corresponding to initial_angle_2

    % Plot trajectory requires plot_trajectory function (assumed available from pa.m or path)
    plot_trajectory(final_solved, G, M, 'trajectory_p2.png'); % Pass G, M

    % Return the final angle and the corresponding solution object
    final_angle = initial_angle_2;
    % final_solved is already updated in the loop
end

% Local function(s) like ode_func, solve_ode_and_get_values, plot_trajectory
% should ideally be defined here if not in separate files or pa.m.
% Assuming they are available from pa.m run previously or on path.
