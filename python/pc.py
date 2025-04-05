import math
import numpy as np

from objects.initial_conditions import InitialConditions
from pa import solve_ode_and_get_values, plot_trajectory, AU
from objects.initial_conditions import InitialConditionsType

def handle_solved(solved, target_time_days: int = 210) -> tuple[float, float]:
    """
    finds the periapsis of first orbit given a solved trajectory
    """

    # Extract the solution
    y_values = solved.y
    x = y_values[0]
    y = y_values[1]
    # Extract the time points
    t_points = solved.t

    # Interpolate the distance over time with splines
    from scipy.interpolate import CubicSpline
    x_interpolated = CubicSpline(t_points, x)
    y_interpolated = CubicSpline(t_points, y)

    # Find value at target time
    target_time = target_time_days * 24 * 3600  # Convert days to seconds

    target_x = x_interpolated(target_time)
    target_y = y_interpolated(target_time)
    return target_x, target_y

def solve_ode(initial_angle: float, initial_speed: float = 23 * 1000):
    initial_angle_rad = 2 * math.pi * (initial_angle / 360)  # 30 degrees
    initial_x = 0  # m
    initial_y = -1 * AU  # m

    initial_conditions = InitialConditions(
        x=initial_x,  # 1 AU
        y=initial_y,
        vx=initial_speed * math.cos(initial_angle_rad),  # 15 km/s
        vy=initial_speed * math.sin(initial_angle_rad),  # 29.78 km/s (orbital speed of Earth)
        type=InitialConditionsType.FIXED_VELOCITY
    )
    solved = solve_ode_and_get_values(initial_conditions)
    new_x, new_y = handle_solved(solved)
    return (new_x, new_y, solved)

def solver(target_x: float = 1, target_y: float = 1, initial_angle: float = 30, initial_speed: float = 26 * 1000, guess_diff: float = 0.1, max_attempts: int = 10):

    attempts = 0

    new_angle = initial_angle
    new_speed = initial_speed

    min_residual = 1e20

    for i in range(MAX_ITER):

        if attempts > 10:
            print("Too many attempts, giving up")
            return None

        # Get initial points with newton rapson

        x, y, _ = solve_ode(new_angle, new_speed)

        # Calculate the residual
        residual_1 = target_x - x
        residual_2 = target_y - y
        residual = math.sqrt(residual_1 ** 2 + residual_2 ** 2) # Check if the residual is within the tolerance

        if residual < 1:
            found = True
            break

        if residual < min_residual:
            min_residual = residual

        print(f"Minimum residual: {min_residual:.6f}", end=" ") print(f"Iteration {i}: Residual = {residual:.6f}", end=" ")
        print(f"Angle: {new_angle:.6f} degrees", end=" ") print(f"Speed: {new_speed:.6f} m/s", end=" ")
        print(f"Residual 1: {residual_1:.6f} m", end=" ")
        print(f"Residual 2: {residual_2:.6f} m")

        # Find the partials using finite difference

        h = 1e1
    
        x_dangle_1, y_dangle_1, _ = solve_ode(new_angle + h, new_speed)
        x_dangle_2, y_dangle_2, _ = solve_ode(new_angle - h, new_speed)

        x_dspeed_1, y_dspeed_1, _ = solve_ode(new_angle, new_speed + h)
        x_dspeed_2, y_dspeed_2, _ = solve_ode(new_angle, new_speed - h)

        # Use f' = (f(x+h) - f(x-h)) / (2h)
        dx_dangle = -(x_dangle_1 - x_dangle_2) / (2 * h)
        dy_dangle = -(y_dangle_1 - y_dangle_2) / (2 * h)
        dx_dspeed = -(x_dspeed_1 - x_dspeed_2) / (2 * h)
        dy_dspeed = -(y_dspeed_1 - y_dspeed_2) / (2 * h)

        # Calculate the Jacobian matrix
        J = np.array([[dx_dangle, dx_dspeed],
                      [dy_dangle, dy_dspeed]])


        # Find residual vector (angle, speed)*T
        residual_vector = np.array([residual_1, residual_2])

        try:
            # Solve for the delta using the Jacobian
            delta = np.linalg.solve(J, residual_vector)

            delta_angle = delta[0]
            delta_speed = delta[1]

            # Update with controlled step size
            new_angle = (new_angle - delta_angle) % 360
            new_speed = new_speed - delta_speed

        except np.linalg.LinAlgError:

            print("Jacobian is singular, trying again")
            initial_angle += np.random.uniform(-guess_diff, guess_diff)
            attempts += 1
            continue


    final_results = solve_ode(new_angle, new_speed)
    new_x, new_y, solved = final_results

    # Print results
    if found:
        print(f"Found solution in {i} iterations")
    else:
        print("Could not find solution within tolerance")

    print(f"Initial angle: {initial_angle} degrees")
    print(f"Final angle: {new_angle} degrees")
    print(f"Initial speed: {initial_speed} m/s")
    print(f"Final speed: {new_speed} m/s")

    print(f"Final Residual: {residual_2}")

    plot_trajectory(solved, filename="trajectory_p3.png")

MAX_ITER = 1000000

if __name__ == "__main__":

    target_x = 0
    target_y = 2.7e10

    guess_diff = 0.1

    initial_angle = 320
    initial_speed = 26000# 24.2 km/s

    solver(target_x, target_y, initial_angle, initial_speed, guess_diff)
