import math
import numpy as np

from objects.initial_conditions import InitialConditions
from pa import solve_ode_and_get_values, plot_trajectory, AU
from objects.initial_conditions import InitialConditionsType

def handle_solved(solved) -> float:
    """
    finds the periapsis of first orbit given a solved trajectory
    """

    # Extract the solution
    y_values = solved.y
    x = y_values[0]
    y = y_values[1]
    # Extract the time points
    t_points = solved.t

    # Calculate the distance from the origin
    r = (x**2 + y**2)**0.5

    # Interpolate the distance over time with splines

    from scipy.interpolate import CubicSpline
    r_interpolated = CubicSpline(t_points, r)

    # Find derivate
    dr_dt = r_interpolated.derivative()

    # Find the first 0 crossing of the derivative
    zero_crossing = np.where(np.diff(np.sign(dr_dt(t_points))))[0]
    if len(zero_crossing) == 0:
        raise ValueError("Could not find periapsis time, try increasing the time span")

    periapsis_time = t_points[zero_crossing[0]]
    periapsis_value = float(r_interpolated(periapsis_time))  # Get the value at the zero crossing

    return periapsis_value

def solve_ode(new_angle: float):
    initial_speed = 23 * 1000  # 15 km/s
    initial_angle_rad = 2 * math.pi * (new_angle / 360)  # 30 degrees
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
    periapsis_value = handle_solved(solved)
    return (periapsis_value, solved)

def solver(target_periapsis: float = 1 * AU, initial_angle: float = 30, guess_diff: float = 0.1, attempts: int = 0):
    if attempts > 10:
        print("Too many attempts, giving up")
        return None

    # Get initial points for secant method
    p1, _ = solve_ode(initial_angle)
    residual_1 = p1 - target_periapsis

    initial_angle_2 = initial_angle + guess_diff
    p2, _ = solve_ode(initial_angle_2)
    residual_2 = p2 - target_periapsis

    # Calculate first new angle using secant method
    try:
        new_angle = initial_angle_2 - residual_2 * (initial_angle_2 - initial_angle) / (residual_2 - residual_1)
    except ZeroDivisionError:
        # If residuals are identical, try a different guess
        return solver(target_periapsis, np.random.uniform(0, 360), guess_diff, attempts + 1)

    # Iterative solver
    found = False
    for i in range(MAX_ITER):
        print(f"Iteration {i}, angle: {new_angle}, residual: {residual_2}")

        # Save current values before updating
        angle_prev = initial_angle_2
        residual_prev = residual_2

        # Solve with new angle
        p_new, solved = solve_ode(new_angle)
        residual_new = p_new - target_periapsis

        # Update for next iteration
        initial_angle = initial_angle_2
        initial_angle_2 = new_angle
        residual_1 = residual_2
        residual_2 = residual_new

        # Check for convergence
        if abs(residual_new) < 1e-4:
            found = True
            break

        # Calculate next angle using secant method
        try:
            new_angle = initial_angle_2 - residual_2 * (initial_angle_2 - initial_angle) / (residual_2 - residual_1)
        except ZeroDivisionError:
            # If secant method fails, try with a small perturbation
            new_angle = initial_angle_2 + np.random.uniform(-1, 1) * guess_diff

        # Ensure angle stays in reasonable range
        while new_angle < 0 or new_angle > 360:
            new_angle = new_angle % 360

    # Print results
    if found:
        print(f"Found solution in {i} iterations")
    else:
        print("Could not find solution within tolerance")

    print(f"Initial angle: {initial_angle} degrees")
    print(f"Final angle: {new_angle} degrees")
    print(f"Final Residual: {residual_2}")

    plot_trajectory(solved, filename="trajectory_p2.png")
    return new_angle, solved

MAX_ITER = 1000

if __name__ == "__main__":
    
    target_periapsis = 2.7e10
    initial_angle = 30
    guess_diff = 0.1
    solver(target_periapsis, initial_angle, guess_diff)
