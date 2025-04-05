from scipy.integrate import solve_ivp
from math import sin, cos
import math
import numpy as np
from objects.initial_conditions import InitialConditions
from constants import AU, G, M

def ode_func(t, y):
    """
    Function to define the system of ODEs.
    y[0] = x position
    y[1] = y position
    y[2] = x velocity
    y[3] = y velocity
    """

    r = (y[0]**2 + y[1]**2)**0.5
    ax = -G * M * y[0] / r**3
    ay = -G * M * y[1] / r**3

    return [y[2], y[3], ax, ay]


def solve_ode_and_get_values(initial_conditions: InitialConditions):
    # Using RK45 method to solve the ODE
    # d2r/dt2 = -G * M / r^2, r = (x^2 + y^2) ** 0.5 (in 2d)
    # Initial conditions
    t_points = np.linspace(0, 2e7, 200000)

    solved = solve_ivp(
        fun=ode_func,
        t_span=(0, 2e7),
        y0=[
            initial_conditions.x,
            initial_conditions.y,
            initial_conditions.vx,
            initial_conditions.vy
        ],
        method='RK45',
        t_eval=t_points,
        vectorized=False,
    )
    return solved

def plot_trajectory(solved, filename: str = "trajectory.png"):
    # Extract the solution
    y_values = solved.y
    x = y_values[0]
    y = y_values[1]
    # Plotting the trajectory

    r = (x**2 + y**2)**0.5
    # Find indexs where the sign of the derivative changes
    dr_dt = np.gradient(r)
    dr_dt_sign = np.sign(dr_dt)
    sign_change = np.where(np.diff(dr_dt_sign))[0]

    # Limit data to the second periapsis
    try:
        x = x[:sign_change[2]]
        y = y[:sign_change[2]]
    except:
        print("Could not find second periapsis, using all data")

    # Find the periapsis
    periapsis_index = np.argmin(r)
    periapsis_x = x[periapsis_index]
    periapsis_y = y[periapsis_index]

    import matplotlib.pyplot as plt
    plt.plot(x, y)
    plt.xlabel('x position (m)')
    plt.ylabel('y position (m)')
    plt.scatter(periapsis_x, periapsis_y, color='red', label=f'Periapsis at r = {r[periapsis_index]:.2e} m')
    plt.legend()
    plt.title('Trajectory of the object')
    plt.grid()
    plt.axis('equal')
    # Save to windows clipboard
    plt.savefig(f'/mnt/c/Users/hegra/Documents/{filename}')


if __name__ == "__main__":
    initial_speed = 15 * 1000  # 15 km/s
    initial_angle_rad = 2 * math.pi * (30 / 360)  # 30 degrees
    initial_x = 0  # m
    initial_y = -1 * AU  # m
    initial_conditions = InitialConditions(
        initial_x,
        initial_y,
        initial_speed * cos(initial_angle_rad),
        initial_speed * sin(initial_angle_rad)
    )

    solved = solve_ode_and_get_values(initial_conditions)
    plot_trajectory(solved)
