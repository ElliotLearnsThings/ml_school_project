import enum

class InitialConditionsType(enum.Enum):
    """Enum for initial conditions type."""
    FULLY_SPECIFIED = 1
    FIXED_VELOCITY = 2
    FIXED_ANGLE = 3

class InitialConditions:
    def __init__(self,
                 x: float,
                 y: float,
                 vx: float,
                 vy: float,
                 type: InitialConditionsType = InitialConditionsType.FULLY_SPECIFIED
                 ):
        self.x = x  # m
        self.y = y  # m
        self.vx = vx  # m/s
        self.vy = vy  # m/s
        self.type = type

    def __repr__(self):
        return f"InitialConditions(x={self.x}, y={self.y}, vx={self.vx}, vy={self.vy})"

    def __str__(self):
        return f"Initial Conditions: x={self.x}, y={self.y}, vx={self.vx}, vy={self.vy}"
