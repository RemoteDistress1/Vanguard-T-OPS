import math

class PhysicsCommand:
    def __init__(self):
        self.g = 9.80665

    def solve_incline(self, mass, theta_deg, mu_k):
        """Verified Dynamics for Inclined Planes"""
        theta = math.radians(theta_deg)
        f_n = mass * self.g * math.cos(theta)
        f_p = mass * self.g * math.sin(theta)
        accel = (f_p - (mu_k * f_n)) / mass
        return {'Normal_Force': f_n, 'Acceleration': accel}