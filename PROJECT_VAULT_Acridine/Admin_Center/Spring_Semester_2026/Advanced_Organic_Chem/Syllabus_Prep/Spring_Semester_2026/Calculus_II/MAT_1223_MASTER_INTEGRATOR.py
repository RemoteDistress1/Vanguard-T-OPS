import numpy as np
from scipy.integrate import quad

class CalculusMaster:
    def verify_volume(self, f, a, b):
        """Disk Method: pi * integral of f(x)^2"""
        val, _ = quad(lambda x: np.pi * (f(x)**2), a, b)
        return val

    def verify_volume_shell(self, r_func, h_func, a, b):
        """Shell Method: 2*pi * integral of r(x)*h(x)"""
        integrand = lambda x: 2 * np.pi * r_func(x) * h_func(x)
        val, _ = quad(integrand, a, b)
        return val