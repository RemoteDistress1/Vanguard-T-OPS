class OChemAnalyst:
    def formal_charge(self, v, b, l):
        """FC = Valence - (Bonds + Lone Electrons)"""
        return v - (b + l)

    def get_dou(self, c, h, n, x):
        """Degree of Unsaturation: (2C + 2 + N - H - X) / 2"""
        return (2*c + 2 + n - h - x) / 2