import numpy as np

class CropDevelopmentSSWM:
    "ref: Vico & Porporato (2013)"
    def __init__(self, sm_model, g_plus, a, xi, s_tilde, b0, Tg):
        self.sm_model = sm_model
        self.g_plus   = g_plus
        self.a        = a
        self.xi       = xi
        self.s_tilde  = s_tilde
        self.b0       = b0
        self.Tg       = Tg

    def g_of_s(self, s):
        """
        Ω(s) from Vico & Porporato (2013) Eq. (2):
          0               if s ≤ s_tilde
          g+ * (s/ξ)^a   if s_tilde < s < ξ
          g+              if s ≥ ξ
        """
        s = np.asarray(s, dtype=float)
        return np.where(
            s <= self.s_tilde, 0.0,
            np.where(
                s >= self.xi, self.g_plus,
                self.g_plus * (s / self.xi) ** self.a
            )
        )

    def mean_g(self):
        """⟨g⟩ = Σ_i Ω(s_i) · p0(s_i)"""
        p0  = np.asarray(self.sm_model.p0, dtype=float)
        s   = np.linspace(0.0, 1.0, len(p0))
        return float(np.sum(self.g_of_s(s) * p0))

    def mu(self, t):
        """Mean biomass trajectory: μ(t) = ⟨g⟩·t + b0"""
        return self.mean_g() * np.asarray(t, dtype=float) + self.b0