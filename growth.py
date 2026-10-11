# ref: (Vico & Porporato, 2013), (Goudriaan & Monteith, 1990)
import numpy as np


class LAIDevelopmentSSWM:
    
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
        0                                      if s <= s_tilde
        g+ * ((s - s_tilde)/(xi - s_tilde))^a if s_tilde < s < xi
        g+                                     if s >= xi
        Refer to Section 2.1.4 of my thesis (for this simulation, s_tilde = s* & xi = s_wilt)
        (Vico & Porporato, 2013)
        """
        s = np.asarray(s, dtype=float)
        denom = self.xi - self.s_tilde
        if denom <= 0:
            raise ValueError("xi must be greater than s_tilde")
        stress_frac = (s - self.s_tilde) / denom
        stress_frac = np.clip(stress_frac, 0.0, 1.0)
        return self.g_plus * stress_frac ** self.a

    def mean_g(self):
        """⟨g⟩ = Σ_i Ω(s_i) · p0(s_i)  (Vico & Porporato 2013)"""
        p0 = np.asarray(self.sm_model.p0, dtype=float)
        s  = np.linspace(0.0, 1.0, len(p0))
        return float(np.sum(self.g_of_s(s) * p0))

    def lai_increment(self):
        """
        Annual LAI increment (m2 m-2 yr-1):

            ΔLAI = <g> · T_GS · (1 - φ)

        Refer to Section 2.2 of my thesis 
        (Goudriaan & Monteith, 1990)
        """
        return self.mean_g() * self.Tg * (1.0 - self.sm_model.phiE)


