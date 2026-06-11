import numpy as np


class CropDevelopmentSSWM_log:
    "ref: Vico & Porporato (2013)"
    def __init__(self, sm_model, g_plus, a, xi, s_tilde, b0, Tg,
                 a_logistic=0.0, LAI_ref=None):
        self.sm_model   = sm_model
        self.g_plus     = g_plus
        self.a          = a
        self.xi         = xi
        self.s_tilde    = s_tilde
        self.b0         = b0
        self.Tg         = Tg
        self.a_logistic = a_logistic
        self.LAI_ref    = LAI_ref

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

    def delta_LAI(self, LAI):
        """
        ΔLAI = ⟨g⟩ · Tg · (1 - LAI/LAI_ref)^a_logistic

        a_logistic = 0 → original Vico: ΔLAI = ⟨g⟩ · Tg
        a_logistic > 0 → concave asymptotic recovery
        """
        mg = self.mean_g()
        if self.LAI_ref is None or self.a_logistic == 0.0:
            return mg * self.Tg
        LAI = min(LAI, self.LAI_ref)
        gap = (1.0 - LAI / self.LAI_ref) ** self.a_logistic
        return self.Tg * mg * gap

    def mu(self, t):
        return self.mean_g() * np.asarray(t, dtype=float) + self.b0