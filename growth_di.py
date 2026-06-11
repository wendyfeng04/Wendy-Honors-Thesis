import numpy as np


class CropDevelopmentSSWM:
    "ref: Vico & Porporato (2013)"
    def __init__(self, sm_model, g_plus, a, xi, s_tilde, b0, Tg):
        self.sm_model = sm_model
        self.g_plus = g_plus
        self.a = a
        self.xi = xi
        self.s_tilde = s_tilde
        self.b0 = b0
        self.Tg = Tg

    def get_ps(self):
        ps_mass = np.asarray(self.sm_model.p0, dtype=float)
        s = np.linspace(0.0, 1.0, len(ps_mass))
        return s, ps_mass / (s[1] - s[0])

    def get_Ps(self):
        ps_mass = np.asarray(self.sm_model.p0, dtype=float)
        s = np.linspace(0.0, 1.0, len(ps_mass))
        return s, np.cumsum(ps_mass)

    def rho(self, s):
        return self.sm_model.loss(s)

    def ps_at_xi(self):
        s, ps = self.get_ps()
        return np.interp(self.xi, s, ps)

    def Ps_at_xi(self):
        s, Ps = self.get_Ps()
        return np.interp(self.xi, s, Ps)

    def g_of_s(self, s):
        s = np.asarray(s, dtype=float)
        out = np.where(s >= self.xi, self.g_plus, self.g_plus * (s / self.xi) ** self.a)
        return out

    def g_minus(self):
        return self.g_of_s(np.array([0.5 * (self.s_tilde + self.xi)]))[0]

    def T_up(self):
        return (1.0 - self.Ps_at_xi()) / (self.rho(self.xi) * self.ps_at_xi())
  
        

    def T_down(self):
        return self.Ps_at_xi() / (self.rho(self.xi) * self.ps_at_xi())
   

    def lambda_plus(self):
        return 1.0 / self.T_up()
    
    

    def lambda_minus(self):
        return 1.0 / self.T_down()
    
    def Ps_at_stilde(self):
        s, Ps = self.get_Ps()
        return np.interp(self.s_tilde, s, Ps)

    #def mean_g(self):
        #lp, lm, gm = self.lambda_plus(), self.lambda_minus(), self.g_minus()
        #return (lp * gm + lm * self.g_plus) / (lp + lm)
    def mean_g(self):
        Ps = self.Ps_at_xi()      # P(S <= xi)
        ps = self.ps_at_xi()      # p_s(xi)
        rho_xi = self.rho(self.xi)
         # Case 1: soil moisture is basically always greater than xi
         # P(S <= xi) = 0, so growth should be maximum
        if np.isfinite(Ps) and Ps <= 1e-8:
            return self.g_plus

         # Case 2: soil moisture is basically always below xi
         # P(S <= xi) = 1, so growth should be stressed growth
        if np.isfinite(Ps) and Ps >= 1.0 - 1e-8:
            return self.g_minus()
        #Case 3: always below s_tilde
        Ps_tilde = self.Ps_at_stilde()
        if np.isfinite(Ps_tilde) and Ps_tilde >= 1.0 - 1e-8:
            return 0.0
        else:
            lp = self.lambda_plus()
            lm = self.lambda_minus()
            gm = self.g_minus()
            return (lp * gm + lm * self.g_plus) / (lp + lm)
        
        


    def mu(self, t):
        return self.mean_g() * np.asarray(t, dtype=float) + self.b0

    def sigma(self, t):
        t = np.asarray(t, dtype=float)
        lp, lm = self.lambda_plus(), self.lambda_minus()
        dg = self.g_plus - self.g_minus()
        term = 2.0 * lp * lm * (np.exp(-(lp + lm) * t) + (lp + lm) * t - 1.0) / (lp + lm) ** 3
        return np.sqrt(term * dg ** 2)
    
    