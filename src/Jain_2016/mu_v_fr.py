import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Parameters
k_p = 250
k = 5e-4
m_t = 500
m_r = 10000
v_p = v_t = v_r = 1e-8
f_t = 0.50
f_r = 0.50
q = 4
V_c = 1.0 #um^3
d_t = 0.1 #as used in Jain_2016(h^-1)
d_r = 0.0 #as used in Jain_2016(h^-1)

# Initial Conditions
P0 = 1e6
T0 = 1e6
R0 = 1e6
y0 = [P0, T0, R0]

t = np.linspace(0,100,100)
f_r_range = np.linspace(0,1.0, 100)

def PTRD(y, t, k_p, q, k, m_t, m_r, v_p, v_t, v_r, f_r, f_t, d_t, d_r):
    P, T, R = y
    V = v_p * P + v_t * T + v_r * R
    K_p = k_p * q
    dPdt = K_p * T - k * R * P / V
    dTdt = (k / m_t) * (R * P / V) * f_t - d_t * T
    dRdt = (k / m_r) * (R * P / V) * f_r - d_r * R
    
    return [dPdt, dTdt, dRdt]


mu_values_q1 = []
for f_r in f_r_range:
    f_t = 1 - f_r
    solution = odeint(PTRD, y0, t, args=(k_p, q, k, m_t, m_r, v_p, v_t, v_r, f_r, f_t, d_t, d_r))
    P = solution[:, 0]
    coeff = np.polyfit(t, np.log(P), 1)
    mu = coeff[0]
    mu_values_q1.append(mu)


max_mu = max(mu_values_q1)
max_f_r = f_r_range[mu_values_q1.index(max_mu)]

print(f'Maximum mu: {max_mu}')
print(f'Corresponding f_r: {max_f_r}')

# Plot the results
plt.figure(figsize=(8,6))
plt.plot(f_r_range, mu_values_q1, label='q = 4')
plt.xlabel('f_r')
plt.ylabel('Mu')
plt.title("Mu vs. f_r for q = 4")
plt.legend()
plt.show()