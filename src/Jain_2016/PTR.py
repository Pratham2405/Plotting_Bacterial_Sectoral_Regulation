import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#Parameters
k_p = 250
k = 5e-4
m_t = 500
m_r = 10000
v_p = v_t = v_r = 1e-8
f_t = 0.50
f_r = 0.50
q = 2
V_c = 1.0 #um^3
d_t = 0.1 #as used in Jain_2016(h^-1)
d_r = 0.0 #as used in Jain_2016(h^-1)

#Initial Conditions
P0 = 1e6
T0 = 1e6
R0 = 1e6
y0 = [P0, T0, R0]

t = np.linspace(0,100,100)



def PTR(y, t, k_p, q, k, m_t, m_r, v_p, v_t, v_r, f_r, f_t, d_t, d_r):
    P, T, R = y
    V = v_p * P + v_t * T + v_r * R
    K_p = k_p * q
    dPdt = K_p * T - k * R * P / V
    dTdt = (k / m_t) * (R * P / V) * f_t - d_t * T
    dRdt = (k / m_r) * (R * P / V) * f_r - d_r * R
    
    return [dPdt, dTdt, dRdt]

solution = odeint(PTR, y0, t, args=(k_p, q, k, m_t, m_r, v_p, v_t, v_r, f_r, f_t, d_t, d_r))

plt.figure(figsize=(8,6))
for i, label in enumerate(['P', 'T', 'R']):
    plt.semilogy(t, solution[:, i], label=f'{label}')
plt.xlabel('time(h)')
plt.ylabel('Populations of Different Species')
plt.title("q = 2")
plt.legend()
plt.show()
