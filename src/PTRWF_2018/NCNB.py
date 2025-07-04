import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Parameters
k_F = 0.04  # h^-1 µm^3
F_e = 6e4  # µm^-3
K_F = k_F * F_e
k_P = 5e-3  # h^-1 µm^3
k_W = 0  # h^-1 µm^3
k = 5e-4  # h^-1 µm^3
m_T = 500
m_R = 10000
f_R = 0.4915
f_T = 0.5085
theta_T = theta_R = theta = 2e7  # µm^-3
theta_W = 1e7  # µm^-3
v_P = v_T = v_R = v_W = v_F = 1e-8  # µm^3

# Initial conditions
P0 = 1e4
T0 = 1e4
R0 = 1e4
W0 = 1e4
F0 = 1e4
y0 = [P0, T0, R0, W0, F0]

# Time points
t = np.linspace(0, 10, 1000)  # Adjusted to 6 hours
F_e_range = np.linspace(0,1e5, 1000)

# Define the system of ODEs
def PTRWFwD(y, t, K_F, k_P, k_W, k, m_T, m_R, theta_T, theta_R, theta_W, v_P, v_T, v_R, v_W, v_F):
    P, T, R, W, F = y
    V = v_P * P + v_T * T + v_R * R + v_W * W + v_F * F
    # f_R = theta_R / (theta_R + W/V)
    # f_T = (W/V) / (theta_T + W/V)
    h_W = theta_W / (theta_W + P/V)

    K_F = k_F * F_e
    K_P = k_P * F / V
    K_W = k_W * F / V

    dPdt = K_P * T - k * R * P / V
    dTdt = (k / m_T) * (R * P / V) * f_T
    dRdt = (k / m_R) * (R * P / V) * f_R
    dWdt = K_W * T * h_W
    dFdt = K_F * T - K_P * T - K_W * T * h_W
    

    return [dPdt, dTdt, dRdt, dWdt, dFdt]

mu_values = []

for F_e in F_e_range:
    print('presim')
    solution = odeint(PTRWFwD, y0, t, args=(K_F, k_P, k_W, k, m_T, m_R, theta_T, theta_R, theta_W, v_P, v_T, v_R, v_W, v_F))
    print('postsim')
    P = solution[:,0]
    coeff = np.polyfit(t, np.log(P), 1)
    mu = coeff[0]
    print(mu)
    mu_values.append(mu)




plt.figure(figsize=(8, 6))
plt.plot(F_e_range, mu_values)
plt.xlabel('F_e')
plt.ylabel('Mu')
plt.title('Growth Laws(PTRWFwD)')
plt.legend()
plt.grid(True, which="both", ls="--")
plt.show()

#Tested OK

