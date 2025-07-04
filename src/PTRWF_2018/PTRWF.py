import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Parameters
k_F = 0.04  # h^-1 µm^3
F_e = 6e4  # µm^-3
K_F = k_F * F_e
k_P = 5e-3  # h^-1 µm^3
k_W = 1e-3  # h^-1 µm^3
k = 5e-4  # h^-1 µm^3
m_T = 500
m_R = 10000
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
t = np.linspace(0, 30, 1000)  # Adjusted to 6 hours

# Define the system of ODEs
def model(y, t, K_F, k_P, k_W, k, m_T, m_R, theta_T, theta_R, theta_W, v_P, v_T, v_R, v_W, v_F):
    P, T, R, W, F = y
    V = v_P * P + v_T * T + v_R * R + v_W * W + v_F * F
    f_R = theta_R / (theta_R + W/V)
    f_T = (W/V) / (theta_T + W/V)
    h_W = theta_W / (theta_W + P/V)

    K_P = k_P * F / V
    K_W = k_W * F / V

    dPdt = K_P * T - k * R * P / V
    dTdt = (k / m_T) * (R * P / V) * f_T
    dRdt = (k / m_R) * (R * P / V) * f_R
    dWdt = K_W * T * h_W
    dFdt = K_F * T - K_P * T - K_W * T * h_W
    

    return [dPdt, dTdt, dRdt, dWdt, dFdt]

solution = odeint(model, y0, t, args=(K_F, k_P, k_W, k, m_T, m_R, theta_T, theta_R, theta_W, v_P, v_T, v_R, v_W, v_F))

plt.figure(figsize=(8, 6))
for i, label in enumerate(['P', 'T', 'R', 'W', 'F']):
    plt.semilogy(t, solution[:, i], label=f'{label}')
plt.xlabel('Time (hours)')
plt.ylabel('Populations')
plt.title('Dynamics of Molecular Pools in a Bacterial Cell with Division')
plt.legend()
plt.grid(True, which="both", ls="--")
plt.show()
