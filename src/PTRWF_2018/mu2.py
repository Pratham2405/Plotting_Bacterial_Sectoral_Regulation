import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Parameters
k_F = 0.04  # h^-1 µm^3
F_e = 6e4  # µm^-3
K_F = k_F * F_e
k_P = 5e-3  # h^-1 µm^3
k_W = 1e-3  # h^-1 µm^3
k = 5e-4  # h^-1 µm^3
m_T = 500
m_R = 10000
theta_T = theta_R = theta = 4e14  # µm^-3
theta_W = 1e7  # µm^-3
v_P = v_T = v_R = v_W = v_F = 1e-8  # µm^3
V_c = 1.0  # µm^3

# Initial conditions
P0 = 1e6
T0 = 1e6
R0 = 1e6
W0 = 1e6
F0 = 1e6
y0 = [P0, T0, R0, W0, F0]

# Time points
t = np.linspace(0, 10, 1000)  # 0 to 10 hours

# External food concentration range
F_e_range = np.linspace(0, 1e6, 1000)  

# Define the system of ODEs
def model(y, t, K_F, k_P, k_W, k, m_T, m_R, theta_T, theta_R, theta_W, v_P, v_T, v_R, v_W, v_F):
    P, T, R, W, F = y
    V = v_P * P + v_T * T + v_R * R + v_W * W + v_F * F
    f_R = theta_R / (theta_R + pow((W/V), 2.0))
    f_T = pow((W/V), 2.0) / (theta_T + pow((W/V), 2.0))
    h_W = theta_W / (theta_W + (P/V))

    K_P = k_P * F / V
    K_W = k_W * F / V

    dPdt = K_P * T - k * R * P / V
    dTdt = (k / m_T) * (R * P / V) * f_T
    dRdt = (k / m_R) * (R * P / V) * f_R
    dWdt = K_W * T * h_W
    dFdt = K_F * T - K_P * T - K_W * T * h_W
    

    return [dPdt, dTdt, dRdt, dWdt, dFdt]

# Initialize array to store growth rates
mu_values = []

# Simulate for each F_e and calculate mu
for F_e in F_e_range:
    # Solve the system of ODEs
    K_F = k_F * F_e
    solution = odeint(model, y0, t, args=(K_F, k_P, k_W, k, m_T, m_R, theta_T, theta_R, theta_W, v_P, v_T, v_R, v_W, v_F))
    
    # Calculate growth rate (mu) at the final time point (or average over time)
    T = solution[:, 1]
    dTdt = np.gradient(T, t)  # Compute the derivative of P
    mu = dTdt[-1] / T[-1]  # Growth rate at the final time point
    mu_values.append(mu)

# Plot mu vs F_e
plt.figure(figsize=(8, 6))
plt.plot(F_e_range, mu_values, label='Growth Rate (mu)', color='b')
plt.xlabel('External Food Concentration [F_e]')
plt.ylabel('Growth Rate (mu)')
plt.title('Growth Rate as a Function of External Food Concentration')
plt.legend()
plt.grid(True)
plt.show()
