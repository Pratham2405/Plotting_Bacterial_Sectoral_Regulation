import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Parameters common across models
k_P = 5e-3  # h^-1 µm^3
k_W = 1e-3  # h^-1 µm^3 for some models
k = 5e-4  # h^-1 µm^3
m_T = 500
m_R = 10000
theta_T = theta_R = theta = 2e7  # µm^-3
theta_W = 1e7  # µm^-3
v_P = v_T = v_R = v_W = v_F = 1e-8  # µm^3

# Time points and ranges
t = np.linspace(0, 10, 1000)  # Time vector
f_e_range = np.linspace(0, 1e5, 50)  # External food range (F_e)

# Initial conditions
P0 = 1e4
T0 = 1e4
R0 = 1e4
W0 = 1e4
F0 = 1e4
y0 = [P0, T0, R0, W0, F0]


# Model 1 (cb.py)
def model_cb(y, t, K_F, k_P, k_W, k, m_T, m_R, theta_T, theta_R, theta_W, v_P, v_T, v_R, v_W, v_F):
    P, T, R, W, F = y
    V = v_P * P + v_T * T + v_R * R + v_W * W + v_F * F
    f_R = theta_R / (theta_R + W / V)
    f_T = (W / V) / (theta_T + W / V)
    h_W = theta_W / (theta_W + P / V)
    K_P = k_P * F / V
    K_W = k_W * F / V
    dPdt = K_P * T - k * R * P / V
    dTdt = (k / m_T) * (R * P / V) * f_T
    dRdt = (k / m_R) * (R * P / V) * f_R
    dWdt = K_W * T * h_W
    dFdt = K_F * T - K_P * T - K_W * T * h_W
    return [dPdt, dTdt, dRdt, dWdt, dFdt]


# Model 2 (cnb.py)
def model_cnb(y, t, K_F, k_P, k_W, k, m_T, m_R, theta_T, theta_R, theta_W, v_P, v_T, v_R, v_W, v_F):
    P, T, R, W, F = y
    V = v_P * P + v_T * T + v_R * R + v_W * W + v_F * F
    f_R = 0.4915  # Fixed ribosomal fraction
    f_T = 1 - f_R
    h_W = 0.1219
    K_P = k_P * F / V
    K_W = k_W * F / V
    dPdt = K_P * T - k * R * P / V
    dTdt = (k / m_T) * (R * P / V) * f_T
    dRdt = (k / m_R) * (R * P / V) * f_R
    dWdt = K_W * T * h_W
    dFdt = K_F * T - K_P * T - K_W * T * h_W
    return [dPdt, dTdt, dRdt, dWdt, dFdt]


# Model 3 (ncnb.py, simplified for this example)
def model_ncnb(y, t, K_F, k_P, k_W, k, m_T, m_R, theta_T, theta_R, theta_W, v_P, v_T, v_R, v_W, v_F):
    P, T, R, W, F = y
    V = v_P * P + v_T * T + v_R * R + v_W * W + v_F * F
    f_R = 0.4915  # Simplified
    f_T = 1 - f_R
    k_W = 0
    h_W = theta_W / (theta_W + P / V)
    K_P = k_P * F / V
    K_W = k_W * F / V
    dPdt = K_P * T - k * R * P / V
    dTdt = (k / m_T) * (R * P / V) * f_T
    dRdt = (k / m_R) * (R * P / V) * f_R
    dWdt = K_W * T * h_W
    dFdt = K_F * T - K_P * T - K_W * T * h_W
    return [dPdt, dTdt, dRdt, dWdt, dFdt]


# Model 4 (new code, unchanged with dynamic f_R)
def model_new(y, t, K_F, k_P, k_W, k, m_T, m_R, theta_W, v_P, v_T, v_R, v_W, v_F, f_R):
    P, T, R, W, F = y
    V = v_P * P + v_T * T + v_R * R + v_W * W + v_F * F
    f_T = 1 - f_R  # Constraint f_R + f_T = 1
    h_W = theta_W / (theta_W + P / V)
    K_P = k_P * F / V
    K_W = k_W * F / V
    dPdt = K_P * T - k * R * P / V
    dTdt = (k / m_T) * (R * P / V) * f_T
    dRdt = (k / m_R) * (R * P / V) * f_R
    dWdt = K_W * T * h_W
    dFdt = K_F * T - K_P * T - K_W * T * h_W
    return [dPdt, dTdt, dRdt, dWdt, dFdt]


# Function to calculate growth rate (mu) for each model
def calculate_mu(model_func, f_e_range, extra_args=None):
    mu_values = []
    for F_e in f_e_range:
        K_F = 0.04 * F_e
        if extra_args is None:
            args = (K_F, k_P, k_W, k, m_T, m_R, theta_T, theta_R, theta_W, v_P, v_T, v_R, v_W, v_F)
        else:
            args = (K_F, k_P, k_W, k, m_T, m_R, theta_W, v_P, v_T, v_R, v_W, v_F) + extra_args

        solution = odeint(model_func, y0, t, args=args)
        P_sol = solution[:, 0]
        coeff_R = np.polyfit(t, np.log(P_sol), 1)  # Logarithmic growth rate (mu)
        slope = coeff_R[0]
        mu_values.append(slope)
    return mu_values


# Calculate mu for all four models
mu_cb = calculate_mu(model_cb, f_e_range)
mu_cnb = calculate_mu(model_cnb, f_e_range)
mu_ncnb = calculate_mu(model_ncnb, f_e_range)
# Calculate mu for the new model without fixing f_R
f_R_values = np.linspace(0.1, 0.9, 5)  # Use a range of f_R values for the new model
mu_new = calculate_mu(model_new, f_e_range, extra_args=(f_R_values[0],))  # Example with f_R = 0.1

# Plotting the results for all models
plt.figure(figsize=(10, 6))
plt.plot(f_e_range, mu_cb, color='r', label='Model CB')
plt.plot(f_e_range, mu_cnb, color='#ADD8E6', label='Model CNB')
plt.plot(f_e_range, mu_ncnb, color='g', label='Model NCNB')
plt.plot(f_e_range, mu_new, color='#00008B', label='Model NCB')
plt.xlabel('Fe')
plt.ylabel('μ')
plt.title('Growth Rate μ over Fe for Different Models')
plt.legend()
plt.grid(True)
plt.show()
