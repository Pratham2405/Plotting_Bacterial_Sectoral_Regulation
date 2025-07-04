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
theta_T = theta_R = theta = 4.5e3  # µm^-3
theta_W = 1e7  # µm^-3
v_P = v_T = v_R = v_W = v_F = 1e-8  # µm^3
V_c = 1.0  # µm^3

# Initial conditions
P0 = 1e4
T0 = 1e4
R0 = 1e4
W0 = 1e4
F0 = 1e4
y0 = [P0, T0, R0, W0, F0]

# Time points
t = np.linspace(0, 10, 1000)  # Adjusted to 6 hours

# Define the system of ODEs
def model(y, t, K_F, k_P, k_W, k, m_T, m_R, theta_T, theta_R, theta_W, v_P, v_T, v_R, v_W, v_F):
    P, T, R, W, F = y
    V = v_P * P + v_T * T + v_R * R + v_W * W + v_F * F
    f_R = theta_R / (theta_R + pow((W/V), 0.5))
    f_T = pow((W/V), 0.5) / (theta_T + pow((W/V), 0.5))
    h_W = theta_W / (theta_W + (P/V))

    K_P = k_P * F / V
    K_W = k_W * F / V

    dPdt = K_P * T - k * R * P / V
    dTdt = (k / m_T) * (R * P / V) * f_T
    dRdt = (k / m_R) * (R * P / V) * f_R
    dWdt = K_W * T * h_W
    dFdt = K_F * T - K_P * T - K_W * T * h_W
    

    return [dPdt, dTdt, dRdt, dWdt, dFdt]

# Simulation with division
def simulate_with_division(y0, t, V_c):
    y = np.array(y0)
    results = [y0]
    for i in range(1, len(t)):
        y = odeint(model, y, [t[i-1], t[i]], args=(K_F, k_P, k_W, k, m_T, m_R, theta_T, theta_R, theta_W, v_P, v_T, v_R, v_W, v_F))[-1]
        V = v_P * y[0] + v_T * y[1] + v_R * y[2] + v_W * y[3] + v_F * y[4]
        if V >= V_c:
            y /= 2  # Halve the populations
        results.append(y)
    return np.array(results)

# Run simulation
solution = simulate_with_division(y0, t, V_c)
P = solution[:,0]
T = solution[:,1]
R = solution[:,2]
W = solution[:,3]
F = solution[:,4]
V = v_P * P + v_T * T + v_R * R + v_W * W + v_F * F
f_R = theta_R / (theta_R + W/V)
max_fr = max(f_R)
print(max_fr)
ss_fr = np.polyfit(t[200:999], f_R[200:999], 1)
coeff = ss_fr[1]
print(coeff)


# Plot results
plt.figure(figsize=(8, 6))
for i, label in enumerate(['P', 'T', 'R', 'W', 'F']):
    plt.semilogy(t, solution[:, i], label=f'{label}')
plt.xlabel('Time (hours)')
plt.ylabel('Populations')
plt.title('Dynamics of Molecular Pools in a Bacterial Cell with Division')
# plt.plot(t, f_R, color= 'b')
# plt.plot(t, [0.4949] * len(t), color='red', linestyle='--')
plt.legend()
plt.grid(True, which="both", ls="--")
plt.show()