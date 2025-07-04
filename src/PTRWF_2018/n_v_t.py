import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Parameters
k_F = 0.04  # h^-1 µm^3
# F_e = 6e4  # µm^-3

k_P = 5e-3  # h^-1 µm^3
# k_W = 1e-3  # h^-1 µm^3
k = 5e-4  # h^-1 µm^3
m_T = 500
m_R = 10000
theta_T = theta_R = theta = 2e7  # µm^-3
theta_W = 1e7  # µm^-3
v_P = v_T = v_R = v_W = v_F = v = 1e-8  # µm^3
V_c = 1.0  # µm^3

# Initial conditions
P0 = 1e6
T0 = 1e6
R0 = 1e6
W0 = 1e6
F0 = 1e6
y0 = [P0, T0, R0, W0, F0]

# Time settings
t = np.linspace(0, 200, 100000)  # Time array

# Function to define the environment switching between two states
def environmenta(t):
    low_food = 1e4
    high_food = 10e4
    period = 20  # Period of environment fluctuation (20 hours per level)
    if t<20:
        return high_food
    elif (t // period) % 2 == 0:
        return high_food
    else:
        return low_food

def environmentb(t):
    low_food = 1e4
    high_food = 10e4
    period = 20  # Period of environment fluctuation (20 hours per level)
    if t<20:
        return high_food
    elif (t // period) % 2 == 0:
        return high_food
    else:
        return low_food

# PTRWF model equations
def PTRWF(y, t, regulated):
    P, T, R, W, F = y

    # External food concentration, which varies with time
    Fe = environmenta(t)

    # Growth rate dynamics
    K_F = k_F * Fe
    V = np.clip(v * (P + T + R + W), 1e-10, None)
    K_P = k_P * F/V
    h_W = theta_W / (theta_W + P/V)
    if regulated:
        # Regulated model where ppGpp affects ribosomal and transporter synthesis
        f_R = theta_R / (theta_R + (W/V))
        f_T = (W/V) / (theta_T + (W/V))
        k_W = 1e-3
    else:
        # Unregulated model with constant fR and fT
        f_R = 0.4915
        f_T = 1 - f_R
        k_W = 1e-7
    K_W = k_W * F/V
    # ODEs
    dPdt = K_P * T - k * R * P / V
    dTdt = (k / m_T) * (R * P / V) * f_T
    dRdt = (k / m_R) * (R * P / V) * f_R
    dWdt = K_W * T * h_W
    dFdt = K_F * T - K_P * T - K_W * T * h_W

    return [dPdt, dTdt, dRdt, dWdt, dFdt]

# Initial conditions vector
y0 = [P0, T0, R0, W0, F0]

# Solve the ODEs for both regulated and unregulated cases
sol_regulated = odeint(PTRWF, y0, t, args=(True,))
sol_unregulated = odeint(PTRWF, y0, t, args=(False,))

# Count the number of generations by integrating the growth rate
def count_generations(solution, P0):
    generations = []
    gen_count = 0
    P_prev = P0
    for P in solution[:, 0]:  # Loop over the precursor population (P)
        if P >= 2 * P_prev:    # If population doubles
            gen_count += 1     # Increment generation count
            P_prev = P         # Reset the previous population to the new level
        generations.append(gen_count)
    return generations

# Get the number of generations for both cases
generations_regulated = count_generations(sol_regulated, P0)
generations_unregulated = count_generations(sol_unregulated, P0)

# Extract environment (Fe) over time for plotting
Fe_values = np.array([environmenta(ti) for ti in t])

# Plot the results with secondary axis for environment
fig, ax1 = plt.subplots(figsize=(10, 6))

# Plot the number of generations on the left y-axis
ax1.plot(t, generations_regulated, label="Regulated Model", color='red')
ax1.plot(t, generations_unregulated, label="Unregulated Model", color='green')
ax1.set_xlabel('Time (hours)')
ax1.set_ylabel('Number of Generations')
ax1.set_title('Number of Generations in Fluctuating Environment')
ax1.legend(loc='upper left')
ax1.grid(True)

# Create a second y-axis for environment
ax2 = ax1.twinx()
ax2.plot(t, Fe_values, label="Environment (Food Conc.)", color='blue', linestyle='--')
ax2.set_ylabel('Food Concentration (Fe)', color='blue')

# Display the plot
plt.show()