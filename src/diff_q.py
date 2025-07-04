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
V_c = 1.0 #um^3
d_t = 0.1 #as used in Jain_2016(h^-1)
d_r = 0.0 #as used in Jain_2016(h^-1)

#Initial Conditions
P0 = 1e6
T0 = 1e6
R0 = 1e6
y0 = [P0, T0, R0]

# Time vector
t = np.linspace(0, 100, 1000)
f_R_range = np.linspace(0, 1.0, 100)
# F_e_range = np.linspace(0, 1e5, 5)  # Using fewer points for clear plotting
q_range = np.linspace(1,10,10)

# Function to compute growth rate
def PTR(y, t, k_p, q, k, m_t, m_r, v_p, v_t, v_r, f_R, d_t, d_r):
    P, T, R = y
    V = v_p * P + v_t * T + v_r * R
    K_p = k_p * q
    dPdt = K_p * T - k * R * P / V
    dTdt = (k / m_t) * (R * P / V) * (1-f_R) - d_t * T
    dRdt = (k / m_r) * (R * P / V) * f_R - d_r * R
    
    return [dPdt, dTdt, dRdt]


# Arrays to store results
optimal_mu_values = []  # To store the max mu for each F_e

# Create a figure to plot multiple lines
plt.figure(figsize=(10, 6))

# Loop over F_e_range to compute mu vs f_R for each F_e
for q in q_range:
    mu_values = []

    # Compute mu for each f_R in f_R_range
    for f_R in f_R_range:
        solution = odeint(PTR, y0, t, args=(k_p, q, k, m_t, m_r, v_p, v_t, v_r, f_R, d_t, d_r))
        R = solution[:, 2]  # Ribosomal concentration
        log_R = np.log(np.clip(R, 1e-6, None))  # Avoid log(0)
        coeff = np.polyfit(t, log_R, 1)  # Linear fit to log of ribosomes to get growth rate
        mu = coeff[0]
        mu_values.append(mu)

    # Plot mu vs f_R for this F_e
    plt.plot(f_R_range, mu_values, label=f'q = {q}')

    # Store the maximum mu for this F_e
    optimal_mu_values.append(max(mu_values))

print(optimal_mu_values)
# Customize the plot for mu vs f_R_range
plt.xlabel('Ribosomal Fraction (f_R)')
plt.ylabel('Growth Rate (μ)')
plt.title('μ vs. f_R for Different Medium Quality(q)')
plt.legend()
plt.grid(True)
plt.show()

# # Now plot μ vs. F_e_range
# plt.figure(figsize=(10, 6))
# plt.plot(f_R_range, optimal_mu_values, marker='o', linestyle='-', color='b', label='Optimal Growth Rate (μ)')
# plt.xlabel('External Food Concentration (F_e)')
# plt.ylabel('Optimal Growth Rate (μ)')
# plt.title('Optimal Growth Rate (μ) vs. External Food Concentration (F_e)')
# plt.grid(True)
# plt.legend()
# plt.show()