# import numpy as np
# from scipy.integrate import odeint
# import matplotlib.pyplot as plt

# #Parameters
# k_p = 250
# k = 5e-4
# m_t = 500
# m_r = 10000
# v_p = v_t = v_r = 1e-8
# f_t = 0.50
# f_r = 0.50
# V_c = 1.0 #um^3
# d_t = 0.1 #as used in Jain_2016(h^-1)
# d_r = 0.0 #as used in Jain_2016(h^-1)

# #Initial Conditions
# P0 = 1e6
# T0 = 1e6
# R0 = 1e6
# y0 = [P0, T0, R0]

# t = np.linspace(0,100,100)
# f_r_range = np.linspace(0,1.0, 100)
# q_range = np.linspace(0, 100, 100)


# def PTRD(y, t, k_p, q, k, m_t, m_r, v_p, v_t, v_r, f_r, f_t, d_t, d_r):
#     P, T, R = y
#     V = v_p * P + v_t * T + v_r * R
#     K_p = k_p * q
#     dPdt = K_p * T - k * R * P / V
#     dTdt = (k / m_t) * (R * P / V) * f_t - d_t * T
#     dRdt = (k / m_r) * (R * P / V) * f_r - d_r * R
    
#     return [dPdt, dTdt, dRdt]

# # def div_sim(y0, t, V_c):
# #     y = y0
# #     results = [y0]
# #     for i in range(1, len(t)):
# #         y = odeint(PTRD, y, [t[i-1], t[i]], args=(k_p, q, k, m_t, m_r, v_p, v_t, v_r, f_r, f_t, d_t, d_r))[-1]
# #         V = v_p * y[0] + v_t * y[1] + v_r * y[2]
# #         if V >= V_c:
# #             y = y/2
# #         results.append(y)
# #     return np.array(results)




# for q in q_range:
#     mu_values = []
#     for f_r in f_r_range:
#         f_t = 1- f_r
#         print('presim')
#         solution = odeint(PTRD, y0, t, args=(k_p, q, k, m_t, m_r, v_p, v_t, v_r, f_r, f_t, d_t, d_r))
#         print('postsim')
#         print(q)
#         P = solution[:, 0]
#         coeff = np.polyfit(t, np.log(P), 1)
#         mu = coeff[0]
#         # dPdt = np.gradient(P, t)
#         # mu = dPdt[-1]/P[-1]
#         print(mu)
#         mu_values.append(mu)
#     dmudfr = np.gradient(mu_values, f_r_range)
#     for i in range(len(dmudfr)):
#         if dmudfr[i] < 0:
#             mu_max = mu_values
#         break


# plt.figure(figsize=(8,6))
# plt.plot(q_range, mu_max)
# plt.xlabel('q')
# plt.ylabel('Mu_max')
# plt.title("Replication of Monod's Law")
# plt.legend()
# plt.show()

# import numpy as np
# from scipy.integrate import odeint
# import matplotlib.pyplot as plt

# # Parameters
# k_p = 250
# k = 5e-4
# m_t = 500
# m_r = 10000
# v_p = v_t = v_r = 1e-8
# f_t = 0.50
# f_r = 0.50
# V_c = 1.0 # um^3
# d_t = 0.1 # as used in Jain_2016(h^-1)
# d_r = 0.0 # as used in Jain_2016(h^-1)

# # Initial Conditions
# P0 = 1e6
# T0 = 1e6
# R0 = 1e6
# y0 = [P0, T0, R0]

# t = np.linspace(0, 100, 100)
# f_r_range = np.linspace(0, 1.0, 100)
# q_range = np.linspace(0, 100, 100)

# def PTRD(y, t, k_p, q, k, m_t, m_r, v_p, v_t, v_r, f_r, f_t, d_t, d_r):
#     P, T, R = y
#     V = v_p * P + v_t * T + v_r * R
#     K_p = k_p * q
#     dPdt = K_p * T - k * R * P / V
#     dTdt = (k / m_t) * (R * P / V) * f_t - d_t * T
#     dRdt = (k / m_r) * (R * P / V) * f_r - d_r * R
    
#     return [dPdt, dTdt, dRdt]

# mu_max_list = []

# for q in q_range:
#     mu_values = []
#     for f_r in f_r_range:
#         f_t = 1 - f_r
#         solution = odeint(PTRD, y0, t, args=(k_p, q, k, m_t, m_r, v_p, v_t, v_r, f_r, f_t, d_t, d_r))
#         P = solution[:, 0]
#         coeff = np.polyfit(t, np.log(P), 1)
#         mu = coeff[0]
#         mu_values.append(mu)
    
#     # Find maximum mu for this q value and store it in mu_max_list
#     mu_max_list.append(max(mu_values))

# # Plotting mu_max vs q
# plt.figure(figsize=(8,6))
# plt.plot(q_range, mu_max_list)
# plt.xlabel('q')
# plt.ylabel('Mu_max')
# plt.title("Replication of Monod's Law")
# plt.legend(['Mu_max'])
# plt.show()
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Parameters
k_p = 250
k = 5e-4
m_t = 500
m_r = 10000
v_p = v_t = v_r = 1e-8
d_t = 0.1 # as used in Jain_2016(h^-1)
d_r = 0.0 # as used in Jain_2016(h^-1)

# Initial Conditions
P0 = 1e6
T0 = 1e6
R0 = 1e6
y0 = [P0, T0, R0]

t_span = (0, 50)
t_eval = np.linspace(0, 50, 50)
f_r_range = np.linspace(0, 1.0, 50)
q_range = np.linspace(0, 100, 50)

def PTRD(t, y, k_p, q, k, m_t, m_r, v_p, v_t, v_r, f_r, f_t, d_t, d_r):
    P, T, R = y
    V = v_p * P + v_t * T + v_r * R
    if V <= 0:
        V = 1e-10
    K_p = k_p * q
    dPdt = K_p * T - k * R * P / V
    dTdt = (k / m_t) * (R * P / V) * f_t - d_t * T
    dRdt = (k / m_r) * (R * P / V) * f_r - d_r * R

    # print(f"t={t}, P={P}, T={T}, R={R}, V={V}, dPdt={dPdt}, dTdt={dTdt}, dRdt={dRdt}")
    
    return [dPdt, dTdt, dRdt]

mu_max_list = []

for q in q_range:
    mu_values = []
    for f_r in f_r_range:
        f_t = 1 - f_r
        sol = solve_ivp(PTRD, t_span, y0,
                        args=(k_p, q, k, m_t, m_r, v_p, v_t, v_r, f_r, f_t, d_t, d_r),
                        t_eval=t_eval,
                        method='BDF', rtol = 1e-3, atol = 1e-6) # Use a robust solver
        
        if sol.success:
            P = sol.y[0]
            if np.all(P > 0): 
                coeff = np.polyfit(t_eval, np.log(P), 1)
                mu_values.append(coeff[0])
    
    if mu_values: 
        mu_max_list.append(max(mu_values))
    else:
        mu_max_list.append(np.nan) 

# Plotting mu_max vs q
plt.figure(figsize=(8,6))
plt.plot(q_range, mu_max_list)
# plt.plot(f_r +  )
plt.xlabel('q')
plt.ylabel('Mu_max')
plt.title("Replication of Monod's Law")
plt.legend(['Mu_max'])
plt.show()
