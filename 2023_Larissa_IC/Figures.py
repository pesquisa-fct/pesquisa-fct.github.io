import matplotlib.pyplot as plt
import numpy as np

data4Lx4L_M1_Re1 = np.loadtxt('../Results/4Lx4L_M1_Re1/asymptotics.dat')

data4Lx4L_M2_Re1 = np.loadtxt('../Results/4Lx4L_M2_Re1/asymptotics.dat')

dataSlope4Lx4L_M1_Re1 = np.loadtxt('../Slopes/Slopes4Lx4L_M1_Re1.dat')


# Extract the columns you want to plot
r_col = 0  # Change this to the index of the column you want to use as x-axis
u_col = 1  # Change this to the index of the column you want to use as y-axis
v_col = 2
p_col = 3
psi_col = 4
omega_col = 5

#4Lx4L, M1, Re=1 - Slope
r_valSlope4Lx4L_M1_Re1 = dataSlope4Lx4L_M1_Re1[:, r_col]
u_valSlope4Lx4L_M1_Re1 = dataSlope4Lx4L_M1_Re1[:, u_col]
v_valSlope4Lx4L_M1_Re1 = dataSlope4Lx4L_M1_Re1[:, v_col]
p_valSlope4Lx4L_M1_Re1 = dataSlope4Lx4L_M1_Re1[:, p_col]
psi_valSlope4Lx4L_M1_Re1 = dataSlope4Lx4L_M1_Re1[:, psi_col]
omega_valSlope4Lx4L_M1_Re1 = dataSlope4Lx4L_M1_Re1[:, omega_col]

#4Lx4L, M1, Re=1
r_val4Lx4L_M1_Re1 = np.log10(abs(data4Lx4L_M1_Re1[:, r_col]))
u_val4Lx4L_M1_Re1 = np.log10(abs(data4Lx4L_M1_Re1[:, u_col]))
v_val4Lx4L_M1_Re1 = np.log10(abs(data4Lx4L_M1_Re1[:, v_col]))
p_val4Lx4L_M1_Re1 = np.log10(abs(data4Lx4L_M1_Re1[:, p_col]))
psi_val4Lx4L_M1_Re1 = np.log10(abs(data4Lx4L_M1_Re1[:, psi_col]))
omega_val4Lx4L_M1_Re1 = np.log10(abs(data4Lx4L_M1_Re1[:, omega_col]))

#4Lx4L, M2, Re=1
r_val4Lx4L_M2_Re1 = np.log10(abs(data4Lx4L_M2_Re1[:, r_col]))
u_val4Lx4L_M2_Re1 = np.log10(abs(data4Lx4L_M2_Re1[:, u_col]))
v_val4Lx4L_M2_Re1 = np.log10(abs(data4Lx4L_M2_Re1[:, v_col]))
p_val4Lx4L_M2_Re1 = np.log10(abs(data4Lx4L_M2_Re1[:, p_col]))
psi_val4Lx4L_M2_Re1 = np.log10(abs(data4Lx4L_M2_Re1[:, psi_col]))
omega_val4Lx4L_M2_Re1 = np.log10(abs(data4Lx4L_M2_Re1[:, omega_col]))


#(u,v)
plt.scatter(r_val4Lx4L_M1_Re1, u_val4Lx4L_M1_Re1, color='red', marker='o', label=r'$M_{1}$')
plt.scatter(r_val4Lx4L_M2_Re1, u_val4Lx4L_M2_Re1, color='blue', marker='d', label=r'$M_{2}$')
plt.plot(r_valSlope4Lx4L_M1_Re1, u_valSlope4Lx4L_M1_Re1, color='black', label=r'$r^{0.5445}$')

plt.scatter(r_val4Lx4L_M1_Re1, v_val4Lx4L_M1_Re1, color='red', marker='o')
plt.scatter(r_val4Lx4L_M2_Re1, v_val4Lx4L_M2_Re1, color='blue', marker='d')
plt.plot(r_valSlope4Lx4L_M1_Re1, v_valSlope4Lx4L_M1_Re1, color='black')

plt.text(-1.1, 0.0, r'$\log\left( u \right)$', fontsize=12, ha='center', va='center', rotation=0)
plt.text(-1.2, -1.1, r'$\log\left( v \right)$', fontsize=12, ha='center', va='center', rotation=0)

plt.xlabel(r'$\log\left( r \right)$', fontsize=12)
plt.legend(fontsize=12)
plt.legend(loc='lower left')

plt.grid(False)
plt.savefig('v.pdf')
plt.show()


