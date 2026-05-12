import numpy as np
import math
import matplotlib.pyplot as plt
import tensorflow as tf

def SolucaoAnalitica(x):
    return (-math.sin(math.pi * x) / (math.pi**2))

def SolucaoNumerica(a, b, ya, yb, x, N):
    Y = np.zeros(N+1)
    Y[0] = ya
    Y[N] = yb
    h = (b-a)/N
    A = np.zeros((N-1, N-1))
    B = np.zeros(N-1)
    for i in range(N-1):
        n=i+1
        A[i][i] = -2.
        B[i] = h*h*math.sin(math.pi * x[n])
        if i != 0:
            A[i][i-1] = 1.
        else:
            B[i] += -Y[0]
        if i != N-2:
            A[i][i+1] = 1.
        else:
            B[i] += -Y[N]
    sol = np.linalg.solve(A,B)
    for i in range(N-1):
        n=i+1
        Y[n] = sol[i]
    return Y


# === CONFIGURAÇÕES ===
tf.random.set_seed(0)
N_f = 200    # pontos internos (collocation)
N_b = 2      # pontos de contorno (x=0 e x=1)

# === DEFINIR O PDE ===
# u_xx = sin(pi x)

# ==============================
# 1) PONTOS DO DOMÍNIO
# ==============================
x_f = tf.random.uniform((N_f,1), minval=0.0, maxval=1.0)

# Pontos de contorno
x_b = tf.constant([[0.0],[1.0]], dtype=tf.float32)
y_b = tf.constant([[0.0],[0.0]], dtype=tf.float32)


# ==============================
# 2) REDE NEURAL
# ==============================
# === ARQUITETURA DA PINN ===
model = tf.keras.Sequential([
    tf.keras.layers.Dense(32, activation='tanh', input_shape=(1,)),
    tf.keras.layers.Dense(32, activation='tanh'),
    tf.keras.layers.Dense(1, activation='linear')
])

# === FUNÇÃO DE TREINO ===
optimizer = tf.keras.optimizers.Adam(learning_rate=0.001)

# ==============================
# 3) FUNÇÃO LOSS
# ==============================
@tf.function
def compute_loss():
    # PDE residual
    with tf.GradientTape(persistent=True) as tape2:
        tape2.watch(x_f)
        with tf.GradientTape() as tape1:
            tape1.watch(x_f)
            y = model(x_f)
        y_x = tape1.gradient(y, x_f)
    y_xx = tape2.gradient(y_x, x_f)

    f = tf.sin(np.pi * x_f)  # lado direito

    loss_pde = tf.reduce_mean((y_xx - f)**2)

    # Boundary loss
    y_pred_b = model(x_b)
    loss_bc = tf.reduce_mean((y_pred_b - y_b)**2)

    return loss_pde + loss_bc

@tf.function
def train_step():
    with tf.GradientTape() as tape:
        loss = compute_loss()
    grads = tape.gradient(loss, model.trainable_variables)
    optimizer.apply_gradients(zip(grads, model.trainable_variables))
    return loss

# ==============================
# 4) TREINAMENTO
# ==============================
# === TREINAMENTO ===
epochs = 5000
for epoch in range(1, epochs+1):
    loss_value = train_step()
    if epoch % 500 == 0:
        print(f"Epoch {epoch}, Loss = {loss_value.numpy():.6e}")

# ==============================
# 5) PLOTS
# ==============================
# === SOLUÇÃO ANALÍTICA ===
x_plot = np.linspace(0,1,200).reshape(-1,1)
#u_exact = np.sin(np.pi*x_plot)/np.pi**2

y_pred = model(tf.constant(x_plot, dtype=tf.float32))
x_plot = tf.cast(x_plot, tf.float32)
#x_plot = tf.expand_dims(x_plot, 1)
#print(y_pred)
#print(x_plot)

a=0
b=1
ya=0
yb=0
N = 200
x = np.linspace(a, b, N+1)
y_exact = np.zeros(N+1)
#Calcula sol. analitica
for i in range(N+1):
    y_exact[i] = SolucaoAnalitica(x[i])

Y = SolucaoNumerica(a, b, ya, yb, x, N)

# === PLOT ===
plt.figure(figsize=(8,4))
plt.scatter(x, Y, marker='o', color='red', label='Numerical')
plt.scatter(x_plot, y_pred, marker='x', color='blue', label='PINN')
plt.plot(x, y_exact, color='black', label="Analytical")
plt.legend()
plt.show()
