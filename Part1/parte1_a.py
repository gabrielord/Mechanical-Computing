from rungekutta import rungeKutta

# Constantes do sistema para o item a)

V = 50 / 3.6            # velocidade do carro [m/s]
k1 = k2 = 2.8*10**7     # constantes de mola das susp dianteira(1) e traseira(2) [N/m]
c1 = c2 = 3*10**4       # constantes de amortecimento das susp dianteira(1) e traseira(2) [kg/s]

"""## Exercício 1.a)"""

#Plota os gráficos de x, dx/dt, d(dx/dt)/dt, theta, dtheta/dt, d(dtheta/dt)/dt em t0 <= t <= tf e condicoes iniciais X_0
#com passo h = 0.0001 usando o Método de Runge-Kutta de 4a ordem
X, dX_dt, t = rungeKutta(h=0.0001, h_magnitude='Gráficos para h = 0.0001 (pequeno)', V=V, k=k1, c=c1)

#Plota os gráficos de x, dx/dt, d(dx/dt)/dt, theta, dtheta/dt, d(dtheta/dt)/dt em t0 <= t <= tf e condicoes iniciais X_0
#com passo h = 0.001 usando o Método de Runge-Kutta de 4a ordem
X, dX_dt, t = rungeKutta(h=0.001, h_magnitude='Gráficos para h = 0.001 (médio)', V=V, k=k1, c=c1)

#Plota os gráficos de x, dx/dt, d(dx/dt)/dt, theta, dtheta/dt, d(dtheta/dt)/dt em t0 <= t <= tf e condicoes iniciais X_0
#com passo h = 0.01 usando o Método de Runge-Kutta de 4a ordem
X, dX_dt, t = rungeKutta(h=0.01, h_magnitude='Gráficos para h = 0.01 (grande)', V=V, k=k1, c=c1)
