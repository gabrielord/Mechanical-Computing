"""## Exercício 1.b)"""

from rungekutta import rungeKutta

vector_V = [30, 70]                     #diferentes velocidades do carro [km/h]
vector_C = [10**3, 5*10**4, 2.5*10**5]  #diferentes coeficientes de amortecimento da suspensao [kg/s]
vector_K = [2*10**4, 5*10**6, 7*10**8]  #diferentes constantes de mola da suspensao [N/m]



#Plota os gráficos de x, dx/dt, d(dx/dt)/dt, theta, dtheta/dt, d(dtheta/dt)/dt em t0 <= t <= tf com condicoes iniciais X_0 e velocidade V = 30 e 70 km/h
#com passo h = 0.0001 usando o Método de Runge-Kutta de 4a ordem 
for v in vector_V:
  X, dX_dt, t = rungeKutta(h=0.0001, h_magnitude=f'Gráficos para V = {v} km/h', V=v, k=2.8*10**7, c=3*10**4)

#Plota os gráficos de x, dx/dt, d(dx/dt)/dt, theta, dtheta/dt, d(dtheta/dt)/dt em t0 <= t <= tf com condicoes iniciais X_0 e
#coeficientes de amortecimento c = 10^3 e 5*10^4 e 2,5*10^5 kg/s com passo h = 0.0001 usando o Método de Runge-Kutta de 4a ordem 
for c in vector_C:
  X, dX_dt, t = rungeKutta(h=0.0001, h_magnitude=f'Gráficos para c = {"{:.2e}".format(c)} kg/s', V=50 / 3.6, k=2.8*10**7, c=c)

#Plota os gráficos de x, dx/dt, d(dx/dt)/dt, theta, dtheta/dt, d(dtheta/dt)/dt em t0 <= t <= tf com condicoes iniciais X_0 e
#constantes de mola k = 2*10^4 e 5*10^6 e 7*10^8 N/m com passo h = 0.0001 usando o Método de Runge-Kutta de 4a ordem 
for k in vector_K:
  X, dX_dt, t = rungeKutta(h=0.0001, h_magnitude=f'Gráficos para k = {"{:.2e}".format(k)} N/m', V=50 / 3.6, k=k, c=3*10**4)