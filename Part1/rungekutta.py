# -*- coding: utf-8 -*-
"""EP1.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1e1dDBQA7KardqJcfNJ0QsF1cVPDie3VJ

<h1>1° Exercício Programa de PMR 3401</h1>
Gabriel Souza Lima - NUSP: <br>
Vitor Vac Bitu Alves - NUSP: 11913833<br>
<h2>Parte 1</h2>
<h3>Método de Runge Kutta</h3>
"""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (20,12)

# Constantes do sistema

me = 20                 # massa do virabrequim [kg]
r = 0.045               # raio de giro do virabrequim [m]
fe = 2100               # frequencia de giro virabrequim [rpm]
we = fe * 2*np.pi / 60  # frequencia de giro virabrequim [rad/s]
M = 1783                # massa do carro [kg]
a = 1220/1000           # distancia da susp dianteira para o CM [m]
b = 1500/1000           # distancia da susp traseira para o CM [m]
Ic = 4000               # momento de inercia do carro [kg/m^2]
e = 0.75                # distancia do virabrequim para o CM [m]
L = 0.5                 # distancia entre picos das oscilacoes do solo [m]
A = 60/1000             # amplitude das oscilacoes do solo [m]
f = 0.35                # altura do virabrequim em relacao ao CM [m]
Fn = me*(we**2)*r       # forca transmitida pelo virabrequim a estrutura [N]
t0 = 0                                                  # instante inicial [s]
tf = 4                                                  # instante final [s]

# Condicoes iniciais do problema

theta_0 = 0.09                                          # rotacao inicial do carro em torno do CM [rad]
x_0 = x_dot_0 = theta_dot_0 = 0                         # posicao, velocidade e velocidade angular iniciais do carro [m][m/s][rad/s]
X_0 = np.array([x_0, x_dot_0, theta_0, theta_dot_0])    # vetor X_0 das condicoes iniciais do sistema [m][m/s][rad][rad/s]

#Funcao da entrada no pneu dianteiro
def func_d1(t:float, w:float, A=A):
  if 0 < t < 2:
    d1 =  A*(1 - np.cos(w*t))
  else:
    d1 = 0
  return d1

#Funcao da entrada no pneu traseiro
def func_d2(t:float, w:float, A=A):
  if 0 < t < 2:
    d2 =  A*(1 + np.cos(w*t))
  else:
    d2 = 0
  return d2

#Derivada da entrada no pneu dianteiro
def func_d1_dot(t:float, w:float, A=A):
  if 0 < t < 2:
    d1_dot =  A*w*np.sin(w*t)
  else:
    d1_dot = 0
  return d1_dot

#Derivada da entrada no pneu traseiro
def func_d2_dot(t:float, w:float, A=A):
  if 0 < t < 2:
    d2_dot = -A*w*np.sin(w*t)
  else:
    d2_dot = 0
  return d2_dot

#Funcao f(t,X(t)) = dX(t)/dt
def func_dX_dt(t:float, X:np.array, V:float, k:float, c:float, a=a, b=b, Fn=Fn, M=M, f=f, Ic=Ic, L=L):
  #k = constante de mola da suspensao [N/m]
  k1 = k2 = k
  #c = constante de amortecimento da suspensao [kg/s]
  c1 = c2 = c  
  
  #Extraindo x, dx/dt, theta, dtheta/dt do vetor X para facilitar o entendimento do codigo
  x = X[0]
  x_dot = X[1]
  theta = X[2]
  theta_dot = X[3]

  #Definindo omega em funcao da velocidade do carro e da distancia entre as "lombadas"
  w = 2*np.pi*V/L

  #Definindo as derivadas d(dx/dt)/dt e d(dtheta/dt)/dt
  x_dot_dot = (-k1*(x - a*theta - func_d1(t, w)) - k2*(x + b*theta - func_d2(t, w)) - c1*(x_dot - a*theta_dot - func_d1_dot(t, w)) - c2*(x_dot + b*theta_dot - func_d2_dot(t, w)) + Fn*np.sin(we*t))/M
  theta_dot_dot = (k1*(x - a*theta - func_d1(t, w))*a - k2*(x + b*theta - func_d2(t, w))*b + c1*(x_dot - a*theta_dot - func_d1_dot(t, w))*a - c2*(x_dot + b*theta_dot - func_d2_dot(t, w))*b - Fn*np.sin(we*t)*e - Fn*np.cos(we*t)*f)/Ic
  
  #Criando um novo vetor dX/dt que contem as derivadas em t de x, dx/dt, theta, dtheta/dt
  dX_dt = np.array([x_dot, x_dot_dot, theta_dot, theta_dot_dot])

  return dX_dt

#Função usada para plotar x e theta e suas derivadas
def plot_x_theta(t, X, dX_dt, h, ax1_title):
  x = X[:,0]
  x_dot = X[:,1]
  x_dot_dot = dX_dt[:,1]
  theta = X[:,2]
  theta_dot = X[:,3]
  theta_dot_dot = dX_dt[:,3]

  fig, (ax1, ax2) = plt.subplots(2,1)
  ax1.plot(t, x*10**2, label=r'$x\times 10^2 (m)$')                         #x é plotado multiplicado por 10^4 para aparecer no gráfico claramente
  ax1.plot(t, x_dot, label=r'$\dot{x}  (m/s)$')            #dx/dt é plotado multiplicado por 10^2 para aparecer no gráfico claramente
  ax1.plot(t, x_dot_dot*10**(-2), label=r'$\ddot{x} \times 10^{-2} (m/s^2)$')                       #d(dx/dt)/dt é plotado sem multiplicação
  legend = ax1.legend(loc='upper right', shadow=False, fontsize='x-large')
  ax1.set_title(ax1_title)

  ax2.plot(t, theta*10**2, label=r'$\theta \times 10^2 (rad)$')             #theta é plotado multiplicado por 10^4 para aparecer no gráfico claramente
  ax2.plot(t, theta_dot, label=r'$\dot{\theta}(rad/s)$')  #dtheta/dt é plotado multiplicado por 10^2 para aparecer no gráfico claramente
  ax2.plot(t, theta_dot_dot*10**(-2), label=r'$\ddot{\theta} \times 10^{-2} (rad/s^2)$')            #d(dtheta/dt)/dt é plotado sem multiplicação
  legend = ax2.legend(loc='upper right', shadow=False, fontsize='x-large')
  plt.show()

#Função que resolve o sistema com o Método de Runge-Kutta de 4a ordem
def rungeKutta(h:float, h_magnitude:str, V:float, k:float, c:float, t0=t0, tf=tf, X0=X_0):
  '''
  x_dot = dx/dt
  x_dot_dot = d(x_dot)/dt
  theta_dot = d(theta)/dt
  theta_dot_dot = d(theta_dot)/dt
  X = [xi, xi_dot, thetai, thetai_dot]
  dX_dt = [dx/dt, d(x_dot)/dt, d(theta)/dt, d(theta_dot_dot)/dt] = [f(t, x, theta, x1, x2, x3, x4)]
  '''  
  n = int((tf - t0)//h)       #n é o número de passos a ser dados entre o intervalo t0 e tf com passo h
  X = np.zeros((n+1,4))       #cria a matriz X para armazenar os resultados do vetor X em cada passo
  X[0] = X0                   #define o valor inicial do vetor X igual a X0 que contém os valores iniciais definidos anteriormente 
  dX_dt = np.zeros((n+1,4))   #cria a matriz dX/dt para armazenar os resultados do vetor dX/dt em cada passo

  t = np.arange(t0, tf, h)    #cria um array com os valores de cada tempo t em cada passo

  #Calcula K1, K2, K3, e K4 para cada passo i
  for i in range(n):
      K1 = func_dX_dt(t[i], X[i], V, k, c)
      K2 = func_dX_dt(t[i] + h/2, X[i] + h/2*K1, V, k, c)
      K3 = func_dX_dt(t[i] + h/2, X[i] + h/2*K2, V, k, c)
      K4 = func_dX_dt(t[i] + h, X[i] + h*K3, V, k, c)

      #Armazena os valores do vetor X para o passo i+1 calculado usando X[i] e os valores de K1 a K4
      X[i+1] = X[i] + h/6.0 * (K1 + 2 * K2 + 2 * K3 + K4)

      #Armazena os valores de dX_dt do passo i
      dX_dt[i] = K1
  #Plota os gráficos de x, dx/dt, d(dx/dt)/dt, theta, dtheta/dt, d(dtheta/dt)/dt em 0 <= t <= 4s
  plot_x_theta(t, X, dX_dt, h, h_magnitude)

  return X, dX_dt, t

