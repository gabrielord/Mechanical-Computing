'''Parte 2.3'''

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (20,12)

from Pontos import cria_Ponto, normalize

# Definindo as constantes do problema
L = 3                         # comprimento do carro [m]
d = 0.5*L                     # distancia do carro para a borda [m]
H = 2*L                       # altura da borda [m]
rho_ar = 1.25                 # densidade do ar [kg/m^3]
gamma_ar = 1.4                # coeficiente de expansão adiabática []
p_atm = 101325                # pressão atmosférica [Pa]
k_ar = 0.026                  # [W/(m.K)]
c_p_ar = 1002                 # calor específico do ar [J/(kg.K)]
T_dentro = 25 + 273                # temperatura dentro do carro [K]
T_fora = 20 + 273                  # temperatura fora do carro [K]
T_motor = 80 + 273                  # temperatura no motor [K]
delta = 0.05                # delta usado entre cada nó da malha para Phi
Lambda = 1.85                 # lambda usado na sobrerrelaxação de Phi
convergence_tolerence = 0.01  # tolerancia no erro da sobrerrelaxação

'''Parte 2.3 a)'''
V = 100/3.6                   # velocidade do vento [m/s]
h_list = [0.025, 0.05, 0.2] # lista com os valores de h utilizados
force_h_list = [] # lista com os valores da força calculados
delta = 0.05

for h in h_list:
    # criando malha com corrente nula
    Points = []
    # linhas (y)
    for linha in range(int(H/delta+1)):
      l = []
      # colunas (x)
      for coluna in range(int(H/delta+1)):
        l.append(cria_Ponto(linha, coluna, L, d, h, H, delta, T_fora, T_dentro, T_motor))
      Points.append(l)
    
    # iterando sobre os pontos para calcular o valor da corrente
    criteria = False
    m = []
    while criteria == False:
      criteria = True
    
      # calcula a nova corrente em cada ponto da malha e checa o critério da tolerância de convergência
      for linha in range(len(Points)):
        for coluna in range(len(Points[0])):
          # armazena a corrente antiga do ponto
          old_current = Points[linha][coluna].outflow_current
          
          # calcula a nova corrente no ponto (coluna, linha)
          Points[linha][coluna].outflow_current = Points[linha][coluna].new_outflow_current(linha, coluna, Points, delta, V, Lambda, L, d, H)
    
          # armazena a nova corrente do ponto
          new_current = Points[linha][coluna].outflow_current
          erro = 0
    
          # ignora pontos onde a nova corrente é zero
          if new_current != 0:
            # calcula o erro entre a corrente antiga no ponto e a nova corrente calculada
            erro = abs((new_current - old_current)/new_current)
          
          # caso o erro seja maior que a tolerância, itera de novo
          if erro > convergence_tolerence:
            criteria = False
    
    # Calcula o campo de velocidades em cada ponto
    # itera sobre as linhas
    for linha in range(len(Points)):
      # itera sobre as colunas
      for coluna in range(len(Points[0])):
        # calcula a velocidade no ponto (coluna, linha)
        Points[linha][coluna].velocity_field(linha, coluna, Points, delta, V, H, L, d, h)
    
    # Calcula a pressão em cada ponto
    # itera sobre as linhas
    for linha in range(len(Points)):
      # itera sobre as colunas
      for coluna in range(len(Points[0])):
        # calcula a pressão no ponto (coluna, linha)
        Points[linha][coluna].calc_pressure(V, rho_ar, gamma_ar)
        
    Pressure_on_top = 0
    Pressure_below = 0
    nodes_on_top = 0
    nodes_below = 0
    
    # iterando somente sobre as linhas que coincidem com a posicao do carro
    for linha in range(len(Points)):
      l = []
      # iterando somente sobre as colunas que coincidem com a posicao do carro
      for coluna in range(len(Points)):
        # calcula a força de pontos próximos do carro (até 0,01m do carro) na matriz de pressões sobre o carro
        if 0 < (Points[linha][coluna].i*delta - L/2 - d)**2 + (Points[linha][coluna].j*delta - h)**2 - (L/2)**2 < delta and Points[linha][coluna].j*delta>=h:
          try:
            theta = np.arctan((Points[linha][coluna].j*delta - h)/(abs(L/2 + d - Points[linha][coluna].i*delta)))
          except:
            theta = np.pi/2
          Pressure_on_top += Points[linha][coluna].pressure * np.sin(theta)
          nodes_on_top += 1
    
        elif 0 <= h - Points[linha][coluna].j*delta <= delta and d <= Points[linha][coluna].i*delta <= L+d:
          Pressure_below += Points[linha][coluna].pressure 
          nodes_below += 1
    
    Force_on_top = Pressure_on_top * (np.pi*L/2 / nodes_on_top) * 1.5
    
    Force_below = Pressure_below * L/nodes_below * 1.5
    # Força > 0 aponta na direção +y
    # Força < 0 aponta na direção -y
    Resultant_force = Force_below - Force_on_top
    
    force_h_list.append(Resultant_force)

    # printando resultados
    print(f'Força resultante para h = {h} m: {Resultant_force} N')

plt.plot(h_list, force_h_list)
plt.title('Variação da Força $lift$ com a distância h')
plt.xlabel('h (m)')
plt.ylabel('Força (N)')
plt.show()
    
'''Parte 2.3 b)'''

V_list = [75/3.6, 100/3.6, 140/3.6]                   # lista com os valores da velocidade do vento [m/s]
h = 0.15 # valor da altura do carro
force_V_list = [] # lista com os valores da força calculados

for v in V_list:
    # criando malha com corrente nula
    Points = []
    # linhas (y)
    for linha in range(int(H/delta+1)):
      l = []
      # colunas (x)
      for coluna in range(int(H/delta+1)):
        l.append(cria_Ponto(linha, coluna, L, d, h, H, delta, T_fora, T_dentro, T_motor))
      Points.append(l)
    
    # iterando sobre os pontos para calcular o valor da corrente
    criteria = False
    m = []
    while criteria == False:
      criteria = True
      # se for true, é preciso verificar a condição de erro na próxima iteração. se for false, não precisa
      flag = True
    
      # calcula a nova corrente em cada ponto da malha e checa o critério da tolerância de convergência
      for linha in range(len(Points)):
        for coluna in range(len(Points[0])):
          # armazena a corrente antiga do ponto
          old_current = Points[linha][coluna].outflow_current
          
          # calcula a nova corrente no ponto (coluna, linha)
          Points[linha][coluna].outflow_current = Points[linha][coluna].new_outflow_current(linha, coluna, Points, delta, v, Lambda, L, d, H)
    
          # armazena a nova corrente do ponto
          new_current = Points[linha][coluna].outflow_current
          erro = 0
    
          # ignora pontos onde a nova corrente é zero
          if new_current != 0:
            # calcula o erro entre a corrente antiga no ponto e a nova corrente calculada
            erro = abs((new_current - old_current)/new_current)
          
          # caso o erro seja maior que a tolerância, itera de novo
          if erro > convergence_tolerence:
            criteria = False
    
    # Calcula o campo de velocidades em cada ponto
    # itera sobre as linhas
    for linha in range(len(Points)):
      # itera sobre as colunas
      for coluna in range(len(Points[0])):
        # calcula a velocidade no ponto (coluna, linha)
        Points[linha][coluna].velocity_field(linha, coluna, Points, delta, v, H, L, d, h)
    

        # Calcula a pressão em cada ponto
    # itera sobre as linhas
    for linha in range(len(Points)):
      # itera sobre as colunas
      for coluna in range(len(Points[0])):
        # calcula a pressão no ponto (coluna, linha)
        Points[linha][coluna].calc_pressure(V, rho_ar, gamma_ar)
        
    Pressure_on_top = 0
    Pressure_below = 0
    nodes_on_top = 0
    nodes_below = 0
    
    # iterando somente sobre as linhas que coincidem com a posicao do carro
    for linha in range(len(Points)):
      l = []
      # iterando somente sobre as colunas que coincidem com a posicao do carro
      for coluna in range(len(Points)):
        # calcula a força de pontos próximos do carro (até 0,01m do carro) na matriz de pressões sobre o carro
        if 0 < (Points[linha][coluna].i*delta - L/2 - d)**2 + (Points[linha][coluna].j*delta - h)**2 - (L/2)**2 < delta and Points[linha][coluna].j*delta>=h:
          try:
            theta = np.arctan((Points[linha][coluna].j*delta - h)/(abs(L/2 + d - Points[linha][coluna].i*delta)))
          except:
            theta = np.pi/2
          Pressure_on_top += Points[linha][coluna].pressure * np.sin(theta)
          nodes_on_top += 1
    
        elif 0 <= h - Points[linha][coluna].j*delta <= delta and d <= Points[linha][coluna].i*delta <= L+d:
          Pressure_below += Points[linha][coluna].pressure 
          nodes_below += 1
    
    Force_on_top = Pressure_on_top * (np.pi*L/2 / nodes_on_top) * 1.5
    
    if nodes_below != 0:
      Force_below = Pressure_below * L/nodes_below * 1.5
    else:
      Force_below = 0
    # Força > 0 aponta na direção +y
    # Força < 0 aponta na direção -y
    Resultant_force = Force_below - Force_on_top
    
    force_V_list.append(Resultant_force)
    
    # printando resultados
    print(f'Força resultante para V = {v*3.6} km/h: {Resultant_force} N')

plt.plot(V_list, force_V_list)
plt.title('Variação da Força $lift$ com a velocidade')
plt.xlabel('V (km/h)')
plt.ylabel('Força (N)')
plt.xticks(ticks = np.array(V_list), labels=np.array(V_list)*3.6)
plt.show()