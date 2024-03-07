import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
plt.rcParams["figure.figsize"] = (20,12)

from Pontos import cria_Ponto, normalize

# Definindo as constantes do problema
V = 100/3.6                   # velocidade do vento [m/s]
h = 0.15                       # altura abaixo do carro [m]
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
delta = 0.25*h                # delta usado entre cada nó da malha para Phi
Lambda = 1.85                 # lambda usado na sobrerrelaxação de Phi
convergence_tolerence = 0.01  # tolerancia no erro da sobrerrelaxação

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

# coloca as correntes em uma outra matriz separada
Current = np.zeros((len(Points), len(Points[0])))
for k in range(len(Points)):
  for p in range(len(Points[0])):
    Current[k][p] = Points[k][p].outflow_current

# gera o gráfico das correntes de ar
fig, ax = plt.subplots(figsize=(12,8))

CS = ax.contour(Current, origin='upper', levels=np.arange(0,150,5), colors='white')
ax.clabel(CS, inline=True, fontsize=10)
ax.contourf(Current, origin='upper')
plt.xticks(ticks=plt.xticks()[0][:-1], labels=delta * np.array(plt.xticks()[0][:-1], dtype=np.float64))
plt.yticks(ticks=plt.yticks()[0][:-1], labels=delta * np.array(plt.yticks()[0][:-1], dtype=np.float64))
plt.title("Linhas de corrente de ar")

plt.show()

# Calcula o campo de velocidades em cada ponto
# itera sobre as linhas
for linha in range(len(Points)):
  # itera sobre as colunas
  for coluna in range(len(Points[0])):
    # calcula a velocidade no ponto (coluna, linha)
    Points[linha][coluna].velocity_field(linha, coluna, Points, delta, V, H, L, d, h)

# encontra o vetor normalizado das velocidades
U_norm = []
V_norm = []
Magnitude = []
for linha in range(0,len(Points),5):
  l = []
  p = []
  k = []
  for coluna in range(0,len(Points[0]),5):
    # print(Points[linha][coluna].u - Points[linha][coluna].v)
    u_norm, v_norm, norm = normalize(Points[linha][coluna].u, Points[linha][coluna].v)
    l.append(u_norm)
    p.append(v_norm)
    k.append(norm)
  U_norm.append(l)
  V_norm.append(p)
  Magnitude.append(k)

# gerando o gráfico das velocidades
fig, ax = plt.subplots(figsize=(12,8))

plt.quiver(U_norm, V_norm, Magnitude, cmap='viridis')
ax.invert_yaxis()
# Adicionar uma barra de cores
plt.colorbar(label='Velocidade (m/s)')
plt.xticks([])
plt.yticks([])
plt.title('Campo de Velocidades')
plt.show()

## Parte I c)

# Calcula a pressão em cada ponto
# itera sobre as linhas
for linha in range(len(Points)):
  # itera sobre as colunas
  for coluna in range(len(Points[0])):
    # calcula a pressão no ponto (coluna, linha)
    Points[linha][coluna].calc_pressure(V, rho_ar, gamma_ar)

# coloca as pressões em uma nova matriz separada
Pressure = np.zeros((len(Points), len(Points[0])))
for linha in range(len(Points)):
  for coluna in range(len(Points[0])):
    Pressure[linha][coluna] = Points[linha][coluna].pressure

# gera o gráfico das correntes de ar
fig, ax = plt.subplots(figsize=(12,8))

cf = ax.contourf(Pressure, origin='upper', cbar=True)
plt.xticks(ticks=plt.xticks()[0][:-1], labels=delta * np.array(plt.xticks()[0][:-1], dtype=np.float64))
plt.yticks(ticks=plt.yticks()[0][:-1], labels=delta * np.array(plt.yticks()[0][:-1], dtype=np.float64))
fig.colorbar(cf)
plt.title("Pressão relativa do ar [Pa]")
plt.show()

## Parte I d)

# coloca as pressões sobre o carro em uma matriz separada
Pressure_per_distance_above = []
position_x = []
min_pressure = 150

# iterando somente sobre as colunas que coincidem com a posição do carro
for coluna in range(int(d/delta),int((d+L)/delta+1)):

  # iterando somente sobre as linhas que coincidem com a posição do carro
  for linha in range(int((H-h)/delta+1),int((H-h-L/2-delta)/delta),-1):
    l = []

    # adiciona as pressões de pontos próximos do carro (até 0,01m do carro) na matriz de pressões sobre o carro e ignora os pontos que fazem parte da carroceria
    if 0 < (Points[linha][coluna].i*delta - L/2 - d)**2 + (Points[linha][coluna].j*delta - h)**2 - (L/2)**2 < 5e-2 and Points[linha][coluna].j*delta>=h and Points[linha][coluna].a!=0:
        Pressure_per_distance_above.append(Points[linha][coluna].pressure)
        position_x.append(Points[linha][coluna].i*delta)

        if Points[linha][coluna].pressure < min_pressure:
          min_pressure = Points[linha][coluna].pressure
          min_coluna, min_linha = coluna, linha

# gera o gráfico da pressão ao longo da parte superior do carro
fig, ax = plt.subplots(figsize=(12,8))
ax.plot(position_x, Pressure_per_distance_above)
ax.annotate(f'Mínima Pressão Relativa: {round(min_pressure,2)} Pa', xy=(round(min_coluna*delta), min_pressure), xytext=((min_coluna-12)*delta,min_pressure+70),
            arrowprops=dict(facecolor='black', shrink=0.05))
plt.xlabel('Coordenada x do nó')
plt.ylabel('Pressão relativa do ar [Pa]')
ax.set_title("Pressão relativa do ar no superior do carro [Pa]")
plt.show()



## Parte I e)

Pressure_on_top = 0
Pressure_below = 0
nodes_on_top = 0
nodes_below = 0

'''Segundo o site: https://www.icarros.com.br/volkswagen/fusca/1967/ficha-tecnica, a largura de um fusca 
é de, aproximadamente, 1.5 metros. Esse valor será usado para calcular a força resultando no carro.'''

# coloca as pressoes sobre o carro em uma matriz separada
min_pressure = 150 # valor inicial da pressao minima (ja que a maxima pressao no mapa < 150)
Pressure_on_car = np.zeros((len(Points),len(Points)))

# iterando somente sobre as linhas que coincidem com a posicao do carro
for linha in range(int((H-h-L/2-delta)/delta),int((H-h)/delta+1)):
  l = []
  # iterando somente sobre as colunas que coincidem com a posicao do carro
  for coluna in range(int(d/delta),int((d+L)/delta+1)):
    # calcula a força de pontos próximos do carro na matriz de pressões sobre o carro
    if 0 < (Points[linha][coluna].i*delta - L/2 - d)**2 + (Points[linha][coluna].j*delta - h)**2 - (L/2)**2 < delta and Points[linha][coluna].j*delta>=h:
      try:
        theta = np.arctan((Points[linha][coluna].j*delta - h)/(abs(L/2 + d - Points[linha][coluna].i*delta)))
      except:
        theta = np.pi/2
      Pressure_on_top += Points[linha][coluna].pressure * np.sin(theta)
      nodes_on_top += 1

      Pressure_on_car[linha][coluna] = Points[linha][coluna].pressure
      # se a pressao no ponto (coluna, linha) for menor do que a pressao minima atual, a nova pressao sera a pressao minima
      if Points[linha][coluna].pressure < min_pressure:
        min_pressure = Points[linha][coluna].pressure

        # armazena a linha e a coluna do ponto de pressao minima
        min_linha = linha
        min_coluna = coluna

    # parte inferior do carro
    elif Points[linha][coluna].j*delta == h and d <= Points[linha][coluna].i*delta <= L+d:
      Pressure_below += Points[linha][coluna].pressure 
      Pressure_on_car[linha][coluna] = Points[linha][coluna].pressure
      nodes_below += 1
        
    
      # se a pressao no ponto (coluna, linha) for menor do que a pressao minima atual, a nova pressao sera a pressao minima
      if Points[linha][coluna].pressure < min_pressure:
        min_pressure = Points[linha][coluna].pressure

        # armazena a linha e a coluna do ponto de pressao minima
        min_linha = linha
        min_coluna = coluna
    
    
  

Force_on_top = Pressure_on_top * (np.pi*L/2 / nodes_on_top) * 1.5
Force_below = Pressure_below * L/nodes_below * 1.5

# como a forca em cima e negativa, ela aponta na direcao contraria da pressao
Resultant_force = Force_below - Force_on_top

# printando resultados
print(f' Força em cima: {Force_on_top} N \n Força embaixo: {Force_below} N \n Força resultante: {Force_below - Force_on_top} N')


# gera gráfico das pressões na carroceria
fig, ax = plt.subplots(figsize=(12,8))
norm = TwoSlopeNorm(vmin=Pressure_on_car.min(), vcenter=0, vmax=Pressure_on_car.max())
# ax=sns.heatmap(Pressure_on_car,cmap='coolwarm')
ax=sns.heatmap(Pressure_on_car,cmap='coolwarm', norm=norm)
ax.collections[0].colorbar.set_label('Pressão (Pa)', labelpad=20, fontsize=14)
ax.annotate(f'Mínima Pressão Relativa: {round(min_pressure,2)} Pa', xy=(min_coluna,min_linha), xytext=(H/(2*delta)-20,(H-h-L/2)/delta-20),
            arrowprops=dict(facecolor='black', shrink=0.05))
plt.xticks([])
plt.yticks([])
plt.title('Pressão relativa à pressão atmosférica na carroceria do carro [Pa]')
plt.show()