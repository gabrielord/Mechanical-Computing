import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
plt.rcParams["figure.figsize"] = (20,12)

from Pontos import cria_Ponto

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
delta = L/8                # delta usado entre cada nó da malha para Phi
Lambda = 1.15                 # lambda usado na sobrerrelaxação de Phi
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

# iterando sobre os pontos para calcular o valor da temperatura
criteria = False
m = []
# for k in range(1):
while criteria == False:
  criteria = True

  # calcula a nova temperatura em cada ponto da malha e checa o critério da tolerância de convergência
  for linha in range(len(Points)):
    for coluna in range(len(Points[0])):
      # armazena a temperatura antiga do ponto
      old_temperature = Points[linha][coluna].temperature
      
      # calcula a nova temperatura no ponto (coluna, linha)
      Points[linha][coluna].temperature = Points[linha][coluna].Temperature(linha, coluna, Points, delta, Lambda, T_dentro, T_fora, T_motor, L, h, d, H, k_ar, rho_ar, c_p_ar)

      # armazena a nova temperatura do ponto
      new_temperature = Points[linha][coluna].temperature
      erro = 0

      # ignora pontos onde a nova temperatura é zero
      if new_temperature != 0:
        # calcula o erro entre a temperatura antiga no ponto e a nova temperatura calculada
        erro = abs((new_temperature - old_temperature)/new_temperature)
      
      # caso o erro seja maior que a tolerância, itera de novo
      if erro > convergence_tolerence:
        criteria = False

# coloca as temperatura em uma outra matriz separada
Temperature_matrix = np.zeros((len(Points), len(Points[0])))
for k in range(len(Points)):
  for p in range(len(Points[0])):
    Temperature_matrix[k][p] = Points[k][p].temperature

Heat_flux_on_top = Heat_flux_below = 0
nodes_on_top = nodes_below = 0
Heat_matrix = np.zeros((len(Points), len(Points[0])))
'''Segundo o site: https://www.icarros.com.br/volkswagen/fusca/1967/ficha-tecnica, a largura de um fusca 
é de, aproximadamente, 1.5 metros. Esse valor será usado para calcular a força resultando no carro.'''

# iterando somente sobre as linhas que coincidem com a posição do carro
for linha in range(int((H-h-L/2-delta)/delta),int((H-h)/delta+1)):
  l = []

  # iterando somente sobre as colunas que coincidem com a posição do carro
  for coluna in range(int(d/delta),int((d+L)/delta+1)):

    if Points[linha][coluna].dist_x != 1 and Points[linha][coluna].dist_y != 1:

      # calcula o fluxo de calor de pontos próximos do carro (até 0,01m do carro) na matriz de pressões sobre o carro
      try:
        theta = np.arctan((Points[linha][coluna].j*delta - h)/(Points[linha][coluna].i*delta - L/2 - d))
      except:
        theta = np.pi/2
        
      # ponto na esquerda do carro
      if Points[linha][coluna].dist_x >= 0:
        dTdx = (T_dentro - Points[linha][coluna-1].temperature)/(delta*(1+Points[linha][coluna].b))

      # a temperatura do contorno é T_motor
      elif (Points[linha][coluna-1].j * delta  <= Points[linha][coluna-1].i * delta * np.tan(np.pi/3) + h - np.tan(np.pi/3) * (L/2 + d)): 
        dTdx = (Points[linha][coluna+1].temperature - T_motor)/(delta*(1+Points[linha][coluna].b))

      # a temperatura do contorno é T_dentro
      else:
        dTdx = (Points[linha][coluna+1].temperature - T_dentro)/(delta*(1+Points[linha][coluna].b))

      # ponto embaixo do carro
      if Points[linha][coluna].dist_y >= 0:
        if (Points[linha][coluna].i*delta <= d + L/2): 
          # a temperatura do contorno é T_dentro
          dTdy = (T_dentro - Points[linha+1][coluna].temperature)/(delta*(1+Points[linha][coluna].a))
        else:
          # a temperatura do contorno é T_motor
          dTdy = (T_motor - Points[linha+1][coluna].temperature)/(delta*(1+Points[linha][coluna].a))          

        Heat_flux_below += -k_ar*( - dTdy)
        nodes_below += 1

      # a temperatura do contorno é T_motor
      elif (Points[linha+1][coluna].j * delta  <= Points[linha+1][coluna].i * delta * np.tan(np.pi/3) + h - np.tan(np.pi/3) * (L/2 + d)): 
        dTdy = (Points[linha-1][coluna].temperature - T_motor)/(delta*(1+Points[linha][coluna].a))
        Heat_flux_on_top += -k_ar*(np.cos(theta)*dTdx + np.sin(theta)*dTdy)
        nodes_on_top += 1

      # a temperatura do contorno é T_dentro
      else:
        dTdy = (Points[linha-1][coluna].temperature - T_dentro)/(delta*(1+Points[linha][coluna].a))
        Heat_flux_on_top += -k_ar*(np.cos(theta)*dTdx + np.sin(theta)*dTdy)
        nodes_on_top += 1
      
if nodes_below != 0:
  total_heat = (Heat_flux_on_top)*(np.pi*L/2/nodes_on_top*1.5) + (Heat_flux_below/nodes_below)*(L/nodes_below*1.5)
else:
  total_heat = (Heat_flux_on_top/nodes_on_top)*(np.pi*L/2/nodes_on_top*1.5)
  
print(f'A taxa de calor retirada do carro é: {round(total_heat,3)} W')

# gera gráfico das pressões na carroceria
fig, ax = plt.subplots(figsize=(12,8))
# norm = TwoSlopeNorm(vmin=0, vcenter=0, vmax=10)
ax=sns.heatmap(Temperature_matrix)
# ax=sns.heatmap(Pressure_on_car,cmap='coolwarm', norm=norm)
ax.collections[0].colorbar.set_label('Temperatura (K)', labelpad=20, fontsize=14)

props = dict(boxstyle='round', facecolor='white', alpha=0.75)

# place a text box in upper left in axes coords
ax.text(0.05, 0.95, f'Quantidade de calor trocada: {round(total_heat,3)} W', transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
plt.xticks([])
plt.yticks([])
plt.title('Temperatura dentro do túnel')
plt.show()