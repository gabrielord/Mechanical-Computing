import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import pandas as pd
plt.rcParams["figure.figsize"] = (20,12)

# classe do objeto Ponto que armazenará as informações deste ponto (corrente, temperatura, etc)
class cria_Ponto:
  def __init__(self, linha, coluna, L, d, h, H, delta, T_fora, T_dentro, T_motor):
    # índice referente à posição x: x = self.i * delta
    self.i = coluna
    # índice referente à posição y: y =  self.j * delta
    self.j = round(H/delta - linha)
    # corrente de ar inicial em cada nó
    self.outflow_current = 75
    # temperatura inicial em cada nó
    self.temperature = 0
    # temperatura na borda esquerda
    if self.i == 0:
      self.temperature = T_fora
    
    
    # Para definir a distância do ponto até o carro, temos os if e else abaixo:
    # ponto na mesma posição x do carro
    if d <= self.i*delta <= (d + L):
      # ponto embaixo do carro e com distância até ele menor que delta
      if 0 < h - delta*self.j < delta:
        dist_y = h - delta*self.j

      # ponto em cima do carro ou dentro do carro
      else:
        dist_y = (h + np.sqrt((L/2)**2 - (self.i*delta - L/2 - d)**2)) - self.j*delta 

      # distância do ponto ao carro maior que o delta
      if abs(dist_y) >= delta:
        dist_y = delta

    # ponto longe do carro (não tem coordenada x igual ao carro)
    else:
      dist_y = delta

    # ponto na mesma posição y do carro
    if h <= self.j*delta <= (L/2 + h):

      # ponto no lado esquerdo do carro
      if self.i*delta < (d + L/2):
        dist_x = (L/2 + d - np.sqrt((L/2)**2 - (self.j*delta - h)**2)) - self.i*delta

      # ponto no lado direito do carro ou no meio dele 
      else:
        dist_x = (L/2 + d + np.sqrt((L/2)**2 - (self.j*delta - h)**2)) - self.i*delta

      # distância do ponto ao carro maior que o delta
      if abs(dist_x) >= delta:
        dist_x = delta

    # ponto longe do carro (não tem coordenada y igual ao carro)
    else:
      dist_x = delta

    # ponto dentro ou na borda do carro
    if ((self.i*delta - L/2 - d)**2 + (self.j*delta - h)**2 <= (L/2)**2 and self.j*delta >= h) or (np.isclose(dist_x,0,5) and np.isclose(dist_y,0,5)):
      # definimos as distancias ao carro como zero (pois estão dentro do carro)
      dist_y = dist_x = 0

      # dentro do carro não há corrente de ar
      self.outflow_current = 0

      # temperatura dentro do carro
      if self.j*delta <= np.tan(np.pi/3)*self.i*delta + h - np.tan(np.pi/3)*(L/2+d):
        self.temperature = T_motor
      else:
        self.temperature = T_dentro

    # na borda inferior não há corrente de ar
    if self.j == 0:
      self.outflow_current = 0

    # distância do ponto ao carro
    self.dist_y = dist_y
    # fator usado nos contornos irregulares: a*delta = distancia em x do ponto até o carro
    self.a = abs(dist_y / delta)

    # distância do ponto ao carro
    self.dist_x = dist_x
    # fator usado nos contornos irregulares: b*delta = distancia em y do ponto até o carro
    self.b = abs(dist_x / delta)

  # função usada para conseguir uma nova iteração da corrente, ou seja, calcula a nova corrente no ponto (coluna, linha)
  def new_outflow_current(self, linha, coluna, Points, delta, V, Lambda, L, d, H):

    # ponto dentro ou na borda do carro
    if self.dist_x == self.dist_y == 0:
      # condição de contorno: corrente dentro do carro = 0
      psi_i_j = 0

    # ponto no interior da malha: fora das bordas e fora do carro
    elif 0 < self.i*delta < (L + 2*d) and 0 < self.j*delta < H:

      # ponto no lado esquerdo e perto do carro, horizontalmente
      if 0 < self.dist_x < delta:
        # armazena psi(i+1,j) e psi(i-1,j), respectivamente
        psi_x_1, psi_x_2 = 0, Points[linha][coluna-1].outflow_current

      # ponto no lado direito e perto do carro, horizontalmente
      elif -1*delta < self.dist_x < 0:
        # armazena psi(i-1,j) e psi(i+1,j), respectivamente
        psi_x_1, psi_x_2 = 0, Points[linha][coluna+1].outflow_current

      # ponto longe do carro, horizontalmente
      else:
        # armazena psi(i-1,j) e psi(i+1,j), respectivamente
        psi_x_1, psi_x_2 = Points[linha][coluna-1].outflow_current, Points[linha][coluna+1].outflow_current

      # ponto embaixo e perto do carro, verticalmente
      if 0 < self.dist_y < delta:
        # armazena psi(i,j+1) e psi(i,j-1), respectivamente
        psi_y_1, psi_y_2 = 0, Points[linha+1][coluna].outflow_current

      # ponto em cima e perto do carro, verticalmente
      elif -1*delta < self.dist_y < 0:
        # armazena psi(i,j-1) e psi(i,j+1), respectivamente
        psi_y_1, psi_y_2 = 0, Points[linha-1][coluna].outflow_current

      # ponto longe do carro, verticalmente
      else:
        # armazena psi(i,j+1) e psi(i,j-1), respectivamente
        psi_y_1, psi_y_2 = Points[linha-1][coluna].outflow_current, Points[linha+1][coluna].outflow_current
        
      # usa os valores de psi próximos a psi(i,j) e os parametros 'a' e 'b' para calcular a nova corrente psi(i,j)
      psi_i_j = 1/(self.a + self.b) * (self.b/(1+self.a) * (psi_y_1 + self.a*psi_y_2) + self.a/(1+self.b) * (psi_x_1 + self.b*psi_x_2))

    # ponto no canto superior esquerdo 
    elif self.i == 0 and self.j*delta == H:
      # psi(i,j) calculado usando as condições de contorno em (0,H): dpsi/dx = 0 e dpsi/dy = V
      psi_i_j = 1/2*(Points[linha][coluna+1].outflow_current + Points[linha+1][coluna].outflow_current + delta*V)

    # ponto no canto superior direito
    elif self.i*delta == self.j*delta == H:
      # psi(i,j) calculado usando as condições de contorno em (H,H): dpsi/dx = 0 e dpsi/dy = V
      psi_i_j = 1/2*(Points[linha][coluna-1].outflow_current + Points[linha+1][coluna].outflow_current + delta*V)

    # ponto na borda inferior
    elif self.j == 0:
      # condição de contorno: na borda inferior corrente = 0
      psi_i_j = 0

    # ponto na borda superior
    elif self.j * delta == H:
      # psi(i,j) calculado usando as condições de contorno em (i,H): dpsi/dy = V
      psi_i_j = 1/4*(Points[linha][coluna+1].outflow_current + Points[linha][coluna-1].outflow_current + 2*Points[linha+1][coluna].outflow_current + 2*delta*V)

    # ponto na borda esquerda
    elif self.i == 0:
      # psi(i,j) calculado usando as condições de contorno em (0,j): dpsi/dx = 0
      psi_i_j = 1/4*(2*Points[linha][coluna+1].outflow_current + Points[linha-1][coluna].outflow_current + Points[linha+1][coluna].outflow_current)

    # ponto na borda direita
    elif self.i * delta == H:
      # psi(i,j) calculado usando as condições de contorno em (H,j): dpsi/dx = 0
      psi_i_j = 1/4*(2*Points[linha][coluna-1].outflow_current + Points[linha-1][coluna].outflow_current + Points[linha+1][coluna].outflow_current)

    # se não estiver em nenhum caso, printar mensagem de erro
    else:
      print("Erro")

    # retorna o novo valor de psi(i,j) calculado usando sobrerrelaxação
    return Lambda * psi_i_j + (1-Lambda) * self.outflow_current

  # função que calcula o campo de velocidade no ponto (coluna, linha)
  def velocity_field(self, linha, coluna, Points, delta, V, H, L, d, h):
    '''
    velocidade = ui + vj
    d(Phi)/dx = -v
    d(Phi)/dy = u
    '''
  
    # ponto dentro ou na borda do carro
    if self.a == self.b == 0:
      # condições de contorno: dentro do carro o campo de velocidades = 0
      self.u = 0
      self.v = 0

    # ponto no canto superior esquerdo ou direito
    elif (self.i == 0 and self.j * delta == H) or (self.i * delta == self.j * delta == H):
      # condições de contorno: nos cantos superiores da malha dpsi/dy = V e dpsi/dx = 0
      self.u = V 
      self.v = 0

    # ponto na borda inferior
    elif self.j == 0:
      # condições de contorno: na borda inferior o campo de velocidades = 0
      self.u = 0
      self.v = 0

    # ponto na borda superior
    elif self.j * delta == H:
      # condição de contorno: na borda superior dpsi/dy = V
      self.u = V 
      # primeira diferença progressiva em x
      self.v = -1*(Points[linha][coluna+1].outflow_current - Points[linha][coluna-1].outflow_current)/(2*delta)

    # ponto na borda esquerda ou direita
    elif self.i == 0 or self.i * delta == H:
      # primeira diferença progressiva em y
      self.u = (Points[linha-1][coluna].outflow_current - Points[linha+1][coluna].outflow_current)/(2*delta)
      # condição de contorno: na borda esquerda ou direita dpsi/dx = 0
      self.v = 0
    
    # ponto no interior do túnel
    else:
        """Para os casos em que o nó está dentro do túnel, vamos calcular primeiro
        a velocidade em y e, depois, em x."""

        # ponto longe do carro horizontalmente
        if self.dist_x == delta:
          # primeira diferença central em x
          self.v = -1*(Points[linha][coluna+1].outflow_current - Points[linha][coluna-1].outflow_current)/(2*delta)

        # ponto perto do carro pela esquerda
        elif 0 < self.dist_x < delta:
          # primeira diferença central em x com contorno irregular
          self.v = -1*(0 - Points[linha][coluna-1].outflow_current)/(delta + abs(self.dist_x))

        # ponto perto do carro pela direita
        else:
          # primeira diferença central em x com contorno irregular
          self.v = -1*(Points[linha][coluna+1].outflow_current - 0)/(delta + abs(self.dist_x))

        # ponto longe do carro verticalmente
        if self.dist_y == delta:
          # primeira diferença central em y
          self.u = (Points[linha-1][coluna].outflow_current - Points[linha+1][coluna].outflow_current)/(2*delta)

        # ponto perto do carro por baixo
        elif 0 < self.dist_y < delta:
          # primeira diferença central em y com contorno irregular
          self.u = (0 - Points[linha+1][coluna].outflow_current)/(delta + abs(self.dist_y))

        # ponto perto do carro por cima
        else:
          # primeira diferença central em y com contorno irregular
          self.u = (Points[linha-1][coluna].outflow_current - 0)/(delta + abs(self.dist_y))
          
  # função que calcula a pressão no ponto (coluna, linha) usando como base o campo de velocidades naquele ponto
  def calc_pressure(self, V, rho_ar, gamma_ar):
    # calcula a pressão
    self.pressure = rho_ar * (gamma_ar - 1)/gamma_ar * (V**2 - np.sqrt(self.u**2+self.v**2)**2)/2
  
  # função que calcula a temperatura em cada ponto da malha
  def Temperature(self, linha, coluna, Points, delta, Lambda, T_dentro, T_fora, T_motor, L, h, d, H, k_ar, rho_ar, c_p_ar):
    # dentro ou na borda do carro
    if self.a == self.b == 0:

      # dentro da zona do motor
      if (self.j * delta  <= self.i * delta * np.tan(np.pi/3) + h - np.tan(np.pi/3) * (L/2 + d)): 
        Tij = T_motor

      # fora da zona do motor
      else:
        Tij = T_dentro

    # canto superior direito
    elif (self.i * delta == self.j * delta == H):
      # condições de contorno: dT/dx = 0 e dT/dy = 0
      Tij = (Points[linha+1][coluna].temperature + Points[linha][coluna-1].temperature)/2
    
    # canto inferior direito
    elif (self.i * delta == H) and (self.j * delta == 0):
      # condições de contorno: dT/dx = 0 e dT/dy = 0
      Tij = (Points[linha-1][coluna].temperature + Points[linha][coluna-1].temperature)/2

    # borda esquerda
    elif self.i == 0:
      # condição de contorno
      Tij = T_fora

    # borda direita
    elif self.i * delta == H:
      # condição de contorno: dT/dx = 0
      Tij = 1/4*(Points[linha-1][coluna].temperature + Points[linha+1][coluna].temperature + 2*Points[linha][coluna-1].temperature)
      
    # borda inferior
    elif self.j == 0:
      # condição de contorno: dT/dy = 0
      Tij = k_ar/(4*k_ar + rho_ar*c_p_ar*self.u*delta)*(2*Points[linha-1][coluna].temperature + Points[linha][coluna+1].temperature + Points[linha][coluna-1].temperature*(1 + delta*rho_ar*c_p_ar*self.u/k_ar))

    # borda superior
    elif self.j * delta == H:
      # condição de contorno: dT/dy = 0
      Tij = k_ar/(4*k_ar + rho_ar*c_p_ar*self.u*delta)*(2*Points[linha+1][coluna].temperature + Points[linha][coluna+1].temperature + Points[linha][coluna-1].temperature*(1 + delta*rho_ar*c_p_ar*self.u/k_ar))

    # interior do túnel
    elif 0 < self.i*delta < (L + 2*d) and 0 < self.j*delta < H:
        # longe do carro horizontalmente
        if self.dist_x == delta:
          # valores usados na segunda derivada em x
          Tx1, Tx2 = Points[linha][coluna+1].temperature, Points[linha][coluna-1].temperature
          # valor usado na primeira derivada em x
          Tx3 = Points[linha][coluna-1].temperature

        # perto do carro pela esquerda
        elif 0 < self.dist_x < delta:
          # valores usados na segunda derivada em x
          Tx1, Tx2 = T_dentro, Points[linha][coluna-1].temperature
          # valor usado na primeira derivada em x
          Tx3 = Points[linha][coluna-1].temperature

        # perto do carro pela direita
        else:
          # a temperatura do contorno é T_motor
          if (Points[linha][coluna-1].j * delta  <= Points[linha][coluna-1].i * delta * np.tan(np.pi/3) + h - np.tan(np.pi/3) * (L/2 + d)): 
            # valores usados na segunda derivada em x 
            Tx1, Tx2 = T_motor, Points[linha][coluna+1].temperature
            # valor usado na primeira derivada em x
            Tx3 = T_motor

          # a temperatura do contorno é T_dentro
          else:
            # valores usados na segunda derivada em x
            Tx1, Tx2 = T_dentro, Points[linha][coluna+1].temperature
            # valore usado na primeira derivada em x
            Tx3 = T_dentro

        # longe do carro verticalmente
        if self.dist_y == delta:
          # valores usados na segunda derivada em y
          Ty1, Ty2 = Points[linha+1][coluna].temperature, Points[linha-1][coluna].temperature

          # velocidade vj > 0
          if self.v>=0:
            # valor usado na primeira derivada em y
            Ty3 = Points[linha+1][coluna].temperature
          
          # velocidade vj < 0
          else:
            # valor usado na primeira derivada em y
            Ty3 = Points[linha-1][coluna].temperature

        # perto do carro por baixo
        elif 0 < self.dist_y < delta:

          # a temperatura do contorno é T_dentro
          if (self.i*delta <= d + L/2): 
            # valores usados na segunda derivada em y
            Ty1, Ty2 = T_dentro, Points[linha+1][coluna].temperature
            
            # velocidade vj > 0
            if self.v>=0:
              # valor usado na primeira derivada em y
              Ty3 = Points[linha+1][coluna].temperature
            else:
              # valor usado na primeira derivada em y
              Ty3 = T_dentro

          # a temperatura do contorno é T_motor
          else:
            # valores usados na segunda derivada em y
            Ty1, Ty2 = T_motor, Points[linha+1][coluna].temperature
            
            # velocidade vj > 0
            if self.v>=0:
              # valor usado na primeira derivada em y
              Ty3 = Points[linha+1][coluna].temperature
            else:
              # valor usado na primeira derivada em y
              Ty3 = T_motor

        # perto do carro por cima
        else:
          # a temperatura do contorno é T_motor
          if (Points[linha+1][coluna].j * delta  <= Points[linha+1][coluna].i * delta * np.tan(np.pi/3) + h - np.tan(np.pi/3) * (L/2 + d)): 
            # valores usados na segunda derivada em y
            Ty1, Ty2 = T_motor, Points[linha-1][coluna].temperature

            # velocidade vj > 0
            if self.v>=0:
              # valor usado na primeira derivada em y
              Ty3 = T_motor
            else:
              # valor usado na primeira derivada em y
              Ty3 = Points[linha-1][coluna].temperature

          else:
            # valores usados na segunda derivada em y
            Ty1, Ty2 = T_dentro, Points[linha-1][coluna].temperature

            # velocidade vj > 0
            if self.v>=0:
              # valor usado na primeira derivada em y
              Ty3 = T_dentro
            # velocidade vj < 0  
            else:
              # valor usado na primeira derivada em y
              Ty3 = Points[linha-1][coluna].temperature

        if self.v >= 0:
          # expressão para v > 0
          Tij = 2*k_ar/(2*k_ar/delta*(1/self.a + 1/self.b) + rho_ar*c_p_ar*(self.u/self.b + self.v/self.a)) * ((Tx1+self.b*Tx2)/(delta*self.b*(self.b+1)) \
                                                                    + (Ty1+self.a*Ty2)/(delta*self.a*(self.a+1)) + rho_ar*c_p_ar/(2*k_ar)*(Tx3*self.u/self.b + Ty3*self.v/self.a))                                                                                       

        else:
          # expressão para v < 0
          Tij = 2*k_ar/(2*k_ar/delta*(1/self.a + 1/self.b) + rho_ar*c_p_ar*(self.u/self.b - self.v/self.a)) * ((Tx1+self.b*Tx2)/(delta*self.b*(self.b+1)) \
                                                                    + (Ty1+self.a*Ty2)/(delta*self.a*(self.a+1)) + rho_ar*c_p_ar/(2*k_ar)*(Tx3*self.u/self.b - Ty3*self.v/self.a))


    # se não estiver em nenhum caso, printar mensagem de erro
    else:
      print("Erro")

    # retorna o novo valor de psi(i,j) calculado usando sobrerrelaxação
    return Lambda*Tij + (1-Lambda)*self.temperature

# função para normalizar vetores
def normalize(v1, v2):
    norm = np.linalg.norm([v1,v2])
    if norm == 0: 
       return v1, v2, norm
    return v1 / norm, v2 / norm, norm