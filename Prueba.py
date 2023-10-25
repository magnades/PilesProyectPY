import numpy as np
from PruebaGenPilotes import generate_piles, generate_links, generar_elementos_pilote, calcular_coordenadas_placa


LinkLim = np.array([126/3,126/3,126/3,246,390])*0.0254  # Reemplaza esto con tus datos reales
LinkSpac= np.array([12,12,12,24,48])*0.0254
MasaPilvc = np.array([2.246,3.124,4.002,4.88])*9.81     # Masa en los nodos desde A,B,....F Vigas cabezales
MasaPilvg = np.array([1.368,5.758])*9.81                # Masa en los nodos desde A,B,....F Vigas Grua
PileDiv = np.array([1, 1, 1, 1, 1])                     # Número de elementos Pile entre cada tipo de link
LongDes = 2.  # Longitud de desarrollo de la barra longitudinal del pilote en la conexión
# Lz = 1.75  # Nivel del eje de la placa en coordenadas globales
Apoyos_Pilvc = np.array([-17.,-19.,-21.,-23.])          # Cota de apoyos para los pilotes en el eje de vigas cabezales
Apoyos_Pilvg = np.array([-15.,-26.])                    # Cota de apoyos para los pilotes en el eje de vigas grua
# Ly = np.array([1.5,2.3,2.3,3.05,3.05,3.05,3.05,3.05,3.05,2.3,2.3,1.5,1.5,1.5])  # Ejemplo de vector Ly
NvTERR = -0.2457  # Nivel del terreno de la primera fila en coordenadas globales
Pend = 2  # Pendiente del terreno
Fil_Masvc = np.array([4,6,8,10])        # Filas de vigas cabezales
Fil_Masvg = np.array([1,13])            # Filas de vigas de grua

OrdenPile=np.array([2,1,1,1,1,2])       #Se coloca en orden desde A,B... que tipo de pilote es (VC o VG)

N=81    #Numero de porticos
# Para Col_Masvc
Col_Masvc = np.arange(1, (N - 1) * 3 + 2, 3)
# Para Col_Masvg
Col_Masvg = np.arange(1, (N - 1) * 3 + 2)

ctr1=1  #(0=EJE 1, 1=TODOS EJES) Se crean todos los ejes de pilotes o solo el primer eje?

Ly=np.array([1.5,2.3,2.3,3.05,3.05,3.05,3.05,3.05,3.05,2.3,2.3,1.5,1.5,1.5])     #Vector de espaciamientos en el sentido Y
Lx=2.5                                                                           #Separación típica en x
Lz=1.75                                                                          #Altura donde esta la placa en coordenadas globales.


Link_Vc, NlinkVc = generate_links(MasaPilvc, Fil_Masvc, Apoyos_Pilvc, LinkSpac, LinkLim, Ly, NvTERR, Pend)
Pile_Vc = generate_piles(MasaPilvc, Apoyos_Pilvc, Link_Vc, NlinkVc, LinkSpac, Lz, LongDes, PileDiv)

Link_Vg, NlinkVg = generate_links(MasaPilvg, Fil_Masvg, Apoyos_Pilvg, LinkSpac, LinkLim, Ly, NvTERR, Pend)
Pile_Vg = generate_piles(MasaPilvg, Apoyos_Pilvg, Link_Vg, NlinkVg, LinkSpac, Lz, LongDes, PileDiv)

Coord = calcular_coordenadas_placa(Lx, Ly, Lz, N)
# print(Coord['x_coords'][0])
resultado = generar_elementos_pilote(ctr1, MasaPilvc, MasaPilvg, Col_Masvc, Col_Masvg, Coord, Fil_Masvc, Fil_Masvg, Pile_Vc, Pile_Vg, OrdenPile)
print(resultado)