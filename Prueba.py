import numpy as np
from PruebaGenPilotes import generate_piles, generate_links

LinkLim = np.array([126/3,126/3,126/3,246,390])*0.0254  # Reemplaza esto con tus datos reales
LinkSpac= np.array([12,12,12,24,48])*0.0254
MasaPilvc = np.array([2.246,3.124,4.002,4.88])*9.81  # Reemplaza esto con tus datos reales
PileDiv = np.array([1, 1, 1, 1, 1])
LongDes = 2.  # Reemplaza esto con tus datos reales
Lz = 1.75  # Reemplaza esto con tus datos reales
Apoyos_Pilvc = np.array([-17.,-19.,-21.,-23.])  # Reemplaza esto con tus datos reales
Ly = np.array([1.5,2.3,2.3,3.05,3.05,3.05,3.05,3.05,3.05,2.3,2.3,1.5,1.5,1.5])  # Ejemplo de vector Ly
NvTERR = -0.2457  # Ejemplo de NvTERR
Pend = 2  # Ejemplo de Pend
Fil_Masvc=np.array([4,6,8,10])


Link_Vc, NlinkVc = generate_links(MasaPilvc, Fil_Masvc, Apoyos_Pilvc, LinkSpac, LinkLim, Ly, NvTERR, Pend)
Pile_Vc = generate_piles(MasaPilvc, Apoyos_Pilvc, Link_Vc, NlinkVc, LinkSpac, Lz, LongDes, PileDiv)

print(Pile_Vc)