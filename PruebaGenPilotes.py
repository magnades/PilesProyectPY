import numpy as np
#
# LinkLim = np.array([126/3,126/3,126/3,246,390])*0.0254  # Reemplaza esto con tus datos reales
# LinkSpac= np.array([12,12,12,24,48])*0.0254
# MasaPilvc = np.array([2.246,3.124,4.002,4.88])*9.81  # Reemplaza esto con tus datos reales
# PileDiv = np.array([1, 1, 1, 1, 1])
# LongDes = 2.  # Reemplaza esto con tus datos reales
# Lz = 1.75  # Reemplaza esto con tus datos reales
# Apoyos_Pilvc = np.array([-17.,-19.,-21.,-23.])  # Reemplaza esto con tus datos reales
# Ly = np.array([1.5,2.3,2.3,3.05,3.05,3.05,3.05,3.05,3.05,2.3,2.3,1.5,1.5,1.5])  # Ejemplo de vector Ly
# NvTERR = -0.2457  # Ejemplo de NvTERR
# Pend = 2  # Ejemplo de Pend
#
# Fil_Masvc=np.array([4,6,8,10])

def generate_links(MasaPil, Fil_Masv, Apoyos_Pil, LinkSpac, LinkLim, Ly, NvTERR, Pend):
    num_pilotes = len(MasaPil)
    Link_V = np.zeros((25, num_pilotes))  #Arreglar el numero 25 que no se cuenta automáticamente del sistema
    NlinkV = np.zeros((len(LinkLim),num_pilotes))

    for i in range(num_pilotes):
        Link_V[0][i] = NvTERR - np.sum(Ly[:Fil_Masv[i]-1]) / Pend - 0.5 * LinkSpac[0]
        H = Link_V[0][i]
        # cont = 2
        # Limit = 1
        # contLk = 1

        cont = 1
        Limit = 0
        contLk = 1

        while H > Apoyos_Pil[i]:
            H = H - LinkSpac[Limit]
            if abs(H - Link_V[0][i] - 0.5 * LinkSpac[0]) - 1e-5 <= np.sum(LinkLim[:Limit+1]):
                Link_V[cont][i] = H
                cont += 1
                contLk += 1
            elif Limit < len(LinkLim)-1:
                H = Link_V[cont - 1][i] - 0.5 * (LinkSpac[Limit] + LinkSpac[Limit + 1])
                Link_V[cont][i] = H
                NlinkV[Limit][i] = contLk
                Limit += 1
                cont += 1
                # contLk = 1
                contLk = 1
            else:
                H = Link_V[cont - 1][i] - LinkSpac[Limit]
                Link_V[cont][i] = H
                cont += 1
                contLk += 1
        Link_V[cont - 1][i] = Apoyos_Pil[i]
        NlinkV[Limit][i] = contLk - 1
        # contLk = 0

    return Link_V, NlinkV

def generate_piles(MasaPil, Apoyos_Pil, Link_V, NlinkV, LinkSpac, Lz, LongDes, PileDiv):
    num_pilotes = len(MasaPil)

    Pile_V = np.zeros((Link_V.shape[0] + 2, Link_V.shape[1]))

    for i in range(num_pilotes):
        Pile_V[0][i] = Lz
        Pile_V[1][i] = Lz - LongDes
        cont = 2
        Limit = 0

        for k in range(int(np.sum(NlinkV[:, i]) - 1)):
            if Link_V[k + 1][i] != 0:
                if k <= np.sum(NlinkV[:Limit+1, i]) - 1:
                    pass
                elif Limit < len(LinkSpac)-1:
                    Limit += 1

                Space = np.linspace(Link_V[k][i], Link_V[k + 1][i], PileDiv[Limit] + 1)

                for l in range(int(len(Space) - 1)):
                    Pile_V[cont][i] = Space[l]
                    cont += 1

        Pile_V[cont][i] = Space[-1]
        Pile_V[cont + 1][i] = Apoyos_Pil[i]

    return Pile_V


# Link_Vc, NlinkVc = generate_links(MasaPilvc, Apoyos_Pilvc, LinkSpac, LinkLim, Ly, NvTERR, Pend)
# Pile_Vc = generate_piles(MasaPilvc, Apoyos_Pilvc, Link_Vc, NlinkVc, LinkSpac, Lz, LongDes, PileDiv)
#
# print(Pile_Vc)

# Puedes hacer lo mismo para los pilotes de VG con sus respectivos datos de entrada.
