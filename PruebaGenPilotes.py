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
    Link_V = np.zeros((26, num_pilotes))  #Arreglar el numero 25 que no se cuenta automáticamente del sistema
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

# def generar_elementos_pilote(ctr1, MasaPilvc, MasaPilvg, Col_Masvc, Col_Masvg, Coord, Fil_Masvc, Fil_Masvg, Pile_Vc, Pile_Vg, OrdenPile):
#     if ctr1 == 0:
#         # m = 1
#         # n = 1
#         # o = 1
#         m = 0
#         n = 0
#         o = 0
#     elif ctr1 == 1:
#         # m = len(Col_Masvc)
#         # n = 1
#         # o = len(Col_Masvg)
#         m = len(Col_Masvc)
#         n = 0
#         o = len(Col_Masvg)
#
#     # GENERADOR DE ELEMENTOS PILE - CONEXION
#     Exp_conex = np.zeros(len(MasaPilvc) + len(MasaPilvg))
#     F_CONEX = np.zeros((806,8))
#     # ContFil = 1
#     # Cont_Masvc = 1
#     # Cont_Masvg = 1
#     ContFil = 0
#     Cont_Masvc = 0
#     Cont_Masvg = 0
#
#     for i in range(len(MasaPilvc) + len(MasaPilvg)):
#         ContF = 0
#         case = OrdenPile[i]
#         if case == 1:  # El pilote está en el eje de la viga Cabezal
#             Z_condition = Pile_Vc[0:, Cont_Masvc] != 0
#             for j in range(m):
#                 valid_indices = np.where(Z_condition)[0]
#                 num_valid_indices = len(valid_indices)
#                 if num_valid_indices > 0:
#                     F_CONEX[ContFil:ContFil + num_valid_indices, 0] = Coord[Fil_Masvc[Cont_Masvc], Col_Masvc[j], 0]
#                     F_CONEX[ContFil:ContFil + num_valid_indices, 1] = Coord[Fil_Masvc[Cont_Masvc], Col_Masvc[j], 1]
#                     F_CONEX[ContFil:ContFil + num_valid_indices, 2] = Pile_Vc[valid_indices + 1, Cont_Masvc]
#                     F_CONEX[ContFil:ContFil + num_valid_indices, 3] = Coord[Fil_Masvc[Cont_Masvc], Col_Masvc[j], 0]
#                     F_CONEX[ContFil:ContFil + num_valid_indices, 4] = Coord[Fil_Masvc[Cont_Masvc], Col_Masvc[j], 1]
#                     F_CONEX[ContFil:ContFil + num_valid_indices, 5] = Pile_Vc[valid_indices, Cont_Masvc]
#                     F_CONEX[ContFil:ContFil + num_valid_indices, 6] = 2
#                     F_CONEX[ContFil:ContFil + num_valid_indices, 7] = 0
#                     ContFil += num_valid_indices
#                     ContF += num_valid_indices
#             Cont_Masvc += 1
#         elif case == 2:  # El pilote está en el eje de la viga de grúa
#             Z_condition = Pile_Vg[0:, Cont_Masvg] != 0
#             for j in range(n, o):
#                 valid_indices = np.where(Z_condition)[0]
#                 num_valid_indices = len(valid_indices)
#                 if num_valid_indices > 0:
#                     # F_CONEX[ContFil:ContFil + num_valid_indices, 0] = Coord[Fil_Masvg[Cont_Masvg], Col_Masvg[j], 0]
#                     # F_CONEX[ContFil:ContFil + num_valid_indices, 1] = Coord[Fil_Masvg[Cont_Masvg], Col_Masvg[j], 1]
#                     # F_CONEX[ContFil:ContFil + num_valid_indices, 2] = Pile_Vg[valid_indices + 1, Cont_Masvg]
#                     # F_CONEX[ContFil:ContFil + num_valid_indices, 3] = Coord[Fil_Masvg[Cont_Masvg], Col_Masvg[j], 0]
#                     # F_CONEX[ContFil:ContFil + num_valid_indices, 4] = Coord[Fil_Masvg[Cont_Masvg], Col_Masvg[j], 1]
#                     # F_CONEX[ContFil:ContFil + num_valid_indices, 5] = Pile_Vg[valid_indices, Cont_Masvg]
#                     F_CONEX[ContFil:ContFil + num_valid_indices, 0] = Coord['x_coords'][Fil_Masvg[Cont_Masvg]]
#                     F_CONEX[ContFil:ContFil + num_valid_indices, 1] = Coord['y_coords'][Col_Masvg[j]]
#                     F_CONEX[ContFil:ContFil + num_valid_indices, 2] = Pile_Vg[valid_indices + 1, Cont_Masvg]
#                     F_CONEX[ContFil:ContFil + num_valid_indices, 3] = Coord['x_coords'][Fil_Masvg[Cont_Masvg]]
#                     F_CONEX[ContFil:ContFil + num_valid_indices, 4] = Coord['y_coords'][Col_Masvg[j]]
#                     F_CONEX[ContFil:ContFil + num_valid_indices, 5] = Pile_Vg[valid_indices, Cont_Masvg]
#                     F_CONEX[ContFil:ContFil + num_valid_indices, 6] = 2
#                     F_CONEX[ContFil:ContFil + num_valid_indices, 7] = 0
#                     ContFil += num_valid_indices
#                     ContF += num_valid_indices
#             Cont_Masvg += 1
#         Exp_conex[i] = ContF - 1
#
#     return Exp_conex


def calcular_coordenadas_placa(Lx, Ly, Lz, N):
    num_filas = len(Ly)
    num_columnas = (N - 1) * 3 + 1

    # Calcular coordenadas en la dirección X
    x_coords = np.arange(0, num_columnas * Lx, Lx)

    # Calcular coordenadas en la dirección Y
    y_coords = np.cumsum(Ly)
    y_coords = np.insert(y_coords, 0, 0)  # Agregar 0 al principio

    # Establecer las coordenadas en la dirección Z
    z_coords = np.full((num_filas + 1, num_columnas), Lz)

    # Combinar coordenadas en una matriz

    Coord = {"x_coords": x_coords,
             "y_coords": y_coords,
             "z_coords": z_coords}

    # El problema es que en este ultimo paso no se genere el cubo de matrices con los valores. Solución averiguar como hacer este cubo de matrices en Python

    return Coord


def generar_elementos_pilote(ctr1, MasaPilvc, MasaPilvg, Col_Masvc, Col_Masvg, Coord, Fil_Masvc, Fil_Masvg, Pile_Vc, Pile_Vg, OrdenPile):


    if ctr1 == 0:
        m = 0
        n = 0
        o = 0
    elif ctr1 == 1:
        m = len(Col_Masvc)
        n = 0
        o = len(Col_Masvg)

    # GENERADOR DE ELEMENTOS PILE - CONEXION

        Exp_conex = np.zeros(len(MasaPilvc) + len(MasaPilvg))
        F_CONEX = np.zeros((806,8))
        ContFil = 0
        Cont_Masvc = 0
        Cont_Masvg = 0

    for i in range(len(MasaPilvc) + len(MasaPilvg)):
        ContF = 0
        case = OrdenPile[i]
        if case == 1:  # El pilote está en el eje de la Viga Cabezal
            for j in range(m):
                # for l in range(1):  # Marca Fila en matrix de coordenada Z
                    if Pile_Vc[1, Cont_Masvc] != 0:
                        F_CONEX[ContFil, 0] = Coord["x_coords"][Col_Masvc[j] - 1]
                        F_CONEX[ContFil, 1] = Coord["y_coords"][Fil_Masvc[Cont_Masvc] - 1]
                        F_CONEX[ContFil, 2] = Pile_Vc[1, Cont_Masvc]
                        F_CONEX[ContFil, 3] = Coord["x_coords"][Col_Masvc[j] - 1]
                        F_CONEX[ContFil, 4] = Coord["y_coords"][Fil_Masvc[Cont_Masvc] - 1]
                        F_CONEX[ContFil, 5] = Pile_Vc[0, Cont_Masvc]
                        F_CONEX[ContFil, 6] = 2  # DIRECCIÓN DEL EJE LOCAL 2 PARALELO AL $$
                        F_CONEX[ContFil, 7] = 0  # ROTACION DEL EJE LOCAL 2 ALRREDEDOR DEL EJE LOCAL 1
                        ContFil += 1
                        ContF += 1
            Cont_Masvc += 1
        elif case == 2:  # El pilote está en el eje de la Viga de grúa
            for j in range(n, o):
                # for l in range(1):  # Marca Fila en matrix de coordenada Z
                    if Pile_Vg[1, Cont_Masvg] != 0:
                        F_CONEX[ContFil, 0] = Coord["x_coords"][Col_Masvg[j] - 1]
                        F_CONEX[ContFil, 1] = Coord["y_coords"][Fil_Masvg[Cont_Masvg] - 1]
                        F_CONEX[ContFil, 2] = Pile_Vg[1, Cont_Masvg]
                        F_CONEX[ContFil, 3] = Coord["x_coords"][Col_Masvg[j] - 1]
                        F_CONEX[ContFil, 4] = Coord["y_coords"][Fil_Masvg[Cont_Masvg] - 1]
                        F_CONEX[ContFil, 5] = Pile_Vg[0, Cont_Masvg]
                        F_CONEX[ContFil, 6] = 2  # DIRECCIÓN DEL EJE LOCAL 2 PARALELO AL $$
                        F_CONEX[ContFil, 7] = 0  # ROTACION DEL EJE LOCAL 2 ALRREDEDOR DEL EJE LOCAL 1
                        ContFil += 1
                        ContF += 1
            Cont_Masvg += 1
        Exp_conex[i] = ContF

    # GENERADOR DE ELEMENTOS PILE - PILOTE

    Exp_pile = np.zeros(len(MasaPilvc) + len(MasaPilvg))
    F_PILE = np.zeros((19666, 8))
    ContFil = 0
    Cont_Masvc = 0
    Cont_Masvg = 0

    for i in range(len(MasaPilvc) + len(MasaPilvg)):
        print(f"i:{i}")
        ContF = 0
        case = OrdenPile[i]
        if case == 1:  # El pilote está en el eje de la Viga Cabezal
            for j in range(m):
                for l in range(1, len(Pile_Vc[:, 0]) - 1):  # Marca Fila en matrix de coordenada Z
                    if Pile_Vc[l + 1, Cont_Masvc] != 0:
                        F_PILE[ContFil, 0] = Coord["x_coords"][Col_Masvc[j] - 1]
                        F_PILE[ContFil, 1] = Coord["y_coords"][Fil_Masvc[Cont_Masvc] - 1]
                        F_PILE[ContFil, 2] = Pile_Vc[l + 1, Cont_Masvc]
                        F_PILE[ContFil, 3] = Coord["x_coords"][Col_Masvc[j] - 1]
                        F_PILE[ContFil, 4] = Coord["y_coords"][Fil_Masvc[Cont_Masvc] - 1]
                        F_PILE[ContFil, 5] = Pile_Vc[l, Cont_Masvc]
                        F_PILE[ContFil, 6] = 2  # DIRECCIÓN DEL EJE LOCAL 2 PARALELO AL $$
                        F_PILE[ContFil, 7] = 0  # ROTACION DEL EJE LOCAL 2 ALRREDEDOR DEL EJE LOCAL 1
                        ContFil += 1
                        ContF += 1

                        print(f"filas vc:¨{ContFil}, i:{i}")

                Cont_Masvc += 1
        elif case == 2:  # El pilote está en el eje de la Viga de grúa
            for j in range(n, o):
                for l in range(1, len(Pile_Vg[:, 0]) - 1):  # Marca Fila en matrix de coordenada Z
                    if Pile_Vg[l + 1, Cont_Masvg] != 0:
                        F_PILE[ContFil, 0] = Coord["x_coords"][Col_Masvg[j] - 1]
                        F_PILE[ContFil, 1] = Coord["y_coords"][Fil_Masvg[Cont_Masvg] - 1]
                        F_PILE[ContFil, 2] = Pile_Vg[l + 1, Cont_Masvg]
                        F_PILE[ContFil, 3] = Coord["x_coords"][Col_Masvg[j] - 1]
                        F_PILE[ContFil, 4] = Coord["y_coords"][Fil_Masvg[Cont_Masvg] - 1]
                        F_PILE[ContFil, 5] = Pile_Vg[l, Cont_Masvg]
                        F_PILE[ContFil, 6] = 2  # DIRECCIÓN DEL EJE LOCAL 2 PARALELO AL $$
                        F_PILE[ContFil, 7] = 0  # ROTACION DEL EJE LOCAL 2 ALRREDEDOR DEL EJE LOCAL 1
                        ContFil += 1
                        ContF += 1

                        print(f"filas vg:¨{ContFil}, i:{i}")
            Cont_Masvg += 1
        Exp_pile[i] = ContF


    filePiles = {"F_CONEX": F_CONEX,
                 "Exp_conex": Exp_conex,
                 "F_PILE": F_PILE,
                 "Exp_pile": Exp_pile
                     }

    return filePiles

def linkFun(NvTERR, num_pile, Fil_Masv, Pend, LinkSpac, Apoyos_Pil, Ly, LinkLim):
    Link_V = []
    NlinkV = []

    Link_V.append(NvTERR - np.sum(Ly[:Fil_Masv[num_pile] - 1]) / Pend - 0.5 * LinkSpac[0])
    H = Link_V[0]

    cont = 1
    Limit = 0
    contLk = 1

    while H > Apoyos_Pil:
        H = H - LinkSpac[Limit]

        if abs(H - Link_V[0] - 0.5 * LinkSpac[0]) - 1e-5 <= np.sum(LinkLim[:Limit + 1]):
            Link_V.append(H)
            cont += 1
            contLk += 1
        elif Limit < len(LinkLim) - 1:
            H = Link_V[cont - 1] - 0.5 * (LinkSpac[Limit] + LinkSpac[Limit + 1])
            Link_V.append(H)
            NlinkV.append(contLk)
            Limit += 1
            cont += 1
            contLk = 1
        else:
            H = Link_V[cont - 1] - LinkSpac[Limit]
            Link_V.append(H)
            cont += 1
            contLk += 1

    Link_V[-1] = Apoyos_Pil
    NlinkV.append(contLk-1)
    Link_V = np.array(Link_V)
    NlinkV = np.array(NlinkV)

    return Link_V, NlinkV

def linkDic(MasaPil, Fil_Masv, Apoyos_Pil, LinkSpac, LinkLim, Ly, NvTERR, Pend):

    VarLink = {}
    Nlink = {}

    for i in range(0,len(MasaPil)):
        varLink, nlink = linkFun(NvTERR, i, Fil_Masv, Pend, LinkSpac, Apoyos_Pil[i], Ly, LinkLim)

        VarLink[f"{i}"] =  varLink
        Nlink[f"{i}"] = nlink

    return VarLink, Nlink

