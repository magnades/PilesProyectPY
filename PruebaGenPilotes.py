import numpy as np


def generate_pilescoord(Apoyos_Pil, Link_V, NlinkV, LinkSpac, Lz, LongDes, PileDiv):
    num_pilotes = Link_V.__len__()

    Pile_V = {}

    for i in range(num_pilotes):
        pile_V = []
        linkVal = Link_V[f"{i}"]
        nlinkVal = NlinkV[f"{i}"]

        pile_V.append(Lz)
        pile_V.append(Lz - LongDes)
        Limit = 0

        for k in range(int(np.sum(nlinkVal) - 1)):
            if linkVal[k + 1] != 0:
                if k <= np.sum(nlinkVal[:Limit+1]) - 1:
                    pass
                elif Limit < len(LinkSpac)-1:
                    Limit += 1

                Space = np.linspace(linkVal[k], linkVal[k + 1], PileDiv[Limit] + 1)

                for l in range(int(len(Space) - 1)):
                    pile_V.append(Space[l])

        pile_V.append(Space[-1])
        pile_V.append(Apoyos_Pil[i])
        Pile_V[f"{i}"] = np.array(pile_V)

    return Pile_V

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

    return Coord


def generate_piles(ctr1, MasaPilvc, MasaPilvg, Col_Masvc, Col_Masvg, Coord, Fil_Masvc, Fil_Masvg, Pile_Vc, Pile_Vg, OrdenPile):

    if ctr1 == 0:
        m = 0
        n = 0
        o = 0
    elif ctr1 == 1:
        m = len(Col_Masvc)
        n = 0
        o = len(Col_Masvg)

    # GENERADOR DE ELEMENTOS PILE - CONEXION
        num_piles = MasaPilvc.__len__() + MasaPilvg.__len__()
        Exp_conex = np.zeros(num_piles)

        F_CONEX = []
        Cont_Masvc = 0
        Cont_Masvg = 0

    for i in range(num_piles):
        ContF = 0
        case = OrdenPile[i]
        if case == 1:  # El pilote está en el eje de la Viga Cabezal
            for j in range(m):
                    if Pile_Vc[f"{Cont_Masvc}"][1] != 0:
                        F_CONEX.append(
                            [Coord["x_coords"][Col_Masvc[j] - 1],
                             Coord["y_coords"][Fil_Masvc[Cont_Masvc] - 1],
                             Pile_Vc[f"{Cont_Masvc}"][1],
                             Coord["x_coords"][Col_Masvc[j] - 1],
                             Coord["y_coords"][Fil_Masvc[Cont_Masvc] - 1],
                             Pile_Vc[f"{Cont_Masvc}"][0],
                             2.,  # DIRECCIÓN DEL EJE LOCAL 2 PARALELO AL $$
                             0. ]) # ROTACION DEL EJE LOCAL 2 ALRREDEDOR DEL EJE LOCAL 1
                        ContF += 1
            Cont_Masvc += 1
        elif case == 2:  # El pilote está en el eje de la Viga de grúa
            for j in range(n, o):
                    if Pile_Vg[f"{Cont_Masvg}"][1] != 0:
                        F_CONEX.append(
                            [Coord["x_coords"][Col_Masvg[j] - 1],
                             Coord["y_coords"][Fil_Masvg[Cont_Masvg] - 1],
                             Pile_Vc[f"{Cont_Masvg}"][1],
                             Coord["x_coords"][Col_Masvg[j] - 1],
                             Coord["y_coords"][Fil_Masvg[Cont_Masvg] - 1],
                             Pile_Vc[f"{Cont_Masvg}"][0],
                             2,  # DIRECCIÓN DEL EJE LOCAL 2 PARALELO AL $$
                             0 ]) # ROTACION DEL EJE LOCAL 2 ALRREDEDOR DEL EJE LOCAL 1
                        ContF += 1
            Cont_Masvg += 1
        Exp_conex[i] = ContF

    # GENERADOR DE ELEMENTOS PILE - PILOTE

    Exp_pile = np.zeros(num_piles)
    F_PILE = []
    Cont_Masvc = 0
    Cont_Masvg = 0

    for i in range(num_piles):
        ContF = 0
        case = OrdenPile[i]
        if case == 1:  # El pilote está en el eje de la Viga Cabezal
            for j in range(m):
                lenPile = len(Pile_Vc[f"{Cont_Masvc}"])
                for l in range(1,lenPile - 1):  # Marca Fila en matrix de coordenada Z
                    if Pile_Vc[f"{Cont_Masvc}"][1] != 0:
                        F_PILE.append(
                            [Coord["x_coords"][Col_Masvc[j] - 1],
                            Coord["y_coords"][Fil_Masvc[Cont_Masvc] - 1],
                            Pile_Vc[f"{Cont_Masvc}"][l + 1],
                            Coord["x_coords"][Col_Masvc[j] - 1],
                            Coord["y_coords"][Fil_Masvc[Cont_Masvc] - 1],
                            Pile_Vc[f"{Cont_Masvc}"][l],
                            2.,  # DIRECCIÓN DEL EJE LOCAL 2 PARALELO AL $$
                            0.])  # ROTACION DEL EJE LOCAL 2 ALRREDEDOR DEL EJE LOCAL 1
                        ContF += 1
            Cont_Masvc += 1

        elif case == 2:  # El pilote está en el eje de la Viga de grúa
            for j in range(n, o):
                lenPile = len(Pile_Vg[f"{Cont_Masvg}"])
                for l in range(1,lenPile - 1):  # Marca Fila en matrix de coordenada Z
                    if Pile_Vg[f"{Cont_Masvg}"][1] != 0:
                        F_PILE.append(
                            [Coord["x_coords"][Col_Masvg[j] - 1],
                            Coord["y_coords"][Fil_Masvg[Cont_Masvg] - 1],
                            Pile_Vg[f"{Cont_Masvg}"][l + 1],
                            Coord["x_coords"][Col_Masvg[j] - 1],
                            Coord["y_coords"][Fil_Masvg[Cont_Masvg] - 1],
                            Pile_Vg[f"{Cont_Masvg}"][l],
                            2.,  # DIRECCIÓN DEL EJE LOCAL 2 PARALELO AL $$
                            0.])  # ROTACION DEL EJE LOCAL 2 ALRREDEDOR DEL EJE LOCAL 1
                        ContF += 1
            Cont_Masvg += 1
        Exp_pile[i] = ContF


    filePiles = {"F_CONEX": np.array(F_CONEX),
                 "Exp_conex": np.array(Exp_conex),
                 "F_PILE": np.array(F_PILE),
                 "Exp_pile": np.array(Exp_pile)
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

def generate_links(MasaPil, Fil_Masv, Apoyos_Pil, LinkSpac, LinkLim, Ly, NvTERR, Pend):

    VarLink = {}
    Nlink = {}

    for i in range(0,len(MasaPil)):
        varLink, nlink = linkFun(NvTERR, i, Fil_Masv, Pend, LinkSpac, Apoyos_Pil[i], Ly, LinkLim)

        VarLink[f"{i}"] =  varLink
        Nlink[f"{i}"] = nlink

    return VarLink, Nlink


def generate_frames(Coord, Col_Masvc, Fil_Masvg):
    # Generación de archivos de frames para vigas de cabezales (F_VCAB)
    Col_Masvc -= 1
    F_VCAB = []
    for i in Col_Masvc:
        for j in range(len(Coord["y_coords"]) - 1):
            frame_data = [
                Coord["x_coords"][i],
                Coord["y_coords"][j],
                Coord["z_coords"][j,i],
                Coord["x_coords"][i],
                Coord["y_coords"][j + 1],
                Coord["z_coords"][j + 1,i],
                3.,
                0
            ]
            F_VCAB.append(frame_data)

    # Generación de archivos de frames para vigas de grúa (F_VGRUA)
    Fil_Masvg -= 1
    F_VGRUA = []
    for i in Fil_Masvg:
        for j in range(len(Coord["x_coords"]) - 1):
            frame_data = [
                Coord["x_coords"][j],
                Coord["y_coords"][i],
                Coord["z_coords"][i, j],
                Coord["x_coords"][j + 1],
                Coord["y_coords"][i],
                Coord["z_coords"][i, j + 1],
                3.,
                0
            ]
            F_VGRUA.append(frame_data)

    # Generación de archivos de frames para vigas de borde (F_VBORDE)

    Fil_VB = [len(Coord["y_coords"])-1]
    F_VBORDE = []

    for i in Fil_VB:
        for j in range(len(Coord["x_coords"])-1):
            frame_data = [
                Coord["x_coords"][j],
                Coord["y_coords"][i],
                Coord["z_coords"][i, j],
                Coord["x_coords"][j + 1],
                Coord["y_coords"][i],
                Coord["z_coords"][i, j + 1],
                3.,
                0.
            ]
            F_VBORDE.append(frame_data)

    Frames = {"F_VCAB": np.array(F_VCAB),
            "F_VGRUA": np.array(F_VGRUA),
            "F_VBORDE": np.array(F_VBORDE)}

    return Frames


