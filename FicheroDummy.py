import numpy as np
from PruebaGenPilotes import linkDic,linkFun, generate_piles
from itertools import zip_longest

NvTERR = -0.2457
Ly=np.array([1.5,2.3,2.3,3.05,3.05,3.05,3.05,3.05,3.05,2.3,2.3,1.5,1.5,1.5])
Fil_Masv = np.array([4,6,8,10])
Pend = 2
LinkSpac= np.array([12,12,12,24,48])*0.0254
LinkLim = np.array([126/3,126/3,126/3,246,390])*0.0254
Apoyos_Pil = np.array([-17.,-19.,-21.,-23.])
MasaPilvc = np.array([2.246,3.124,4.002,4.88])*9.81
Lz=1.75
LongDes = 2.
PileDiv = np.array([1, 1, 1, 1, 1])
# num_pilotes = 4
#
# Link_V = np.zeros((1, 4))
# NlinkV = np.zeros((len(LinkLim),num_pilotes))

# def linkFun(NvTERR, num_pile, Fil_Masv, Pend, LinkSpac, Apoyos_Pil, Ly, LinkLim):
#     Link_V = []
#     NlinkV = []
#
#     Link_V.append(NvTERR - np.sum(Ly[:Fil_Masv[num_pile] - 1]) / Pend - 0.5 * LinkSpac[0])
#     H = Link_V[0]
#
#     cont = 1
#     Limit = 0
#     contLk = 1
#
#     while H > Apoyos_Pil:
#         H = H - LinkSpac[Limit]
#
#         if abs(H - Link_V[0] - 0.5 * LinkSpac[0]) - 1e-5 <= np.sum(LinkLim[:Limit + 1]):
#             Link_V.append(H)
#             cont += 1
#             contLk += 1
#         elif Limit < len(LinkLim) - 1:
#             H = Link_V[cont - 1] - 0.5 * (LinkSpac[Limit] + LinkSpac[Limit + 1])
#             Link_V.append(H)
#             NlinkV.append(contLk)
#             Limit += 1
#             cont += 1
#             contLk = 1
#         else:
#             H = Link_V[cont - 1] - LinkSpac[Limit]
#             Link_V.append(H)
#             cont += 1
#             contLk += 1
#
#     Link_V[-1] = Apoyos_Pil
#     NlinkV.append(contLk-1)
#     Link_V = np.array(Link_V)
#     NlinkV = np.array(NlinkV)
#
#     return Link_V, NlinkV
#
# def linkDic(MasaPil, Fil_Masv, Apoyos_Pil, LinkSpac, LinkLim, Ly, NvTERR, Pend):
#     num_pil = len(MasaPil)
#     # VarLink = np.empty((0, num_pil))
#     # Nlink = np.empty((0, num_pil))
#     VarLink = {}
#     Nlink = {}
#
#     for i in range(0,len(MasaPil)):
#         varLink, nlink = linkFun(NvTERR, i, Fil_Masv, Pend, LinkSpac, Apoyos_Pil[i], Ly, LinkLim)
#
#         VarLink[f"{i}"] =  varLink
#         Nlink[f"{i}"] = nlink
#
#     # VarLink = np.array(list(VarLink))
#     # Nlink = np.array(list(Nlink))
#     return VarLink, Nlink

# def generate_piles(Apoyos_Pil, Link_V, NlinkV, LinkSpac, Lz, LongDes, PileDiv):
#     num_pilotes = Link_V.__len__()
#
#     Pile_V = {}
#
#     for i in range(num_pilotes):
#         print(f"{i}")
#         pile_V = []
#         linkVal = Link_V[f"{i}"]
#         nlinkVal = NlinkV[f"{i}"]
#
#         pile_V.append(Lz)
#         pile_V.append(Lz - LongDes)
#         cont = 2
#         Limit = 0
#
#         for k in range(int(np.sum(nlinkVal) - 1)):
#             if linkVal[k + 1] != 0:
#                 if k <= np.sum(nlinkVal[:Limit+1]) - 1:
#                     pass
#                 elif Limit < len(LinkSpac)-1:
#                     Limit += 1
#
#                 Space = np.linspace(linkVal[k], linkVal[k + 1], PileDiv[Limit] + 1)
#
#                 for l in range(int(len(Space) - 1)):
#                     pile_V.append(Space[l])
#                     cont += 1
#
#         pile_V.append(Space[-1])
#         pile_V.append(Apoyos_Pil[i])
#         Pile_V[f"{i}"] = np.array(pile_V)
#
#     return Pile_V




VarLink, Nlink = linkDic(MasaPilvc, Fil_Masv, Apoyos_Pil, LinkSpac, LinkLim, Ly, NvTERR, Pend)

Pile_Vc = generate_piles(Apoyos_Pil, VarLink, Nlink, LinkSpac, Lz, LongDes, PileDiv)


print(Pile_Vc)