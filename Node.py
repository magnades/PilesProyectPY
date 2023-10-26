import numpy as np


class Node():
    """
    class: representing a single Node
    """

    def __init__(self, x, y, z, u=0, v=0, w=0):
        self.pos      = np.array([x,y,z])
        self.index    = -1
        self.disp     = np.array([u,v,w])
        self.fixity   = [False, False, False]
        self.force    = np.zeros(3)
        self._hasLoad = False

    def __str__(self):
        s = \
"""Node {}:
   x:{}   y:{}  z:{}
   fix:{} fix:{} fix:{}
   Px:{}  Py:{} Pz:{}
   u:{}   v:{}  w:{}""".format( self.index,
                           self.pos[0],   self.pos[1], self.pos[2],
                           *self.fixity,
                           self.force[0], self.force[1], self.force[2],
                           self.disp[0],  self.disp[1],  self.disp[2])
        return s

    def __repr__(self):
        return "Node({},{},{},u={},v={},w={})".format(self.pos[0],self.pos[1],self.pos[2],self.disp[0],self.disp[1],self.disp[2])

