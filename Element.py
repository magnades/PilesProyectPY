import numpy as np
from Node import *
from Material import *

class Element():
    """
    class: representing a single element
    NameElement (NameNode0, NameNode1, NameMaterial)
    NameNode are elements of class Node()
    NameMaterial is element of class Material()
    Nodes and elements are defined previously
    """

    def __init__(self, node0, node1, material):
        self.nodes    = [node0, node1]
        self.material = material
        self.force    = 0.0
        self.Forces   = [ np.zeros(2), np.zeros(2) ]
        self.Kt       = [ [np.zeros((2,2)), np.zeros((2,2))], [np.zeros((2,2)), np.zeros((2,2))] ]

    def __str__(self):
        s = \
"""{}: node {} to node {}:
   material properties: {} """.format(self.__class__,
                            self.nodes[0].index, self.nodes[1].index,
                            repr(self.material))
        return s

    def __repr__(self):
        return "{}({},{},{})".format(self.__class__,
                                     repr(self.nodes[0]),
                                     repr(self.nodes[1]),
                                     repr(self.material))

