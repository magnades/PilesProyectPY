from Element import *
import numpy as np

class genPile(Element):
    """
    class: Generator of a single pile with the same crosssection given the following variables:

    - An element previously defined of class Element.py
    - Coordinates of head node
    - Coordinates of bottom node
    - Coordinates of soil node
    """

    def __init__(self, node0, node1, material):
        self.nodes = [node0, node1]
        self.material = material
#        self.force = 0.0
#        self.Forces = [np.zeros(2), np.zeros(2), np.zeros(2)]
#        self.Kt = [[np.zeros((2, 2)), np.zeros((2, 2)), np.zeros((2, 2))],
#                   [np.zeros((2, 2)), np.zeros((2, 2)), np.zeros((2, 2))],
#                   [np.zeros((2, 2)), np.zeros((2, 2)), np.zeros((2, 2))]]

    def __str__(self):
        s = \
            """Truss: node {} to node {}:
               material properties: {}  strain:{}   stress:{}  
               internal force: {}
               Pe: [ {} {} ]""".format(self.nodes[0].index, self.nodes[1].index,
                                       repr(self.material), self.material.getStrain(),
                                       self.material.getStress(),
                                       self.force, *self.Forces[1])
        return s

    def __repr__(self):
        return "Truss({},{},{})".format(repr(self.nodes[0]),
                                        repr(self.nodes[1]),
                                        repr(self.material))

