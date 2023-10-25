class Material():
    """
    class: representing a basic material properties
            E, A, nu, fy

    """

    def __init__(self, params={'E':1.0, 'A':1.0, 'nu':0.0, 'fy':1.0e30}):
        self.parameters = params

        # make sure all necessary parameters exist
        if 'E' not in self.parameters:
            self.parameters['E']  = 1.0
        if 'A' not in self.parameters:
            self.parameters['A']  = 1.0
        if 'nu' not in self.parameters:
            self.parameters['nu'] = 0.2
        if 'fy' not in self.parameters:
            self.parameters['fy'] = 1.0e30

    def __str__(self):

        par = self.parameters
        s= "E:{}  A:{} nu:{} fy:{}".format(par['E'], par['A'], par['nu'], par['fy'])
        return s

    def __repr__(self):
        s = "{}({})".format(self.__class__.__name__, self.parameters)
        return s

