class Fragment:
    def __init__(self, mass: float, intensity: float, smiles: str = ''):
        self.mass = mass
        self.intensity = intensity
        self.smiles = smiles

    def __str__(self):
        return 'Fragment(mass={0}, intensity={1}, SMILES=\'{2}\')'.format(self.mass, self.intensity, self.smiles)

    def __repr__(self):
        return self.__str__()