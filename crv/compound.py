class Compound:
    def __init__(self, name: str = '', smiles: str = '', energies: list = []):
        self.name = name
        self.smiles = smiles
        self.energies = energies

    def __str__(self):
        return 'Name: {0}\n'.format(self.name) + \
               'SMILES: {0}\n'.format(self.smiles) + \
               'Energies: {0}\n'.format(self.energies) + \
               '\n'

    def __repr__(self):
        return self.__str__()
