class Compound:
    def __init__(self, name: str = '', smiles: str = '', energies: list = []):
        self.name = name
        self.smiles = smiles
        self.energies = energies

    # Needed for JSON serialization.
    def get_dictionary(self) -> dict:
        dictionary = dict()
        dictionary['Name'] = self.name
        dictionary['SMILES'] = self.smiles
        dictionary['Energies'] = []
        for energy in self.energies:
            dictionary['Energies'].append(energy.get_dictionary())
        return dictionary

    def __str__(self):
        return 'Name: {0}\n'.format(self.name) + \
               'SMILES: {0}\n'.format(self.smiles) + \
               'Energies: {0}\n'.format(self.energies) + \
               '\n'

    def __repr__(self):
        return self.__str__()
