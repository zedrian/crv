from crv.fragment import Fragment


class FragmentMatch:
    def __init__(self, fragment_smiles: Fragment, compound_names: list = []):
        self.fragment_smiles = fragment_smiles
        self.compound_names = compound_names

    # Needed for JSON serialization.
    def get_dictionary(self) -> dict:
        dictionary = dict()
        dictionary['Fragment'] = self.fragment_smiles
        dictionary['Compounds'] = self.compound_names
        return dictionary

    # Needed for JSON deserialization.
    @classmethod
    def from_dictionary(cls, dictionary: dict):
        fragment_smiles = dictionary.get('Fragment')
        compound_names = dictionary.get('Compounds')
        return FragmentMatch(fragment_smiles, compound_names)

    def __str__(self):
        return 'FragmentMatch(fragment=\'{0}\', compounds={1})'.format(self.fragment_smiles, self.compound_names)

    def __repr__(self):
        return self.__str__()
