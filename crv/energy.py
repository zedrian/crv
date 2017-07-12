from crv.fragment import Fragment


class Energy:
    def __init__(self, name: str, fragments: list):
        self.name = name
        self.fragments = fragments

    # Needed for JSON serialization.
    def get_dictionary(self) -> dict:
        dictionary = dict()
        dictionary['Name'] = self.name
        dictionary['Fragments'] = []
        for fragment in self.fragments:
            dictionary['Fragments'].append(fragment.get_dictionary())
        return dictionary

    # Needed for JSON deserialization
    @classmethod
    def from_dictionary(cls, dictionary: dict):
        name = dictionary.get('Name')
        fragment_dictionaries = dictionary.get('Fragments')
        fragments = []
        for fragment in fragment_dictionaries:
            fragments.append(Fragment.from_dictionary(fragment))
        return Energy(name, fragments)

    def __str__(self):
        return 'Name: {0}\n'.format(self.name) + \
               'Fragments: {0}\n'.format(self.fragments) + \
               '\n'

    def __repr__(self):
        return self.__str__()