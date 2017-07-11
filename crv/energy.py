class Energy:
    def __init__(self, name: str, fragments: list):
        self.name = name
        self.fragments = fragments

    def __str__(self):
        return 'Name: {0}\n'.format(self.name) + \
               'Fragments: {0}\n'.format(self.fragments) + \
               '\n'

    def __repr__(self):
        return self.__str__()
