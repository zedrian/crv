from fragment import Fragment


class CID:
    def __init__(self, id: str, energy: str, rt: float, masses: list, intensities: list):
        self.id = id
        self.energy = energy
        self.rt = rt

        self.fragments = []

        # check that numbers of masses and intensities are the same
        if len(masses) != len(intensities):
            print('=============================================')
            print('CID with ID={0} has {1} masses and {2} intensities.'.format(id, len(masses), len(intensities)))
            exit()

        # fill list of fragments
        for i in range(0, len(masses)):
            self.fragments.append(Fragment(masses[i], intensities[i]))

    def __str__(self):
        return 'CID (id={0}, energy={1}, RT={2}, fragments={3})'.format(self.id, self.energy, self.rt, self.fragments)

    def __repr__(self):
        return self.__str__()
