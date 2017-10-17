class MS2Spectrum:
    def __init__(self, id: str, mz: float, rt: float, energy: str, masses: tuple, intensities: tuple):
        self.id = id
        self.mz = mz
        self.rt = rt
        self.energy = energy
        self.masses = list(masses)
        self.intensities = list(intensities)

    def __str__(self):
        return 'MS2 spectrum (id={0}, mz={1}, RT={2}, energy={3}, masses={4}, intensities={5})'.format(self.id, self.mz, self.rt, self.energy, self.masses, self.intensities)

    def __repr__(self):
        return self.__str__()
