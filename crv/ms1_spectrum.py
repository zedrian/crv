class MS1Spectrum:
    def __init__(self, id: str, mz_min: float, mz_max: float, rt: float, masses: list, intensities: list):
        self.id = id
        self.mz_min = mz_min
        self.mz_max = mz_max
        self.rt = rt
        self.masses = masses
        self.intensities = intensities

    def __str__(self):
        return 'MS1 spectrum (id={0}, mz={1}-{2}, RT={3}, masses={4}, intensities={5})'.format(self.id, self.mz_min, self.mz_max, self.rt, self.masses, self.intensities)

    def __repr__(self):
        return self.__str__()
