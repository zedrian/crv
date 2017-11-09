class ExperimentalCompound:
    class Match:
        def __init__(self, mass: float, rt: float, precursor: float, area: float, score: float):
            self.mass = mass
            self.rt = rt
            self.precursor = precursor
            self.area = area
            self.score = score
            self.ms1_spectrum = None
            self.ms2_spectra = []

    def __init__(self, name: str):
        self.name = name
        self.matches = []

    def add_match(self, mass: float, rt: float, precursor: float, area: float, score: float):
        self.matches.append(ExperimentalCompound.Match(mass, rt, precursor, area, score))
