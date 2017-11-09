from cid import CID


class MS2Spectrum:
    def __init__(self, mz: float):
        self.mz = mz
        self.cids = {}

    def add_cid(self, id: str, energy: str, rt: float, masses: list, intensities: list):
        self.cids[energy] = CID(id, energy, rt, masses, intensities)

    def __str__(self):
        return 'MS2 spectrum (mz={0}, CIDs={1})'.format(self.mz, self.cids)

    def __repr__(self):
        return self.__str__()
