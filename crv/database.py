from json import dumps, loads
from hashlib import sha512
from os import listdir
from os.path import isfile, join

from crv.compound import Compound
from crv.energy import Energy
from crv.fragment import Fragment
from crv.fragment_match import FragmentMatch
from crv.parameters import Parameters
from crv.utility import show_progress


# Class that is needed to store compounds and corresponding metadata.
class Database:
    # Is trying to load from previously generated database.json file.
    # Uses SHA-512 digest to track changes in folder.
    def __init__(self, parameters: Parameters):
        root_folder = 'databases'
        self.folder_name = join(root_folder, parameters.database_name)
        self.digest = None
        self.compounds = []
        self.fragment_matches = []

        database_file_name = join(self.folder_name, 'database.json')

        old_json_data = None
        if isfile(database_file_name):
            print('Previously initialized database found.')

            old_database_file = open(database_file_name)
            json_string = ''.join(old_database_file.readlines())
            old_database_file.close()

            old_json_data = loads(json_string)
            old_digest = old_json_data.get('SHA-512')
            new_digest = compute_digest(self.folder_name)

            if old_digest != new_digest:
                print('Current SHA-512 digest differs from stored in database. Re-initializing.')
                can_load_from_json = False
            else:
                print('SHA-512 digest is actual. Can load database from cached file.')
                can_load_from_json = True

            self.digest = new_digest
        else:
            print('Data base was not initialized before. Initializing.')
            can_load_from_json = False

        if can_load_from_json:
            self.load(old_json_data)
        else:
            self.initialize(parameters)
            self.save()


    def initialize(self, parameters: Parameters):
        self.digest = compute_digest(self.folder_name)

        compounds_file_name = join(self.folder_name, parameters.database_compounds_file_name)
        self.compounds = load_compound_skeletons(compounds_file_name)
        print('Compounds to load: {0}'.format(len(self.compounds)))

        check_missing_compound_files(self.compounds, self.folder_name)
        load_compound_energies(self.compounds, self.folder_name, parameters.minimal_intensity)

        self.fragment_matches = compute_fragment_matches(self.compounds)

    def load(self, json_data):
        digest = json_data.get('SHA-512')
        self.digest = digest

        # process compounds

        compound_dictionaries = json_data.get('Compounds')

        progress_label = 'Loading compounds: '
        show_progress(progress_label, 0.0)
        compound_index = 0

        for compound in compound_dictionaries:
            self.compounds.append(Compound.from_dictionary(compound))

            compound_index += 1
            show_progress(progress_label, float(compound_index) / float(len(compound_dictionaries)))

        # process fragment matches

        fragment_match_dictionaries = json_data.get('Fragment matches')

        progress_label = 'Loading fragment matches: '
        show_progress(progress_label, 0.0)
        match_index = 0

        for match in fragment_match_dictionaries:
            self.fragment_matches.append(FragmentMatch.from_dictionary(match))

            match_index += 1
            show_progress(progress_label, float(match_index) / float(len(fragment_match_dictionaries)))

        print('Database successfully loaded ({0} compounds, {1} fragments).'.format(len(self.compounds),
                                                                                    len(self.fragment_matches)))

    def save(self):
        json_data = dict()
        json_data['SHA-512'] = self.digest

        json_data['Compounds'] = []
        for compound in self.compounds:
            json_data['Compounds'].append(compound.get_dictionary())

        json_data['Fragment matches'] = []
        for match in self.fragment_matches:
            json_data['Fragment matches'].append(match.get_dictionary())

        json_string = dumps(json_data, indent=2)
        file = open(join(self.folder_name, 'database.json'), 'w')
        file.write(json_string)
        file.close()

        print('Database successfully saved.')


# Loads compound names and SMILES from file and returns list of compounds.
def load_compound_skeletons(file_name: str) -> list:
    if not isfile(file_name):
        print('Failed to load compounds from \'{0}\': file not exists.'.format(file_name))
        exit()

    compounds = []
    file = open(file_name)
    for line in file:
        components = line[:-1].split(';')
        compounds.append(Compound(name=components[0], smiles=components[1]))
    file.close()

    return compounds


# Computes total digest for all files in folder except database.json.
def compute_digest(folder_name) -> str:
    hash = sha512()
    file_names = [join(folder_name, f) for f in listdir(folder_name) if 'database.json' not in f]
    for name in file_names:
        file = open(name)
        data = file.read()
        encoded = data.encode()
        hash.update(encoded)
        file.close()
    return hash.hexdigest()


# Checks energy files existence for each compound from list.
# If any, shows list of missing compound files and exits.
def check_missing_compound_files(compounds: list, folder_name: str):
    missing_files = []

    for compound in compounds:
        file_name = join(folder_name, compound.name)
        if not isfile(file_name):
            missing_files.append(compound.name)

    if len(missing_files) != 0:
        print('Cannot find energy files for {0} compounds:'.format(len(missing_files)))
        for name in missing_files:
            print('- {0}'.format(name))
        print('You can either provide proper energy files or delete them from compounds list.')
        exit()

    print('All needed energy files exist in database.')


# Loads compound energies from corresponding raw energy files.
# All needed files should exist.
def load_compound_energies(compounds: list, folder_name: str, minimal_intensity: float):
    progress_label = 'Loading compound energies: '
    show_progress(progress_label, 0.0)

    compound_index = 0
    for compound in compounds:
        energy_file = open(join(folder_name, compound.name))
        energies = []

        smiles_started = False
        for line in energy_file.readlines():
            line = line[:-1]
            if line == '':
                smiles_started = True
                continue
            if smiles_started:
                components = line.split(' ')
                fragment_index = int(components[0])
                smiles = components[2]
                for energy in energies:
                    if fragment_index < len(energy.fragments):
                        energy.fragments[fragment_index].smiles = smiles
            else:
                if 'energy' in line:
                    energies.append(Energy(name=line.replace('en', 'En'), fragments=[]))
                else:
                    components = line.split(' ')
                    fragment = Fragment(mass=float(components[0]), intensity=float(components[1]))
                    if fragment.intensity < minimal_intensity:
                        continue
                    energies[-1].fragments.append(fragment)
        energy_file.close()
        compound.energies = energies

        compound_index += 1
        show_progress(progress_label, float(compound_index) / float(len(compounds)))


# Creates list of FragmentMatch objects describing lists of compounds that contain each fragment.
# Returns such list sorted by frequency (unique fragments first).
def compute_fragment_matches(compounds: list) -> list:
    matches = []

    progress_label = 'Computing fragment matches: '
    show_progress(progress_label, 0.0)

    compound_index = 0
    for compound in compounds:
        # we don't want to process multiple times the same fragment from different energies
        processed_fragment_smiles = []
        for energy in compound.energies:
            for fragment in energy.fragments:
                if fragment.smiles in processed_fragment_smiles:
                    continue

                # try to find existing match with such fragment
                match_found = False
                for match in matches:
                    if match.fragment_smiles == fragment.smiles:
                        # match found, store compound
                        match.compound_names.append(compound.name)
                        match_found = True
                        break

                if not match_found:
                    # create new match
                    match = FragmentMatch(fragment.smiles, [compound.name])
                    matches.append(match)

                processed_fragment_smiles.append(fragment.smiles)

        compound_index += 1
        show_progress(progress_label, float(compound_index) / float(len(compounds)))

    # sort matches (unique fragments first)
    matches.sort(key=lambda x: len(x.compound_names))

    return matches
