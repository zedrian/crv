from json import dumps, loads
from hashlib import sha512
from os import listdir
from os.path import isfile, join

from crv.compound import Compound
from crv.energy import Energy
from crv.fragment import Fragment
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

    def initialize(self, parameters: Parameters):
        self.digest = compute_digest(self.folder_name)

        compounds_file_name = join(self.folder_name, parameters.database_compounds_file_name)
        self.compounds = load_compound_skeletons(compounds_file_name)
        print('Compounds to load: {0}'.format(len(self.compounds)))

        check_missing_compound_files(self.compounds, self.folder_name)
        load_compound_energies(self.compounds, self.folder_name, parameters.minimal_intensity)

    def load(self, json_data):
        digest = json_data.get('SHA-512')
        self.digest = digest

        compound_dictionaries = json_data.get('Compounds')

        progress_label = 'Loading compounds: '
        show_progress(progress_label, 0.0)
        compound_index = 0

        for compound in compound_dictionaries:
            self.compounds.append(Compound.from_dictionary(compound))

            compound_index += 1
            show_progress(progress_label, float(compound_index) / float(len(compound_dictionaries)))
        print('Database successfully loaded ({0} compounds).'.format(len(self.compounds)))

    def save(self):
        json_data = dict()
        json_data['SHA-512'] = self.digest

        json_data['Compounds'] = []
        for compound in self.compounds:
            json_data['Compounds'].append(compound.get_dictionary())

        json_string = dumps(json_data, indent=2)
        file = open(join(self.folder_name, 'database.json'), 'w')
        file.write(json_string)
        file.close()


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
