from datetime import datetime
from os import makedirs
from os.path import isfile, join
from sys import argv

from crv.compound import Compound
from crv.database import Database
from crv.energy import Energy
from crv.fragment import Fragment
from crv.parameters import Parameters, generate_default_parameters_file
from crv.utility import show_progress


def show_usage():
    script_file_name = argv[0].split('/')[-1]
    print('Usage: {0} parameters_json_file [mode]'.format(script_file_name))
    print('If parameters_json_file is not presented, CRV automatically generates \'parameters.json\' with default '
          'values.')
    print('If mode is not presented, CRV shows list of available modes and user should choose one among them or exit.')
    print('')


# Generates folder for current session in next format: sessions/YYYY-MM-DD_HH-MM-SS
def generate_session_folder() -> str:
    root_folder = 'sessions'
    time = datetime.now()
    session_name = '{0}-{1}-{2}_{3}-{4}-{5}'.format(time.year, time.month, time.day, time.hour, time.minute,
                                                    time.second)
    session_folder_name = join(root_folder, session_name)
    makedirs(session_folder_name)

    return session_folder_name


# Loads compound names and SMILES from file and returns list of compounds.
def load_compounds(file_name: str) -> list:
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


# Checks energy files existence for each compound from list.
# If any, shows list of missing compound files and exits.
def check_missing_compound_files(compounds: list, database_folder_name: str):
    missing_files = []

    for compound in compounds:
        file_name = join(database_folder_name, compound.name)
        if not isfile(file_name):
            missing_files.append(compound.name)

    if len(missing_files) != 0:
        print('Cannot find energy files for {0} compounds:'.format(len(missing_files)))
        for name in missing_files:
            print('- {0}'.format(name))
        print('You can either provide proper energy files or delete them from compounds list.')
        exit()

    print('All needed energy files exist in database.')


def load_compound_energies(compounds: list, database_folder_name: str, minimal_intensity: float):
    progress_label = 'Loading compound energies: '
    show_progress(progress_label, 0.0)

    compound_index = 0
    for compound in compounds:
        energy_file = open(join(database_folder_name, compound.name))
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


def write_compound_energies_to_msp(compounds: list, energies_file_name: str):
    progress_label = 'Writing compound energies to .msp: '
    show_progress(progress_label, 0.0)

    file = open(energies_file_name, 'w')
    compound_index = 0
    for compound in compounds:
        for energy in compound.energies:
            file.write('Name: {0}\n'.format(compound.name))
            file.write('ID: {0}\n'.format(compound.name))
            file.write('Comment: {0}\n'.format(energy.name))
            file.write('Num peaks: {0}\n'.format(len(energy.fragments)))

            for fragment in energy.fragments:
                file.write('{0} {1}\n'.format(fragment.mass, fragment.intensity))
            file.write('\n')

        compound_index += 1
        show_progress(progress_label, float(compound_index) / float(len(compounds)))
    file.close()


def write_candidates_list(compounds: list, candidates_file_name: str, energies_file_name: str):
    progress_label = 'Writing candidates list: '
    show_progress(progress_label, 0.0)

    file = open(candidates_file_name, 'w')
    compound_index = 0
    for compound in compounds:
        file.write('{0} {1} {2}\n'.format(compound.name, compound.smiles, energies_file_name))

        compound_index += 1
        show_progress(progress_label, float(compound_index) / float(len(compounds)))
    file.close()


def write_training_list(compounds: list, training_list_file_name: str):
    progress_label = 'Writing training list: '
    show_progress(progress_label, 0.0)

    file = open(training_list_file_name, 'w')
    file.write('{0}\n'.format(len(compounds)))

    compound_index = 0
    for compound in compounds:
        file.write('{0} {1} {2}\n'.format(compound.name, compound.smiles, 0))

        compound_index += 1
        show_progress(progress_label, float(compound_index) / float(len(compounds)))
    file.close()


# Constructs compounds list and all needed files from database for current session.
# Returns list of compounds.
def construct_database(parameters: Parameters, session_folder_name: str) -> list:
    root_folder = 'databases'
    database_folder_name = join(root_folder, parameters.database_name)

    compounds_file_name = join(database_folder_name, parameters.database_compounds_file_name)
    compounds = load_compounds(compounds_file_name)

    check_missing_compound_files(compounds, database_folder_name)
    load_compound_energies(compounds, database_folder_name, parameters.minimal_intensity)

    energies_file_name = join(session_folder_name, parameters.generated_energies_file_name)
    write_compound_energies_to_msp(compounds, energies_file_name)

    candidates_file_name = join(session_folder_name, parameters.generated_candidates_file_name)
    write_candidates_list(compounds, candidates_file_name, energies_file_name)

    training_list_file_name = join(session_folder_name, parameters.generated_training_list_file_name)
    write_training_list(compounds, training_list_file_name)


if __name__ == '__main__':
    available_modes = ['experiment', 'statistics', 'train']
    parameters = None
    mode = None
    if len(argv) == 1:
        # if script is running without parameters, show usage and generate file with default values
        show_usage()
        generate_default_parameters_file()
        exit()
    else:
        parameters = Parameters.from_json_file(argv[1])
        if len(argv) > 2:
            mode = argv[2]          
            if not mode in available_modes:
                print('Mode is not supported: {0}.'.format(mode))
            else:
                print('Mode: {0}'.format(mode))
        else:
            print('No mode provided from command line.')
    if not mode in available_modes:
        print('Select mode:')
    

    # print('Current parameters:')
    # print(parameters)
    # print('')

    # session_folder_name = generate_session_folder()

    database = Database(parameters)
