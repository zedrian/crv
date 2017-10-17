from collections import namedtuple
from datetime import datetime
from os import listdir, makedirs, system
from os.path import isfile, join, isdir
import subprocess
from sys import stdout
from pandas import read_csv

from crv.mz_data import MzData


# smart progress bar
def show_progress(label, width, percentage):
    progress = '['
    for i in range(0, width):
        if i / width < percentage:
            progress += '#'
        else:
            progress += ' '
    progress += '] {0:.1%}'.format(percentage)
    print('\r' + label + progress, end='')
    stdout.flush()


def load_compounds_list(file_name):
    compounds_list = []
    with open(file_name) as file:
        for line in file:
            components = line[:-1].split(';')
            name = components[0]
            smile = components[1]
            compound = {'name': name, 'smile': smile}
            compounds_list.append(compound)
    print('Compounds found: {0}'.format(len(compounds_list)))

    return compounds_list


def select_database(root_folder, ask_user: bool):
    database_names = [f for f in listdir(root_folder) if isdir(join(root_folder, f))]
    index = 1
    print('Available databases:')
    for name in database_names:
        print('{0}. {1}'.format(index, name))
        index += 1

    if ask_user:
        index = int(input('Select database: '))
    else:
        index = 1

    file_name = join(root_folder, database_names[index - 1], 'compounds-list.txt')
    if not isfile(file_name):
        print('Cannot find compounds-list.txt in selected database.')
        exit()

    return database_names[index - 1]


def get_compounds_list_file_name(database_folder_name, ask_user: bool):
    if ask_user:
        file_name = input('Enter compounds list file name (or press ENTER to use default compounds-list.txt from database): ')
    else:
        file_name = ''
    if file_name == '':
        file_name = join(database_folder_name, 'compounds-list.txt')
    return file_name


def construct_database(strain_name: str, ask_user: bool = True, minimal_intensity: float = 0.0):
    database_root_folder = 'databases'
    database_name = select_database(database_root_folder, ask_user)

    database_folder_name = join(database_root_folder, database_name)
    compounds_list_file_name = get_compounds_list_file_name(database_folder_name, ask_user)

    compounds = load_compounds_list(compounds_list_file_name)

    missing_compound_files = []
    for compound in compounds:
        compound_file_name = join(database_root_folder, database_name, compound['name'])
        if not isfile(compound_file_name):
            missing_compound_files.append(compound['name'])
    if len(missing_compound_files) != 0:
        print('Cannot find energy files for {0} compounds: '.format(len(missing_compound_files)))
        for name in missing_compound_files:
            print('- {0}'.format(name))
        print('You can either provide proper energy files or delete them from compounds list.')
        exit()

    session_root_folder = 'sessions'
    time = datetime.now()
    session_name = '{0}-{1}-{2}_{3}-{4}-{5}__{6}'.format(time.year, time.month, time.day, time.hour, time.minute,
                                                         time.second, strain_name)
    session_folder_name = join(session_root_folder, session_name)
    makedirs(session_folder_name)

    msp_file_name = join(session_folder_name, 'energies.msp')
    if minimal_intensity != 0.0:
        msp_file_name = msp_file_name.replace('energies.msp', 'energies-{0}.msp'.format(minimal_intensity))
    msp_file = open(msp_file_name, 'w')
    compound_index = 0
    for compound in compounds:
        energy_file = open(join(database_folder_name, compound['name']))
        energies = []

        smile_started = False
        for line in energy_file.readlines():
            line = line[:-1]
            if line == '':
                smile_started = True
                continue
            if smile_started:
                components = line.split(' ')
                fragment_index = int(components[0])
                smile = components[2]
                for energy in energies:
                    if fragment_index < len(energy['fragments']):
                        energy['fragments'][fragment_index]['smile'] = smile
            else:
                if 'energy' in line:
                    energy = {'name': line.replace('en', 'En'), 'fragments': []}
                    energies.append(energy)
                else:
                    components = line.split(' ')
                    fragment = {'mass': float(components[0]), 'intensity': float(components[1]), 'smile': ''}
                    if minimal_intensity != 0.0:
                        if fragment['intensity'] < minimal_intensity:
                            continue
                    energies[-1]['fragments'].append(fragment)

        for energy in energies:
            msp_file.write('Name: {0}\n'.format(compound['name']))
            # msp_file.write('ID: {0}\n'.format(compound_index))
            msp_file.write('ID: {0}\n'.format(compound['name']))
            msp_file.write('Comment: {0}\n'.format(energy['name']))
            msp_file.write('Num peaks: {0}\n'.format(len(energy['fragments'])))
            for fragment in energy['fragments']:
                msp_file.write('{0} {1}\n'.format(fragment['mass'], fragment['intensity']))
            msp_file.write('\n')
        compound['energies'] = energies

        compound_index += 1
    msp_file.close()

    candidates_file_name = join(session_folder_name, 'candidates.txt')
    if minimal_intensity != 0.0:
        candidates_file_name = candidates_file_name.replace('candidates.txt',
                                                            'candidates-{0}.txt'.format(minimal_intensity))
    candidates_file = open(candidates_file_name, 'w')
    compound_index = 0
    for compound in compounds:
        # candidates_file.write('{0} {1} {2}\n'.format(compound_index, compound['smile'], msp_file_name))
        candidates_file.write('{0} {1} {2}\n'.format(compound['name'], compound['smile'], msp_file_name))
        compound_index += 1
    candidates_file.close()

    training_list_file_name = join(session_folder_name, 'training-list.txt')
    if minimal_intensity != 0.0:
        training_list_file_name = training_list_file_name.replace('training-list.txt',
                                                                  'training-list-{0}.txt'.format(minimal_intensity))
    training_list_file = open(training_list_file_name, 'w')
    training_list_file.write('{0}\n'.format(len(compounds)))
    for compound in compounds:
        training_list_file.write('{0} {1} {2}\n'.format(compound['name'], compound['smile'], 0))
    training_list_file.close()

    return session_folder_name, compounds, candidates_file_name


def filter_fragments(sample_energies, database_energies):
    if len(sample_energies) == 0:
        return sample_energies

    database_keys = ['Energy0', 'Energy1', 'Energy2']
    sample_keys = ['10', '20', '40']
    Spectrum = namedtuple('Spectrum', 'id mz rt energy masses intensities')
    Fragment = namedtuple('Fragment', 'mass intensity')
    delta_mass = 0.001

    filtered_energies = []

    methods = ['filter peaks that are close to the same from database', 'get 80% of intensity']
    method = 'get 80% of intensity'

    if method == 'filter peaks that are close to the same from database':
        for sample_energy in sample_energies:
            filtered_energy = {}
            for i in range(0, len(sample_energy)):
                sample_spectrum = sample_energy[sample_keys[i]]
                database_fragments = [e for e in database_energies if e['name'] == database_keys[i]][0]['fragments']

                nearest_database_peaks = []
                for sample_peak_index in range(0, len(sample_spectrum.masses)):
                    nearest_database_peak = database_fragments[0]
                    sample_peak = Fragment(mass=sample_spectrum.masses[sample_peak_index],
                                           intensity=sample_spectrum.intensities[sample_peak_index])
                    for database_peak in database_fragments:
                        if abs(database_peak['mass'] - sample_peak.mass) < abs(
                                        nearest_database_peak['mass'] - sample_peak.mass):
                            nearest_database_peak = database_peak
                    if abs(nearest_database_peak['mass'] - sample_peak.mass) < delta_mass:
                        nearest_database_peaks.append({'database': nearest_database_peak, 'sample': sample_peak})

                filtered_fragments = []
                for database_peak in database_fragments:
                    nearest_peaks = []
                    for peak in nearest_database_peaks:
                        if peak['database'] == database_peak:
                            nearest_peaks.append(peak['sample'])
                    if len(nearest_peaks) == 0:
                        continue

                    nearest_peak = nearest_peaks[0]
                    for peak in nearest_peaks:
                        if abs(peak.mass - database_peak['mass']) < abs(nearest_peak.mass - database_peak['mass']):
                            nearest_peak = peak
                    filtered_fragments.append(nearest_peak)

                masses = [f.mass for f in filtered_fragments]
                intensities = [f.intensity for f in filtered_fragments]
                filtered_energy[sample_keys[i]] = Spectrum(id=sample_spectrum.id, mz=sample_spectrum.mz,
                                                           rt=sample_spectrum.rt, energy=sample_spectrum.energy,
                                                           masses=masses, intensities=intensities)
            filtered_energies.append(filtered_energy)
    else:
        for sample_energy in sample_energies:
            filtered_energy = {}
            for i in range(0, len(sample_energy)):
                sample_spectrum = sample_energy[sample_keys[i]]

                sorted_fragments = []
                for fragment_index in range(0, len(sample_spectrum.masses)):
                    sorted_fragments.append(Fragment(mass=sample_spectrum.masses[fragment_index],
                                                     intensity=sample_spectrum.intensities[fragment_index]))
                sorted_fragments = sorted(sorted_fragments, key=lambda x: x.intensity, reverse=True)

                total_intensity = 0.0
                for fragment in sorted_fragments:
                    total_intensity += fragment.intensity

                target_intensity = 0.8 * total_intensity

                accumulated_intensity = 0.0
                filtered_fragments = []
                for fragment in sorted_fragments:
                    filtered_fragments.append(fragment)
                    accumulated_intensity += fragment.intensity
                    if accumulated_intensity >= target_intensity:
                        break

                masses = [f.mass for f in filtered_fragments]
                intensities = [f.intensity for f in filtered_fragments]
                filtered_energy[sample_keys[i]] = Spectrum(id=sample_spectrum.id,
                                                           mz=sample_spectrum.mz, rt=sample_spectrum.rt,
                                                           energy=sample_spectrum.energy, masses=masses,
                                                           intensities=intensities)
                filtered_energies.append(filtered_energy)

    return filtered_energies


def parse_sample_data(xml_file_name, csv_file_name, compounds_without_ms2_spectra_file_name, spectra_file_name,
                      candidates_list):
    ionization_cases_list_file_name = 'ionization-cases-list.ssv'
    ionization_cases_list = read_csv(ionization_cases_list_file_name, index_col=False, sep=';')

    Spectrum = namedtuple('Spectrum', 'id mz rt energy masses intensities')
    print('Parsing file: {0}'.format(xml_file_name))

    mz_data = MzData()
    mz_data.load(xml_file_name)

    print('Processing compounds from file: {0}'.format(csv_file_name))
    all_compounds_list = read_csv(csv_file_name, header=2, index_col=False,
                                  usecols=['Name', 'Mass', 'RT', 'Area', 'Score', 'Precursor'])
    Compound = namedtuple('Compound', 'name mass rt area score precursor matches')

    # collect compounds
    compounds_list = []
    for i in range(0, len(all_compounds_list)):
        precursor = all_compounds_list['Precursor'][i]
        mass = all_compounds_list['Mass'][i]
        matches = []
        for spectrum in mz_data.ms2_spectra:
            for ionization_case_index in range(0, len(ionization_cases_list)):
                if abs(mass + ionization_cases_list['DeltaMass'][
                    ionization_case_index] - spectrum.mz) <= mass / 100000.0:
                    if spectrum.energy == '10':
                        matches.append({'10': spectrum})
                    else:
                        matches[-1][spectrum.energy] = spectrum

        name = all_compounds_list['Name'][i]
        mass = all_compounds_list['Mass'][i]
        duplicate_rt = all_compounds_list['RT'][i]
        area = all_compounds_list['Area'][i]
        score = all_compounds_list['Score'][i]

        if precursor == mass:
            print('-- compound with precursor equal to mass found: {0}'.format(name))
            continue

        # try to filter fragments
        for candidate in candidates_list:
            if candidate['name'].lower().replace(' ', '_') == name.lower().replace(' ', '_'):
                matches = filter_fragments(matches, candidate['energies'])
                break

        compound = Compound(name=name, mass=mass, rt=duplicate_rt, area=area, score=score, precursor=precursor,
                            matches=matches)
        compounds_list.append(compound)

    # process compound duplicates
    print('Processing duplicates.')
    processed_compounds = []
    for compound in compounds_list:
        # do not process compounds that are already processed (check by name)
        if compound.name in [c.name for c in processed_compounds]:
            continue

        # get list of compounds with the same name
        duplicates = [c for c in compounds_list if c.name == compound.name]
        if len(duplicates) > 1:
            # inform user
            print('Duplicates found: {0}'.format(compound.name))
            for i in range(0, len(duplicates)):
                print('{0}. precursor={1:.3f}, mass={2:.3f}, RT={3:.3f}, area={4:.3f}, score={5:.3f}'.format(
                    i + 1, duplicates[i].precursor, duplicates[i].mass, duplicates[i].rt,
                    duplicates[i].area, duplicates[i].score))

            # collect their matches
            matches = []
            for duplicate in duplicates:
                for match in duplicate.matches:
                    matches.append(match)

            # precursor value should not be equal to mass
            for duplicate in duplicates:
                if duplicate.precursor == duplicate.mass:
                    print('Throwing away duplicate with precursor equal to mass: {0}'.format(duplicate.precursor))
                    duplicates.remove(duplicate)

            # each match should go to duplicate with closest RT
            matches_for_duplicate_by_rt = {}
            for duplicate in duplicates:
                matches_for_duplicate_by_rt[duplicate.rt] = []
            for match in matches:
                nearest_rt = duplicates[0].rt
                for duplicate_rt in matches_for_duplicate_by_rt:
                    match_rt = match['10'].rt
                    if abs(match_rt - duplicate_rt) < abs(match_rt - nearest_rt):
                        nearest_rt = duplicate_rt
                matches_for_duplicate_by_rt[nearest_rt].append(match)

            # construct new compounds with resorted matches
            for duplicate in duplicates:
                compound = Compound(name=duplicate.name, mass=duplicate.mass, rt=duplicate.rt, area=duplicate.area,
                                    score=duplicate.score, precursor=duplicate.precursor,
                                    matches=matches_for_duplicate_by_rt[duplicate.rt])
                processed_compounds.append(compound)
        else:
            processed_compounds.append(compound)
    compounds_list = processed_compounds

    # write compounds to files
    compound_without_spectra_file = open(compounds_without_ms2_spectra_file_name, 'w')
    compound_without_spectra_file.write('Name;Mass;RT;Area;Score;Precursor;PMD\n')
    spectra_file = open(spectra_file_name, 'w')
    for compound in compounds_list:
        if len(compound.matches) == 0:
            print('-- compound without spectra found: {0}'.format(compound.name))
            compound_without_spectra_file.write(
                '{0};{1};{2};{3};{4};{5};{6}\n'.format(compound.name, compound.mass, compound.rt, compound.area,
                                                       compound.score, compound.precursor,
                                                       abs(compound.precursor - compound.mass)))
        else:
            print(compound)
            for match in compound.matches:
                id = match['10'].id
                for energy in ['10', '20', '40']:
                    comment = 'Energy'
                    if energy == '10':
                        comment += '0'
                    elif energy == '20':
                        comment += '1'
                    else:
                        comment += '2'
                    if not energy in match:
                        continue
                    spectra_file.write('Name: {0} {1}\n'.format(match[energy].mz, compound.name))
                    spectra_file.write('ID: {0}\n'.format(id))
                    spectra_file.write('Comment: {0}\n'.format(comment))
                    spectra_file.write('Num peaks: {0}\n'.format(len(match[energy].masses)))
                    for index in range(0, len(match[energy].masses)):
                        spectra_file.write(
                            '{0} {1}\n'.format(match[energy].masses[index], match[energy].intensities[index]))
                    spectra_file.write('\n')
    compound_without_spectra_file.close()
    spectra_file.close()

    root = None
    tree = None

    print('Compounds list constructed.')
    print('Spectra saved to \'{0}\'.'.format(spectra_file_name))
    print('Compounds without MS2 spectra saved to \'{0}\'.'.format(compounds_without_ms2_spectra_file_name))
    return compounds_list


def receive_cfm_answers(compounds_list, candidates_list, candidates_file_name, spectra_file_name):
    label = 'CFM-ID analysis processing: '
    show_progress(label, 40, 0.0)

    command = 'cfm-id/cfm-id-precomputed.exe {0} {1} {2} -1 10 0.0005 DotProduct'

    cfm_answers = []
    progress = 0.0
    for compound in compounds_list:
        cfm_answer = {'name': compound.name, 'RT': compound.rt, 'compound_id': '', 'alternatives': {}}
        for match in compound.matches:
            id = match['10'].id
            current_command = command.format(spectra_file_name, id, candidates_file_name)
            cfm = subprocess.Popen(current_command.split(), stdout=subprocess.PIPE)
            output, error = cfm.communicate()
            lines = output.decode('utf-8').split('\n')
            first_answer_index = 0
            for i in range(0, len(lines)):
                if 'ID: ' in lines[i]:
                    first_answer_index = i + 1
                    break
            answer_lines = lines[first_answer_index:-1]
            for line in answer_lines:
                components = line.split(' ')
                candidate_name = components[2]
                candidate_score = float(components[1])
                if candidate_name.lower().replace(' ', '_') == compound.name.lower().replace(' ', '_'):
                    cfm_answer['compound_id'] = candidate_name
                if not candidate_name in cfm_answer['alternatives']:
                    cfm_answer['alternatives'][candidate_name] = []
                cfm_answer['alternatives'][candidate_name].append(candidate_score)
        if not compound.name in cfm_answer['alternatives']:
            cfm_answer['matches_ratio'] = 0.0
        else:
            cfm_answer['matches_ratio'] = len(cfm_answer['alternatives'][compound.name]) / len(compound.matches)
        cfm_answers.append(cfm_answer)
        progress += 1 / len(compounds_list)
        show_progress(label, 40, progress)
    print()
    return cfm_answers


def write_approved_compounds_list(cfm_answers, file_name, compounds_list,
                                  candidates_list):  # return dictionary: key is approved compound name, value is dictionary: key is parameter name, value is parameter value
    def get_ms2_score_for_answer(answer):
        if not answer['name'] in answer['alternatives']:
            return 0.0
        best_ms2_score = 0.0
        for score in answer['alternatives'][answer['name']]:
            best_ms2_score = max(best_ms2_score, score)
        return best_ms2_score

    # remove answers with duplicating names - leave only answers with best MS2 score
    processed_answers = []
    for answer in cfm_answers:
        # don't process twice
        if answer['name'] in [a['name'] for a in processed_answers]:
            continue

        duplicates = [a for a in cfm_answers if a['name'] == answer['name']]
        if len(duplicates) == 1:
            processed_answers.append(answer)
            continue

        # find answer with best MS2 score
        best_ms2_score = get_ms2_score_for_answer(duplicates[0])
        best_answer = duplicates[0]
        for duplicate in duplicates[1:]:
            if get_ms2_score_for_answer(duplicate) > best_ms2_score:
                best_ms2_score = get_ms2_score_for_answer(duplicate)
                best_answer = duplicate

        # store it
        processed_answers.append(best_answer)

    cfm_answers = processed_answers

    from json import loads
    database_file = open('databases/phenazines/database.json')
    json_string = ''.join(database_file.readlines())
    database_file.close()

    json_data = loads(json_string)
    fragment_matches_data = json_data['Fragment matches']

    possible_score_threshold = 0.05
    approved_score_threshold = 0.35

    approved_compounds_parameters = {}

    file = open(file_name, 'w')
    file.write(
        'Name;Mass;RT;Area;Score;Precursor;Best MS2 score;Matches ratio;Best alternative;Best alternative MS2 score;Found fragments;Unique compound fragments;Found unique fragments\n')
    for answer in cfm_answers:
        if answer['matches_ratio'] == 0.0:
            continue
        candidate_name = answer['name']
        best_ms2_score = 0.0
        for score in answer['alternatives'][candidate_name]:
            best_ms2_score = max(best_ms2_score, score)
        if best_ms2_score < possible_score_threshold:
            continue
        ms2_ratio = answer['matches_ratio']
        best_alternative_name = ''
        best_alternative_score = 0.0

        for alternative in answer['alternatives']:
            if alternative.lower().replace(' ', '_') == candidate_name.lower().replace(' ', '_'):
                continue
            current_alternative_best_score = 0.0
            for score in answer['alternatives'][alternative]:
                current_alternative_best_score = max(current_alternative_best_score, score)
            if current_alternative_best_score > best_alternative_score:
                best_alternative_name = alternative
                best_alternative_score = current_alternative_best_score
        for compound in compounds_list:
            if compound.name.lower().replace(' ', '_') == candidate_name.lower().replace(' ', '_'):
                # add to approved_compounds_parameters if not already added
                if not compound.name in approved_compounds_parameters:
                    approved_compounds_parameters[compound.name] = {}
                    approved_compounds_parameters[compound.name]['Mass'] = compound.mass
                    approved_compounds_parameters[compound.name]['RT'] = compound.rt
                    approved_compounds_parameters[compound.name]['Area'] = compound.area
                    approved_compounds_parameters[compound.name]['Score'] = compound.score
                file.write(
                    '{0};{1};{2};{3};{4};{5};{6};{7};{8};{9};'.format(compound.name, compound.mass, compound.rt,
                                                                      compound.area, compound.score,
                                                                      compound.precursor, best_ms2_score,
                                                                      ms2_ratio, best_alternative_name,
                                                                      best_alternative_score))
                found_fragment_smiles = []
                for candidate in candidates_list:
                    if candidate['name'].lower().replace(' ', '_') == candidate_name.lower().replace(' ', '_'):
                        for match in compound.matches:
                            # TODO: refactor this shit!
                            compound_keys = ['10', '20', '40']
                            for energy_index in range(0, 3):
                                if not compound_keys[energy_index] in match:
                                    continue
                                compound_energy = match[compound_keys[energy_index]]
                                if len(compound_energy) == 0:
                                    continue
                                candidate_energy = candidate['energies'][energy_index]
                                for mass in compound_energy.masses:
                                    for fragment in candidate_energy['fragments']:
                                        # 10. - is mass_tolerance variable
                                        if abs(mass - fragment['mass']) < 10. / 1000000. * max(mass,
                                                                                               fragment['mass']):
                                            smile = fragment['smile']
                                            if not smile in found_fragment_smiles:
                                                found_fragment_smiles.append(smile)
                        for smile in found_fragment_smiles:
                            file.write('{0},'.format(smile))
                        break
                file.write(';')

                # write unique compound fragments
                unique_fragments = []
                for fragment_match in fragment_matches_data:
                    containing_compounds = fragment_match['Compounds']
                    if len(containing_compounds) > 1:
                        continue
                    if containing_compounds[0] == compound.name:
                        unique_fragments.append(fragment_match['Fragment'])
                        file.write('{0},'.format(fragment_match['Fragment']))
                file.write(';')

                # write found unique fragments
                for fragment in unique_fragments:
                    if fragment in found_fragment_smiles:
                        file.write('{0},'.format(fragment))
                # file.write('')

                file.write('\n')
                break
    print('Approved compounds list saved to \'{0}\'.'.format(file_name))
    return approved_compounds_parameters


def remove_wrong_lines_from_file(all_compounds_file_name: str):
    file = open(all_compounds_file_name)
    lines = file.readlines()
    file.close()

    fixed_lines = lines[0:3]
    for line in lines[3:]:
        first_symbol = line[0]
        if not first_symbol in [str(i) for i in range(0, 10)]:
            print('Wrong line found: {0}'.format(line), end='')
        else:
            fixed_lines.append(line)

    file = open(all_compounds_file_name, 'w')
    for line in fixed_lines:
        file.write(line)
    file.close()


def folder_contains_strains(data_folder_name: str) -> bool:
    subfolders = [join(data_folder_name, f) for f in listdir(data_folder_name) if isdir(join(data_folder_name, f))]

    # get first subfolder
    subfolder = subfolders[0]

    # if subfolder contains mzdata file, it is a sample folder => data folder is a single strain data
    mzdata_files = [f for f in listdir(subfolder) if isfile(join(subfolder, f)) and '.mzdata.xml' in f]
    if len(mzdata_files) != 0:
        return False

    return True


def process_single_strain(data_folder_name: str, ask_user: bool):
    strain_name = data_folder_name.replace('\\', '/').split('/')[-1].replace(' ', '-')
    print('Processing strain: {0}'.format(strain_name))

    session_folder_name, candidates_list, candidates_file_name = construct_database(strain_name, ask_user, 0.0001)

    sample_folder_names = [f for f in listdir(data_folder_name) if isdir(join(data_folder_name, f))]
    print('Folders with sample data found:')
    for name in sample_folder_names:
        print('- {0}'.format(name))

    sample_records = []

    makedirs(join(session_folder_name, 'results'))
    for sample_folder_name in sample_folder_names:
        sample_record = {}
        sample_record['name'] = sample_folder_name

        print('Processing sample: {0}'.format(sample_folder_name))
        current_sample_data_folder_name = join(data_folder_name, sample_folder_name)
        current_sample_session_folder_name = join(session_folder_name, 'results', sample_folder_name)
        makedirs(current_sample_session_folder_name)

        compounds_without_ms2_spectra_file_name = join(current_sample_session_folder_name,
                                                       'compounds-without-ms2-spectra.ssv')
        spectra_file_name = join(current_sample_session_folder_name, 'spectra.msp')

        xml_file_name = [f for f in listdir(current_sample_data_folder_name) if '.mzdata.xml' in f][0]
        xml_file_name = join(current_sample_data_folder_name, xml_file_name)

        all_compounds_file_name = [f for f in listdir(current_sample_data_folder_name) if '.csv' in f][0]
        all_compounds_file_name = join(current_sample_data_folder_name, all_compounds_file_name)

        remove_wrong_lines_from_file(all_compounds_file_name)

        compounds_list = parse_sample_data(xml_file_name, all_compounds_file_name,
                                           compounds_without_ms2_spectra_file_name, spectra_file_name,
                                           candidates_list)


        # isotope analysis, motherfucker!
        mzdata = MzData()
        mzdata.load(xml_file_name)
        all_compounds_list = read_csv(all_compounds_file_name, header=2, index_col=False, usecols=['Name', 'Mass', 'RT'])

        # construct list of compounds
        compounds = []
        for i in range(0, len(all_compounds_list)):
            name = all_compounds_list['Name'][i]
            precursor_like_mass = all_compounds_list['Mass'][i]
            rt = all_compounds_list['RT'][i]
            compounds.append({'name': name, 'mass': precursor_like_mass, 'rt': rt})

        # load ionization cases
        ionization_cases_list_file_name = 'ionization-cases-list.ssv'
        ionization_cases_list = read_csv(ionization_cases_list_file_name, index_col=False, sep=';')

        # construct list of precursors
        for compound in compounds:
            for ionization_index in range(0, len(ionization_cases_list)):
                precursor_like_mass = compound['mass']
                delta_mass = ionization_cases_list['DeltaMass'][ionization_index]
                if not 'precursors' in compound:
                    compound['precursors'] = []
                compound['precursors'].append(precursor_like_mass + delta_mass)

        # look for MS1 spectra with the same mass
        for compound in compounds:
            # find spectrum with nearest RT
            nearest_spectrum_delta_rt = 10000
            nearest_spectrum_index = -1
            for spectrum_index in range(0, len(mzdata.ms1_spectra)):
                spectrum = mzdata.ms1_spectra[spectrum_index]
                delta_rt = abs(spectrum.rt - compound['rt'])
                if delta_rt < nearest_spectrum_delta_rt:
                    nearest_spectrum_delta_rt = delta_rt
                    nearest_spectrum_index = spectrum_index

            spectrum = mzdata.ms1_spectra[nearest_spectrum_index]

            # go through precursors trying to find almost the same mass
            for precursor_mass in compound['precursors'][::-1]:
                for mass_index in range(0, len(spectrum.masses)):
                    precursor_like_mass = spectrum.masses[mass_index]
                    if abs(precursor_like_mass - precursor_mass) <= 5e-6 * precursor_mass:
                        print('Precursor with almost same mass found:')
                        print('Compound: {0}'.format(compound['name']))
                        print('Precursor mass: {0}'.format(precursor_mass))
                        print('Spectrum masses: {0}'.format(spectrum.masses))

                        # find mass with maximal intensity
                        maximal_intensity = 0
                        maximal_intensity_mass = None
                        for i in range(0, len(spectrum.masses)):
                            precursor_like_mass = spectrum.masses[i]
                            intensity = spectrum.intensities[i]
                            if intensity > maximal_intensity:
                                maximal_intensity = intensity
                                maximal_intensity_mass = precursor_like_mass
                        print('Mass with maximal intensity: {0} ({1})'.format(maximal_intensity_mass, maximal_intensity))

                        # find intensity for precursor-like mass
                        precursor_like_mass_intensity = spectrum.intensities[mass_index]
                        print('Precursor-like mass: {0} ({1})'.format(precursor_like_mass, precursor_like_mass_intensity))

                        # find M+1
                        m_plus_1_mass = precursor_mass + 1.0009
                        print('M+1: {0}'.format(m_plus_1_mass))
                        for i in range(0, len(spectrum.masses)):
                            m_plus_1_like_mass = spectrum.masses[i]
                            if abs(m_plus_1_like_mass - m_plus_1_mass) <= 10e-6 * m_plus_1_mass:
                                m_plus_1_like_mass_intensity = spectrum.intensities[i]
                                print('M+1-like mass: {0} ({1})'.format(m_plus_1_like_mass, m_plus_1_like_mass_intensity))

                        # find M+2
                        m_plus_2_mass = precursor_mass + 1.0009 * 2
                        print('M+2: {0}'.format(m_plus_2_mass))
                        for i in range(0, len(spectrum.masses)):
                            m_plus_2_like_mass = spectrum.masses[i]
                            if abs(m_plus_2_like_mass - m_plus_2_mass) <= 10e-6 * m_plus_2_mass:
                                m_plus_2_like_mass_intensity = spectrum.intensities[i]
                                print('M+2-like mass: {0} ({1})'.format(m_plus_2_like_mass, m_plus_2_like_mass_intensity))

                        print()
                        break

        cfm_answers = receive_cfm_answers(compounds_list, candidates_list, candidates_file_name, spectra_file_name)
        print(cfm_answers)

        approved_compounds_list_file_name = join(current_sample_session_folder_name, 'approved-compounds-list.ssv')
        sample_record['approved compounds parameters'] = write_approved_compounds_list(cfm_answers,
                                                                                       approved_compounds_list_file_name,
                                                                                       compounds_list,
                                                                                       candidates_list)
        sample_records.append(sample_record)
        print('------------------------------')

    # process collected compounds parameters

    # 1. get list of all approved compounds
    all_approved_compound_names = []
    for record in sample_records:
        for compound_name in record['approved compounds parameters']:
            if not compound_name in all_approved_compound_names:
                all_approved_compound_names.append(compound_name)

    # 2. get list of collected parameters
    parameter_names = []
    for record in sample_records:
        if len(record['approved compounds parameters']) > 0:
            # small hack to get first compound name
            first_compound_name = ''
            for name in record['approved compounds parameters']:
                first_compound_name = name
                break

            for parameter_name in record['approved compounds parameters'][first_compound_name]:
                parameter_names.append(parameter_name)
            break
    parameter_to_show_plot_names = ['Area']

    # 3. get list of sample names
    sample_names = [record['name'] for record in sample_records]

    # 4. collect and show lists of parameter values for each parameter for each compound
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plot
    from math import log10

    colors = ['r', 'g', 'b', 'y', 'm', 'c']

    for parameter_name in parameter_names:
        if parameter_name in parameter_to_show_plot_names:
            total_figure = plot.figure()
            total_ax = total_figure.add_subplot(111, projection='3d')
            total_ax.set_title(parameter_name)
            total_ax.set_xlabel('Samples')
            total_ax.set_zlabel('log {0}'.format(parameter_name))
            total_ax.set_xticks([i for i in range(0, len(sample_names))])
            total_ax.set_xticklabels(sample_names)
            total_ax.set_yticks([i for i in range(0, len(all_approved_compound_names))])
            total_ax.set_yticklabels(all_approved_compound_names)

        compound_data_per_sample = {}
        compound_index = 0
        for compound_name in all_approved_compound_names:
            compound_data_per_sample[compound_name] = []
            for sample_name in sample_names:
                # get corresponding sample record
                sample_record = [r for r in sample_records if r['name'] == sample_name][0]

                # get value that corresponds to compound if it is in approved list
                if compound_name in sample_record['approved compounds parameters']:
                    parameters = sample_record['approved compounds parameters'][compound_name]
                    compound_data_per_sample[compound_name].append(parameters[parameter_name])
                else:
                    compound_data_per_sample[compound_name].append(0)

            xs = [i for i in range(0, len(sample_names))]
            ys = compound_data_per_sample[compound_name]
            log_ys = []
            for y in ys:
                if y > 0:
                    log_ys.append(log10(y))
                else:
                    log_ys.append(0)
            if parameter_name in parameter_to_show_plot_names:
                total_ax.bar(xs, log_ys, zs=compound_index, zdir='y', color=colors[compound_index % len(colors)],
                             alpha=0.8)
            compound_index += 1

        print('{0}: {1}'.format(parameter_name, compound_data_per_sample))

        if parameter_name in parameter_to_show_plot_names:
            for compound_name in all_approved_compound_names:
                per_compound_figure = plot.figure()
                per_compound_ax = per_compound_figure.add_subplot(111)
                per_compound_ax.set_title(compound_name)
                per_compound_ax.set_xlabel('Samples')
                per_compound_ax.set_ylabel(parameter_name)
                per_compound_ax.set_xticks([i for i in range(0, len(sample_names))])
                per_compound_ax.set_xticklabels(sample_names)
                xs = [i for i in range(0, len(sample_names))]
                ys = compound_data_per_sample[compound_name]
                per_compound_ax.bar(xs, ys, color='b', alpha=1.)
                plot_file_name = join(session_folder_name, 'results',
                                      '{0}_{1}.png'.format(parameter_name, compound_name))
                per_compound_figure.savefig(plot_file_name, dpi=600, facecolor='w', orientation='landscape',
                                            format='png')

        parameter_file_name = join(session_folder_name, 'results', '{0}.ssv'.format(parameter_name))
        file = open(parameter_file_name, 'w')
        file.write('Compound')
        for sample_name in sample_names:
            file.write(';{0}'.format(sample_name))
        file.write('\n')
        for compound_name in all_approved_compound_names:
            file.write(compound_name)
            for parameter_value in compound_data_per_sample[compound_name]:
                file.write(';{0}'.format(parameter_value))
            file.write('\n')
        file.close()
        print('{0} data saved to {1}.'.format(parameter_name, parameter_file_name))
    print('==============================================================================')


def main():
    data_folder_name = input('Enter data folder name: ')

    process_several_strains = folder_contains_strains(data_folder_name)

    if not process_several_strains:
        print('Data folder contains only one strain.')
        process_single_strain(data_folder_name, True)
    else:
        print('Data folder contains several strains.')
        strain_folder_names = [join(data_folder_name, f) for f in listdir(data_folder_name) if isdir(join(data_folder_name, f))]
        for strain_folder_name in strain_folder_names:
            process_single_strain(strain_folder_name, False)


if __name__ == '__main__':
    main()
