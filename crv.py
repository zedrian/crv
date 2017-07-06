from os import listdir, makedirs
from os.path import isfile, join, isdir
from datetime import datetime
from pandas import DataFrame, read_csv
from xml.etree.ElementTree import parse
from base64 import standard_b64decode
from struct import unpack
from collections import namedtuple
import subprocess
from sys import stdout


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


def select_database(root_folder):
    database_names = [f for f in listdir(root_folder) if isdir(join(root_folder, f))]
    index = 1
    print('Available databases:')
    for name in database_names:
        print('{0}. {1}'.format(index, name))
        index += 1
    index = int(input('Select database: '))

    file_name = join(root_folder, database_names[index - 1], 'compounds-list.txt')
    if not isfile(file_name):
        print('Cannot find compounds-list.txt in selected database.')
        exit()

    return database_names[index-1]


def get_compounds_list_file_name(database_folder_name):
    file_name = input('Enter compounds list file name (or press ENTER to use default compounds-list.txt from database): ')
    if file_name == '':
        file_name = join(database_folder_name, 'compounds-list.txt')
    return file_name


def construct_database(minimal_intensity: float = 0.0):
    database_root_folder = 'databases'
    database_name = select_database(database_root_folder)

    database_folder_name = join(database_root_folder, database_name)
    compounds_list_file_name = get_compounds_list_file_name(database_folder_name)

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
    session_name = '{0}-{1}-{2}_{3}-{4}-{5}'.format(time.year, time.month, time.day, time.hour, time.minute, time.second)
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
        candidates_file_name = candidates_file_name.replace('candidates.txt', 'candidates-{0}.txt'.format(minimal_intensity))
    candidates_file = open(candidates_file_name, 'w')
    compound_index = 0
    for compound in compounds:
        # candidates_file.write('{0} {1} {2}\n'.format(compound_index, compound['smile'], msp_file_name))
        candidates_file.write('{0} {1} {2}\n'.format(compound['name'], compound['smile'], msp_file_name))
        compound_index += 1
    candidates_file.close()

    training_list_file_name = join(session_folder_name, 'training-list.txt')
    if minimal_intensity != 0.0:
        training_list_file_name = training_list_file_name.replace('training-list.txt', 'training-list-{0}.txt'.format(minimal_intensity))
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
    Spectrum = namedtuple('Spectrum', 'id mz energy masses intensities')
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
                    sample_peak = Fragment(mass=sample_spectrum.masses[sample_peak_index], intensity=sample_spectrum.intensities[sample_peak_index])
                    for database_peak in database_fragments:
                        if abs(database_peak['mass'] - sample_peak.mass) < abs(nearest_database_peak['mass'] - sample_peak.mass):
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
                filtered_energy[sample_keys[i]] = Spectrum(id=sample_spectrum.id, mz=sample_spectrum.mz, energy=sample_spectrum.energy, masses=masses, intensities=intensities)
            filtered_energies.append(filtered_energy)
    else:
        for sample_energy in sample_energies:
            filtered_energy = {}
            for i in range(0, len(sample_energy)):
                sample_spectrum = sample_energy[sample_keys[i]]

                sorted_fragments = []
                for fragment_index in range(0, len(sample_spectrum.masses)):
                    sorted_fragments.append(Fragment(mass=sample_spectrum.masses[fragment_index], intensity=sample_spectrum.intensities[fragment_index]))
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
                filtered_energy[sample_keys[i]] = Spectrum(id=sample_spectrum.id, mz=sample_spectrum.mz,
                                                           energy=sample_spectrum.energy, masses=masses,
                                                           intensities=intensities)
            filtered_energies.append(filtered_energy)

    return filtered_energies


def parse_sample_data(xml_file_name, csv_file_name, compounds_without_ms2_spectra_file_name, spectra_file_name, candidates_list):
    ionization_cases_list_file_name = 'ionization-cases-list.ssv'
    ionization_cases_list = read_csv(ionization_cases_list_file_name, index_col=False, sep=';')

    Spectrum = namedtuple('Spectrum', 'id mz energy masses intensities')
    print('Parsing file: {0}'.format(xml_file_name))
    tree = parse(xml_file_name)
    root = tree.getroot()
    spectrumList = root.find('spectrumList')

    spectrum_list = spectrumList.findall('spectrum')
    print('Spectra found: {0}'.format(len(spectrum_list)))

    spectrum_with_precursors_list = []
    for spectrum in spectrum_list:
        desc = spectrum.find('spectrumDesc')
        if desc.find('precursorList') is not None:
            spectrum_with_precursors_list.append(spectrum)
    print('Spectra with precursors found: {0}'.format(len(spectrum_with_precursors_list)))

    spectrum_data_list = []
    for spectrum in spectrum_with_precursors_list:
        id = spectrum.attrib['id']
        precursorList = spectrum.find('spectrumDesc').find('precursorList')
        precursor_list = precursorList.findall('precursor')
        if len(precursor_list) > 1:
            print('WARNING: spectrum with {0} precursors found: {1}'.format(len(precursor_list), id))
        for precursor in precursor_list:
            mz = 0.0
            cvParam_list = precursor.find('ionSelection').findall('cvParam')
            for cvParam in cvParam_list:
                if cvParam.attrib['name'] == 'MassToChargeRatio':
                    mz = float(cvParam.attrib['value'])
            energy = 0
            cvParam_list = precursor.find('activation').findall('cvParam')
            for cvParam in cvParam_list:
                if cvParam.attrib['name'] == 'CollisionEnergy':
                    energy = cvParam.attrib['value']

        mzArrayBinary = spectrum.find('mzArrayBinary').find('data')
        mzArrayBinary_data = mzArrayBinary.text
        mzArrayBinary_decoded = standard_b64decode(mzArrayBinary_data)
        mzArrayBinary_count = len(mzArrayBinary_decoded) // 8
        mzArray = unpack('<{0}d'.format(mzArrayBinary_count), mzArrayBinary_decoded)

        intenArrayBinary = spectrum.find('intenArrayBinary').find('data')
        intenArrayBinary_data = intenArrayBinary.text
        intenArrayBinary_decoded = standard_b64decode(intenArrayBinary_data)
        intenArray = unpack('<{0}f'.format(mzArrayBinary_count), intenArrayBinary_decoded)

        masses = []
        intensities = []
        for i in range(0, len(mzArray)):
            if mzArray[i] < 1.00001 * mz:
                masses.append(mzArray[i])
                intensities.append(intenArray[i])

        spectrum_data_list.append(Spectrum(id=id, mz=mz, energy=energy, masses=masses, intensities=intensities))

    print('Processing compounds from file: {0}'.format(csv_file_name))
    all_compounds_list = read_csv(csv_file_name, header=2, index_col=False, usecols=['Name', 'Mass', 'RT', 'Area', 'Score', 'Precursor'])
    Compound = namedtuple('Compound', 'name mass rt area score precursor matches')

    compound_without_spectra_file = open(compounds_without_ms2_spectra_file_name, 'w')
    compound_without_spectra_file.write('Name;Mass;RT;Area;Score;Precursor;PMD\n')

    spectra_file = open(spectra_file_name, 'w')
    compounds_list = []
    for i in range(0, len(all_compounds_list)):
        precursor = all_compounds_list['Precursor'][i]
        mass = all_compounds_list['Mass'][i]
        matches = []
        for spectrum in spectrum_data_list:
            for ionization_case_index in range(0, len(ionization_cases_list)):
                if abs(mass + ionization_cases_list['DeltaMass'][ionization_case_index] - spectrum.mz) <= mass / 100000.0:
                    if spectrum.energy == '10':
                        matches.append({'10': spectrum})
                    else:
                        matches[-1][spectrum.energy] = spectrum

        name = all_compounds_list['Name'][i]
        mass = all_compounds_list['Mass'][i]
        rt = all_compounds_list['RT'][i]
        area = all_compounds_list['Area'][i]
        score = all_compounds_list['Score'][i]

        # try to filter fragments
        for candidate in candidates_list:
            if candidate['name'] == name:
                matches = filter_fragments(matches, candidate['energies'])

        compound = Compound(name=name, mass=mass, rt=rt, area=area, score=score, precursor=precursor, matches=matches)
        if len(matches) == 0:
            print('-- compound without spectra found: {0}'.format(name))
            compound_without_spectra_file.write('{0};{1};{2};{3};{4};{5};{6}\n'.format(name, mass, rt, area, score, precursor, abs(precursor-mass)))
        else:
            print(compound)
            for match in compound.matches:
                id = match['10'].id
                for energy in ['10', '20', '40']:
                    comment = 'Energy'
                    if energy == '10': comment += '0'
                    elif energy == '20': comment += '1'
                    else: comment += '2'
                    if not energy in match:
                        continue
                    spectra_file.write('Name: {0} {1}\n'.format(match[energy].mz, compound.name))
                    spectra_file.write('ID: {0}\n'.format(id))
                    spectra_file.write('Comment: {0}\n'.format(comment))
                    spectra_file.write('Num peaks: {0}\n'.format(len(match[energy].masses)))
                    for index in range(0, len(match[energy].masses)):
                        spectra_file.write('{0} {1}\n'.format(match[energy].masses[index], match[energy].intensities[index]))
                    spectra_file.write('\n')
            compounds_list.append(compound)
    compound_without_spectra_file.close()
    spectra_file.close()
    print('Compounds list constructed.')
    print('Spectra saved to \'{0}\'.'.format(spectra_file_name))
    print('Compounds without MS2 spectra saved to \'{0}\'.'.format(compounds_without_ms2_spectra_file_name))
    return compounds_list


def receive_cfm_answers(compounds_list, candidates_list, candidates_file_name, spectra_file_name):
    label = 'CFM-ID analysis processing: '
    show_progress(label, 40, 0.0)

    command = 'cfm-id/cfm-id-precomputed.exe {0} {1} {2} -1 10 0.0005'# DotProduct'

    cfm_answers = []
    progress = 0.0
    for compound in compounds_list:
        cfm_answer = {'name': compound.name, 'compound_id': '', 'alternatives': {}}
        for match in compound.matches:
            id = match['10'].id
            current_command = command.format(spectra_file_name, id, candidates_file_name)
            cfm = subprocess.Popen(current_command.split(), stdout=subprocess.PIPE)
            output, error = cfm.communicate()
            lines = output.decode('utf-8').split('\n')
            first_answer_index = 0
            for i in range(0, len(lines)):
                if 'ID: ' in lines[i]:
                    first_answer_index = i+1
                    break
            answer_lines = lines[first_answer_index:-1]
            for line in answer_lines:
                components = line.split(' ')
                candidate_name = components[2]
                candidate_score = float(components[1])
                if candidate_name == compound.name:
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


def write_approved_compounds_list(cfm_answers, file_name, compounds_list, candidates_list):
    possible_score_threshold = 0.05
    approved_score_threshold = 0.35

    file = open(file_name, 'w')
    file.write('Name;Mass;RT;Area;Score;Precursor;Average MS2 score;Matches ratio;Nearest alternative;Nearest alternative MS2 score;Fragments\n')
    for answer in cfm_answers:
        if answer['matches_ratio'] == 0.0:
            continue
        candidate_name = answer['name']
        average_ms2_score = 0.0
        for score in answer['alternatives'][candidate_name]:
            average_ms2_score += score
        average_ms2_score /= len(answer['alternatives'][candidate_name])
        if average_ms2_score < possible_score_threshold:
            continue
        ms2_ratio = answer['matches_ratio']
        nearest_alternative_name = ''
        nearest_alternative_score = 0.0
        for alternative in answer['alternatives']:
            if alternative == candidate_name:
                continue
            current_alternative_average_score = 0.0
            for score in answer['alternatives'][alternative]:
                current_alternative_average_score += score
            current_alternative_average_score /= len(answer['alternatives'][alternative])
            if abs(average_ms2_score - current_alternative_average_score) < abs(average_ms2_score - nearest_alternative_score):
                nearest_alternative_name = alternative
                nearest_alternative_score = current_alternative_average_score
        for compound in compounds_list:
            if compound.name == candidate_name:
                file.write('{0};{1};{2};{3};{4};{5};{6};{7};{8};{9}'.format(compound.name, compound.mass, compound.rt,
                                                                       compound.area, compound.score,
                                                                       compound.precursor, average_ms2_score,
                                                                       ms2_ratio, nearest_alternative_name, nearest_alternative_score))
                if average_ms2_score > approved_score_threshold:
                    file.write(';\n')
                else:
                    for candidate in candidates_list:
                        if candidate['name'] == candidate_name:
                            corresponding_smiles = []
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
                                            if abs(mass - fragment['mass']) < 10.00:
                                                smile = fragment['smile']
                                                if not smile in corresponding_smiles:
                                                    corresponding_smiles.append(smile)
                            for smile in corresponding_smiles:
                                file.write(';{0}'.format(smile))
                            break
                    file.write('\n')
                break
    print('Approved compounds list saved to \'{0}\'.'.format(file_name))


if __name__ == '__main__':
    session_folder_name, candidates_list, candidates_file_name = construct_database(0.0001)

    data_folder_name = input('Enter data folder name: ')
    sample_folder_names = [f for f in listdir(data_folder_name) if isdir(join(data_folder_name, f))]
    print('Folders with sample data found:')
    for name in sample_folder_names:
        print('- {0}'.format(name))

    makedirs(join(session_folder_name, 'results'))
    for sample_folder_name in sample_folder_names:
        print('Processing sample: {0}'.format(sample_folder_name))
        current_sample_data_folder_name = join(data_folder_name, sample_folder_name)
        current_sample_session_folder_name = join(session_folder_name, 'results', sample_folder_name)
        makedirs(current_sample_session_folder_name)

        compounds_without_ms2_spectra_file_name = join(current_sample_session_folder_name, 'compounds-without-ms2-spectra.ssv')
        spectra_file_name = join(current_sample_session_folder_name, 'spectra.msp')

        xml_file_name = [f for f in listdir(current_sample_data_folder_name) if '.mzdata.xml' in f][0]
        xml_file_name = join(current_sample_data_folder_name, xml_file_name)

        all_compounds_file_name = [f for f in listdir(current_sample_data_folder_name) if '.csv' in f][0]
        all_compounds_file_name = join(current_sample_data_folder_name, all_compounds_file_name)
        compounds_list = parse_sample_data(xml_file_name, all_compounds_file_name, compounds_without_ms2_spectra_file_name, spectra_file_name, candidates_list)

        cfm_answers = receive_cfm_answers(compounds_list, candidates_list, candidates_file_name, spectra_file_name)
        print(cfm_answers)

        approved_compounds_list_file_name = join(current_sample_session_folder_name, 'approved-compounds-list.ssv')
        write_approved_compounds_list(cfm_answers, approved_compounds_list_file_name, compounds_list, candidates_list)
        print('------------------------------')
