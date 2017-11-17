from os.path import join
import subprocess

from crv.cid import CID
from crv.database import Database
from crv.experimental_compound import ExperimentalCompound
from crv.utility import show_progress


class CFMIDFilter:
    class Answer:
        class SpectrumAnswer:
            def __init__(self, id: str, match_index: int, score: float, best_alternative_name: str, best_alternative_score: float):
                self.id = id
                self.match_index = match_index
                self.score = score
                self.best_alternative_name = best_alternative_name
                self.best_alternative_score = best_alternative_score

            def __str__(self):
                return 'Spectrum answer (id={0}, match index={1}, score={2}, best alternative=\'{3}\', best alternative score={4})'.format(self.id, self.match_index, self.score, self.best_alternative_name, self.best_alternative_score)

            def __repr__(self):
                return self.__str__()

        def __init__(self, name: str):
            self.name = name
            self.spectrum_answers = []

        def add_spectrum_answer(self, id: str, match_index: int, score: float, best_alternative_name: str, best_alternative_score: float):
            self.spectrum_answers.append(CFMIDFilter.Answer.SpectrumAnswer(id, match_index, score, best_alternative_name, best_alternative_score))

        def __str__(self):
            return 'Answer (compound=\'{0}\', spectrum answers={1})'.format(self.name, self.spectrum_answers)

        def __repr__(self):
            return self.__str__()

    class ApprovedCompound:
        def __init__(self, spectrum_answer, compound: ExperimentalCompound):
            self.name = compound.name
            self.cfmid_score = spectrum_answer.score
            self.best_alternative_name = spectrum_answer.best_alternative_name
            self.best_alternative_score = spectrum_answer.best_alternative_score

            match = compound.matches[spectrum_answer.match_index]
            self.mass = match.mass
            self.rt = match.rt
            self.precursor = match.precursor
            self.area = match.area
            self.score = match.score

    def __init__(self, database: Database, session_folder_name: str):
        print('Constructing filter: CFM-ID')

        # TODO: add precision and score function as parameters
        self.command = 'cfm-id/cfm-id-precomputed.exe {0} {1} {2} -1 10 0.0005 DotProduct'

        # construct file with energies
        self.energies_file_name = join(session_folder_name, 'energies.msp')
        construct_energies_file(database.compounds, self.energies_file_name)

        # construct file with info about candidates
        self.candidates_file_name = join(session_folder_name, 'candidates.txt')
        construct_candidates_file(database.compounds, self.candidates_file_name, self.energies_file_name)

        print('CFM-ID filter constructed.')

    def apply(self, compound: ExperimentalCompound, session_folder_name: str) -> Answer:
        # if there are not matches with MS2 spectra, just return empty answer
        ms2_spectra_presented = False
        for match in compound.matches:
            if len(match.ms2_spectra) != 0:
                ms2_spectra_presented = True
                break
        if not ms2_spectra_presented:
            return CFMIDFilter.Answer(compound.name)

        # construct file with experimental energy
        experimental_energy_file_name = join(session_folder_name, '{0}.msp'.format(compound.name))
        construct_experimental_energy_file(compound, experimental_energy_file_name)

        answer = self.construct_answer(compound, experimental_energy_file_name)

        # remove all matches except match with best result?
        # remove all spectra except spectrum with best result?

        return answer

    def construct_answer(self, compound: ExperimentalCompound, energy_file_name: str) -> Answer:
        answer = CFMIDFilter.Answer(compound.name)

        # call CFM-ID and get its answer for each spectrum
        for match_index in range(0, len(compound.matches)):
            match = compound.matches[match_index]

            for spectrum in match.ms2_spectra:
                id = spectrum.cids['10'].id

                # construct command and run CFM-ID
                current_command = self.command.format(energy_file_name, id, self.candidates_file_name)
                cfm = subprocess.Popen(current_command.split(), stdout=subprocess.PIPE)
                output, error = cfm.communicate()

                # process output
                lines = output.decode('utf-8').split('\n')
                # is it correct?
                if len(lines) <= 3:
                    print()
                    print('=======================================================')
                    print('Something got wrong when calling CFM-ID.')
                    print('Command: {0}'.format(current_command))
                    print('Response:')
                    print(lines)
                    exit()

                # find line with first candidate
                first_candidate_line_index = None
                for i in range(0, len(lines)):
                    if 'ID: {0}'.format(id) in lines[i]:
                        first_candidate_line_index = i + 1
                        break
                # remove useless lines
                lines = lines[first_candidate_line_index:-1]

                # process first candidate line
                candidate_name, candidate_score = process_candidate_line(lines[0])
                # is first candidate a target compound?
                if candidate_name == compound.name:
                    # best alternative is second candidate
                    best_alternative_name, best_alternative_score = process_candidate_line(lines[1])
                    answer.add_spectrum_answer(id, match_index, candidate_score, best_alternative_name,
                                               best_alternative_score)
                else:
                    # best alternative is first candidate
                    best_alternative_name = candidate_name
                    best_alternative_score = candidate_score

                    # go through other lines looking for target compound
                    target_compound_found = False
                    for line in lines[1:]:
                        candidate_name, candidate_score = process_candidate_line(line)
                        if candidate_name == compound.name:
                            target_compound_found = True
                            answer.add_spectrum_answer(id, match_index, candidate_score, best_alternative_name,
                                                       best_alternative_score)
                            break

                    # target should not found, otherwise there is a trouble with compound naming
                    if not target_compound_found:
                        print()
                        print('============================================')
                        print('Score was not found for compound: {0}'.format(compound.name))
                        print('Possible reason: this compound has different name in experimental data and database.')
                        exit()

        return answer


def construct_energies_file(compounds: list, file_name: str):
    label = 'Constructing energies file: '
    progress = 0
    show_progress(label, progress)

    file = open(file_name, 'w')

    # write info about each compound in format that CFM-ID can understand
    for compound in compounds:
        for energy in compound.energies:
            file.write('Name: {0}\n'.format(compound.name))
            file.write('ID: {0}\n'.format(compound.name))
            file.write('Comment: {0}\n'.format(energy.name))
            file.write('Num peaks: {0}\n'.format(len(energy.fragments)))
            for fragment in energy.fragments:
                file.write('{0} {1}\n'.format(fragment.mass, fragment.intensity))
            file.write('\n')

        progress += 1
        show_progress(label, progress / len(compounds))

    file.close()


def construct_candidates_file(compounds: list, file_name: str, energies_file_name: str):
    label = 'Constructing candidates file: '
    progress = 0
    show_progress(label, progress)

    file = open(file_name, 'w')

    # write info about each compound in format that CFM-ID can understand
    for compound in compounds:
        file.write('{0} {1} {2}\n'.format(compound.name, compound.smiles, energies_file_name))

        progress += 1
        show_progress(label, progress / len(compounds))

    file.close()


def get_comment_for_cid(cid: CID):
    if cid.energy == '10':
        return 'Energy0'
    if cid.energy == '20':
        return 'Energy1'
    if cid.energy == '40':
        return 'Energy2'

    print('============================')
    print('Unsupported CID energy: {0}'.format(cid.energy))
    exit()


def construct_experimental_energy_file(compound: ExperimentalCompound, file_name: str):
    file = open(file_name, 'w')
    for match in compound.matches:
        for spectrum in match.ms2_spectra:
            id = spectrum.cids['10'].id  # same ID will be used when calling CFM-ID
            for energy in spectrum.cids:
                cid = spectrum.cids[energy]
                comment = get_comment_for_cid(cid)
                file.write('Name: {0} {1} (RT={2})\n'.format(compound.name, match.mass, match.rt))
                file.write('ID: {0}\n'.format(id))
                file.write('Comment: {0}\n'.format(comment))
                file.write('Num peaks: {0}\n'.format(len(cid.fragments)))
                for fragment in cid.fragments:
                    file.write('{0} {1}\n'.format(fragment.mass, fragment.intensity))
                file.write('\n')
    file.close()


def process_candidate_line(line: str) -> tuple:
    components = line.split(' ')
    score = float(components[1])
    name = components[2]
    return name, score
