from os.path import isfile
from json import dumps, loads
from json.decoder import JSONDecodeError

from crv.utility import value_if_not_none


# Class that should hold all parameters needed for current CRV session (one execution) such as mass tolerance,
# database name and others.
# Uses JSON as storage format.
# You can look for default values in constructor as arguments for `value_if_not_none()` calls.
# If you want to add new parameters you should edit constructor, `from_json_file()` and `get_parameters_dictionary()`.
class Parameters:
    def __init__(self, mass_tolerance_ppm: float = None, scoring_function: float = None, minimal_score: float = None,
                 approved_score: float = None, database_name: str = None, database_compounds_file_name: str = None,
                 experimental_data_folder_name: str = None, experimental_compounds_file_name: str = None,
                 experimental_mzdata_file_name: str = None, generated_energies_file_name: str = None,
                 generated_candidates_file_name: str = None, generated_training_list_file_name: str = None,
                 minimal_intensity: float = None, postprocess_predicted_spectra: bool = None, cfm_id_type: str = None,
                 ionization_mode: str = None, ions_file_name: str = None, positive_ions: list = None,
                 negative_ions: list = None, neutral_losses: list = None):
        self.mass_tolerance_ppm = value_if_not_none(mass_tolerance_ppm, 10.0)
        self.scoring_function = value_if_not_none(scoring_function, 'DotProduct')
        self.minimal_score = value_if_not_none(minimal_score, 0.1)
        self.approved_score = value_if_not_none(approved_score, 0.4)
        self.database_name = value_if_not_none(database_name, '')
        self.database_compounds_file_name = value_if_not_none(database_compounds_file_name, 'compounds-list.txt')
        self.experimental_data_folder_name = value_if_not_none(experimental_data_folder_name, '')
        self.experimental_compounds_file_name = value_if_not_none(experimental_compounds_file_name, 'all compounds.csv')
        self.experimental_mzdata_file_name = value_if_not_none(experimental_mzdata_file_name, 'mzdata.xml')
        self.generated_energies_file_name = value_if_not_none(generated_energies_file_name, 'energies.msp')
        self.generated_candidates_file_name = value_if_not_none(generated_candidates_file_name, 'candidates.txt')
        self.generated_training_list_file_name = value_if_not_none(generated_training_list_file_name,
                                                                   'training-list.txt')
        self.minimal_intensity = value_if_not_none(minimal_intensity, 0.001)
        self.postprocess_predicted_spectra = value_if_not_none(postprocess_predicted_spectra, False)
        self.cfm_id_type = value_if_not_none(cfm_id_type, 'cfm-id-precomputed')
        self.ionization_mode = value_if_not_none(ionization_mode, 'ESI MS/MS')
        self.ions_file_name = value_if_not_none(ions_file_name, 'ions.txt')
        self.positive_ions = value_if_not_none(positive_ions, ['+H', '+NH4'])
        self.negative_ions = value_if_not_none(negative_ions, ['-H', '+CL'])
        self.neutral_losses = value_if_not_none(neutral_losses, ['H2O'])

    @classmethod
    def from_json_file(cls, file_name: str = ''):
        if isfile(file_name):
            file = open(file_name)
            json_string = ''.join(file.readlines())
            file.close()

            parameters_dictionary = dict()
            try:
                parameters_dictionary = loads(json_string)
            except JSONDecodeError:
                print('Failed to load parameters from file \'{0}\' as JSON.'.format(file_name))
                print('If you have troubles with constructing parameters file, try to use default.')
                exit()

            mass_tolerance_ppm = parameters_dictionary.get('Mass tolerance (PPM)')
            scoring_function = parameters_dictionary.get('Scoring function (\'DotProduct\' / \'Jaccard\')')
            minimal_score = parameters_dictionary.get('Minimal score')
            approved_score = parameters_dictionary.get('Approved score')
            database_name = parameters_dictionary.get('Database')
            database_compounds_file_name = parameters_dictionary.get('Database compounds file')
            experimental_data_folder_name = parameters_dictionary.get('Experimental data folder')
            experimental_compounds_file_name = parameters_dictionary.get('Experimental compounds file')
            experimental_mzdata_file_name = parameters_dictionary.get('Experimental .mzdata file')
            generated_energies_file_name = parameters_dictionary.get('Generated energies .msp file')
            generated_candidates_file_name = parameters_dictionary.get('Generated candidates file')
            generated_training_list_file_name = parameters_dictionary.get('Generated training list file')
            minimal_intensity = parameters_dictionary.get('Minimal intensity')
            postprocess_predicted_spectra = parameters_dictionary.get('Post-process predicted spectra')
            cfm_id_type = parameters_dictionary.get('CFM-ID type (\'cfm-id\' / \'cfm-id-precomputed\')')
            ionization_mode = parameters_dictionary.get('Ionization mode (???)')
            ions_file_name = parameters_dictionary.get('Ions file')
            positive_ions = parameters_dictionary.get('Positive ions')
            negative_ions = parameters_dictionary.get('Negative ions')
            neutral_losses = parameters_dictionary.get('Neutral losses')
            return Parameters(mass_tolerance_ppm, scoring_function, minimal_score, approved_score, database_name,
                              database_compounds_file_name, experimental_data_folder_name,
                              experimental_compounds_file_name, experimental_mzdata_file_name,
                              generated_energies_file_name, generated_candidates_file_name,
                              generated_training_list_file_name, minimal_intensity, postprocess_predicted_spectra,
                              cfm_id_type, ionization_mode, ions_file_name, positive_ions, negative_ions,
                              neutral_losses)
        else:
            print('Failed to load parameters from file \'{0}\': file not found.'.format(file_name))
            exit()

    def get_parameters_dictionary(self) -> dict():
        parameters_dictionary = dict()
        parameters_dictionary['Mass tolerance (PPM)'] = self.mass_tolerance_ppm
        parameters_dictionary['Scoring function (\'DotProduct\' / \'Jaccard\')'] = self.scoring_function
        parameters_dictionary['Minimal score'] = self.minimal_score
        parameters_dictionary['Approved score'] = self.approved_score
        parameters_dictionary['Database'] = self.database_name
        parameters_dictionary['Database compounds file'] = self.database_compounds_file_name
        parameters_dictionary['Experimental data folder'] = self.experimental_data_folder_name
        parameters_dictionary['Experimental compounds file'] = self.experimental_compounds_file_name
        parameters_dictionary['Experimental .mzdata file'] = self.experimental_mzdata_file_name
        parameters_dictionary['Generated energies .msp file'] = self.generated_energies_file_name
        parameters_dictionary['Generated candidates file'] = self.generated_candidates_file_name
        parameters_dictionary['Generated training list file'] = self.generated_training_list_file_name
        parameters_dictionary['Minimal intensity'] = self.minimal_intensity
        parameters_dictionary['Post-process predicted spectra'] = self.postprocess_predicted_spectra
        parameters_dictionary['CFM-ID type (\'cfm-id\' / \'cfm-id-precomputed\')'] = self.cfm_id_type
        parameters_dictionary['Ionization mode (???)'] = self.ionization_mode
        parameters_dictionary['Ions file'] = self.ions_file_name
        parameters_dictionary['Positive ions'] = self.positive_ions
        parameters_dictionary['Negative ions'] = self.negative_ions
        parameters_dictionary['Neutral losses'] = self.neutral_losses
        return parameters_dictionary

    def __str__(self) -> str:
        parameters_dictionary = self.get_parameters_dictionary()
        parameter_strings = ['{0}: {1}'.format(key, parameters_dictionary[key]) for key in parameters_dictionary]
        return '\n'.join(parameter_strings)

    def write_to_json_file(self, file_name):
        file = open(file_name, 'w')
        file.write(dumps(self.get_parameters_dictionary(), indent=2))
        file.close()


def generate_default_parameters_file():
    file_name = 'parameters.json'

    default_parameters = Parameters()
    default_parameters.write_to_json_file(file_name)
    print('Default parameters file {0} generated.'.format(file_name))
