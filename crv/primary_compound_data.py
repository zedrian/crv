from pandas import read_csv


class PrimaryCompoundData:
    def __init__(self):
        # self.compounds = []
        self.data = None

    def load(self, file_name: str):
        print('Loading primary compound data from file: {0}'.format(file_name))

        # check at first are there non-data lines in the head of file
        file = open(file_name)
        lines = file.readlines()
        file.close()
        if lines[0].startswith('#'):
            # file needs processing

            # remove first two lines
            lines = lines[2:]

            # split lines into components by comma
            components = [line.split(',') for line in lines]

            # remove last 10 components from each line
            components = [c[:-10] for c in components]

            # join them back to lines
            fixed_lines = [','.join(c) for c in components]

            # write to file with .fixed postfix
            fixed_file_name = file_name.replace('.csv', '.fixed.csv')
            file = open(fixed_file_name, 'w')
            # write header first
            file.write('{0}\n'.format(fixed_lines[0]))
            # write rest lines if they start from number
            for line in fixed_lines[1:]:
                if line == '':
                    continue  # no need in empty lines
                if line[0] in '0123456789':
                    file.write('{0}\n'.format(line))
            file.close()
            file_name = fixed_file_name

        # now we should be sure that file is correct, read it
        self.data = read_csv(file_name, index_col=False, usecols=['Name', 'Mass', 'RT', 'Area', 'Score', 'Precursor'])
