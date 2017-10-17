from base64 import standard_b64decode
from struct import unpack
from xml.etree.ElementTree import parse, ElementTree

from crv.ms1_spectrum import MS1Spectrum
from crv.ms2_spectrum import MS2Spectrum


class MzData:
    def __init__(self):
        self.file_name = None
        self.ms1_spectra = []
        self.ms2_spectra = []

    def load(self, file_name: str):
        print('Extracting MZ data from file: {0}'.format(file_name))
        self.file_name = file_name
        tree = parse(file_name)
        root = tree.getroot()

        spectrumList_node = root.find('spectrumList')
        spectra_nodes = spectrumList_node.findall('spectrum')
        print('Spectra found: {0}'.format(len(spectra_nodes)))

        for node in spectra_nodes:
            self.load_spectrum_from_xml(node)
        print('MS1 spectra loaded: {0}'.format(len(self.ms1_spectra)))
        print('MS2 spectra loaded: {0}'.format(len(self.ms2_spectra)))

    def load_spectrum_from_xml(self, xml: ElementTree):
        id = xml.attrib['id']

        spectrumInstrument_node = xml.find('spectrumDesc').find('spectrumSettings').find('spectrumInstrument')
        level = spectrumInstrument_node.attrib['msLevel']  # '1' or '2'

        if level == '1':
            # load mz
            mz_min = float(spectrumInstrument_node.attrib['mzRangeStart'])
            mz_max = float(spectrumInstrument_node.attrib['mzRangeStop'])

            # load RT
            rt = None
            cvParam_nodes = spectrumInstrument_node.findall('cvParam')
            for cvParam_node in cvParam_nodes:
                if cvParam_node.attrib['name'] == 'TimeInMinutes':
                    rt_attribute = cvParam_node.attrib['value']
                    # check if RT attribute consists of two components with '-' between them
                    if '-' in rt_attribute:
                        # RT is presented by two values: min-max, get median
                        components = rt_attribute.split('-')
                        rt_min = float(components[0])
                        rt_max = float(components[1])
                        rt = (rt_min + rt_max) / 2.
                    else:
                        # RT is a single value
                        rt = float(rt_attribute)
                    break

            # load masses
            mzArrayBinary_node = xml.find('mzArrayBinary').find('data')
            mz_data = mzArrayBinary_node.text
            if mz_data is None:
                print('WARNING: MS1 spectrum without mz data found: {0}'.format(id))
                return
            mz_decoded = standard_b64decode(mz_data)
            mz_count = len(mz_decoded) // 8
            masses = unpack('<{0}d'.format(mz_count), mz_decoded)

            # load intensities
            intenArrayBinary_node = xml.find('intenArrayBinary').find('data')
            intensities_data = intenArrayBinary_node.text
            intensities_decoded = standard_b64decode(intensities_data)
            intensities = unpack('<{0}f'.format(mz_count), intensities_decoded)

            self.ms1_spectra.append(MS1Spectrum(id, mz_min, mz_max, rt, masses, intensities))
        else:
            # load RT
            rt = None
            cvParam_nodes = spectrumInstrument_node.findall('cvParam')
            for cvParam_node in cvParam_nodes:
                if cvParam_node.attrib['name'] == 'TimeInMinutes':
                    rt_attribute = cvParam_node.attrib['value']
                    # check if RT attribute consists of two components with '-' between them
                    if '-' in rt_attribute:
                        # RT is presented by two values: min-max, get median
                        components = rt_attribute.split('-')
                        rt_min = float(components[0])
                        rt_max = float(components[1])
                        rt = (rt_min + rt_max) / 2.
                    else:
                        # RT is a single value
                        rt = float(rt_attribute)
                    break
            if rt is None:
                print('Cannot find RT value for MS2 spectrum: {0}'.format(id))
                exit()

            # load mz and energy level
            precursorList_node = xml.find('spectrumDesc').find('precursorList')
            precursor_nodes = precursorList_node.findall('precursor')
            if len(precursor_nodes) > 1:
                print('WARNING: MS2 spectrum with {0} precursors found: {1}'.format(len(precursor_nodes), id))
            mz = None
            energy = None
            for precursor_node in precursor_nodes:
                # find mz
                cvParam_nodes = precursor_node.find('ionSelection').findall('cvParam')
                for cvParam_node in cvParam_nodes:
                    if cvParam_node.attrib['name'] == 'MassToChargeRatio':
                        mz = float(cvParam_node.attrib['value'])
                        break

                # find energy
                cvParam_nodes = precursor_node.find('activation').findall('cvParam')
                for cvParam_node in cvParam_nodes:
                    if cvParam_node.attrib['name'] == 'CollisionEnergy':
                        energy = cvParam_node.attrib['value']

            # load masses
            mzArrayBinary_node = xml.find('mzArrayBinary').find('data')
            mz_data = mzArrayBinary_node.text
            if mz_data is None:
                print('WARNING: MS2 spectrum without mz data found: {0}'.format(id))
                return
            mz_decoded = standard_b64decode(mz_data)
            mz_count = len(mz_decoded) // 8
            masses = unpack('<{0}d'.format(mz_count), mz_decoded)

            # load intensities
            intenArrayBinary_nodes = xml.find('intenArrayBinary').find('data')
            intensities_data = intenArrayBinary_nodes.text
            intensities_decoded = standard_b64decode(intensities_data)
            intensities = unpack('<{0}f'.format(mz_count), intensities_decoded)

            # use only masses that are smaller than 1.00001 * mz?
            # masses = []
            # intensities = []
            # for i in range(0, len(masses)):
            #     if masses[i] < 1.00001 * mz:
            #         masses.append(masses[i])
            #         intensities.append(intensities[i])

            self.ms2_spectra.append(MS2Spectrum(id, mz, rt, energy, masses, intensities))

    def __str__(self):
        return '\n'.join([
            'MZ data from file: {0}'.format(self.file_name),
            'MS1 spectra: {0}'.format(len(self.ms1_spectra)),
            '\n'.join([str(spectrum) for spectrum in self.ms1_spectra]),
            '---',
            'MS2 spectra: {0}'.format(len(self.ms2_spectra)),
            '\n'.join([str(spectrum) for spectrum in self.ms2_spectra]),
            '---'
        ])
