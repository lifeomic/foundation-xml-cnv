from unittest import TestCase
from src.convert import extract_copy_numbers
from src.convert import calculate_status
from src.convert import gather_attributes

expected_results = {
    'CopyNumbers': [{
        'status': 'amplification',
        'sample_id': 'SA-1612348',
        'start_position': '58093932',
        'end_position': '58188144',
        'copy_number': 44.0,
        'gene': 'CDK4',
        'chromosome': 'chr12',
        'attributes': {'number-of-exons': '7 of 7'}
    }, {
        'status': 'gain',
        'sample_id': 'SA-1612348',
        'start_position': '41853880',
        'end_position': '41956362',
        'copy_number': 6.0,
        'gene': 'CCND3',
        'chromosome': 'chr6',
        'attributes': {'number-of-exons': '5 of 5'}
    }, {
        'status': 'loss',
        'sample_id': 'SA-1612348',
        'start_position': '128706589',
        'end_position': '128801451',
        'copy_number': 41.0,
        'gene': 'MYC',
        'chromosome': 'chr8',
        'attributes': {'number-of-exons': '5 of 5'}
    }, {
        'status': 'partial_loss',
        'sample_id': 'SA-1612348',
        'start_position': '37138078',
        'end_position': '37141867',
        'copy_number': 6.0,
        'gene': 'PIM1',
        'chromosome': 'chr6',
        'attributes': {'number-of-exons': '7 of 7'}
    }, {
        'status': 'gain',
        'sample_id': 'SA-1612348',
        'start_position': '117859738',
        'end_position': '117878968',
        'copy_number': 7.0,
        'gene': 'RAD21',
        'chromosome': 'chr8',
        'attributes': {'partial amplification': True, 'number-of-exons': '13 of 13'}
    }]
}

foundation_source_dict = {
    'FinalReport': {
        'PMI': {
            'LastName': 'doe',
            'FirstName': 'jane',
            'MRN': '1234',
            'Gender': 'Female',
            'DOB': '1970-01-01',
            'CollDate': '2000-01-01',
            'SubmittedDiagnosis': 'Cancer'
        },
        'Sample': {
            'TestType': 'Test 1'
        }
    },
    'variant-report': {
        'samples': {
            'sample': {
                '@name': 'sample1'
            }
        },
        'short-variants': {
            'short-variant': [
                {
                    '@gene': 'gene1',
                    '@cds-effect': '229C&gt;A',
                    '@functional-effect': 'missense',
                    '@allele-fraction': 0.488,
                    '@position': 'chr1:100',
                    '@depth': 200,
                    '@transcript': 'NM_001',
                    '@status': 'known',
                    '@protein-effect': 'R77S'
                }
            ]
        },
        'copy-number-alterations': {
            'copy-number-alteration': [
                {
                    '@gene': 'CDK4',
                    '@position': 'chr12:58093932-58188144',
                    '@copy-number': '44',
                    '@equivocal': 'false',
                    '@ratio': 11.63,
                    '@status': 'known',
                    '@type': 'amplification',
                    '@number-of-exons': '7 of 7',
                    'dna-evidence': {
                        '@sample': 'SA-1612348'
                    }
                },
                {
                    '@gene': 'CCND3',
                    '@position': 'chr6:41853880-41956362',
                    '@copy-number': '6',
                    '@equivocal': 'true',
                    '@ratio': 2.17,
                    '@status': 'known',
                    '@type': 'amplification',
                    '@number-of-exons': '5 of 5',
                    'dna-evidence': {
                        '@sample': 'SA-1612348'
                    }
                },
                {
                    '@gene': 'MYC',
                    '@position': 'chr8:128706589-128801451',
                    '@copy-number': '41',
                    '@equivocal': 'false',
                    '@ratio': 10.34,
                    '@status': 'known',
                    '@type': 'loss',
                    '@number-of-exons': '5 of 5',
                    'dna-evidence': {
                        '@sample': 'SA-1612348'
                    }
                },
                {
                    '@gene': 'PIM1',
                    '@position': 'chr6:37138078-37141867',
                    '@copy-number': '6',
                    '@equivocal': 'true',
                    '@ratio': 2.14,
                    '@status': 'known',
                    '@type': 'loss',
                    '@number-of-exons': '7 of 7',
                    'dna-evidence': {
                        '@sample': 'SA-1612348'
                    }
                },
                {
                    '@gene': 'RAD21',
                    '@position': 'chr8:117859738-117878968',
                    '@copy-number': '7',
                    '@equivocal': 'true',
                    '@ratio': 2.69,
                    '@status': 'known',
                    '@type': 'partial amplification',
                    '@number-of-exons': '13 of 13',
                    'dna-evidence': {
                        '@sample': 'SA-1612348'
                    }
                }
            ]
        },
        'rearrangements': {
            'rearrangement': [
                {
                    '@status': 'known',
                    '@targeted-gene': 'CDK4',
                    '@type': 'truncation',
                    '@pos1': 'ch17:29557687-29887856',
                    '@pos2': 'ch6:66426718-66427149'
                }
            ]
        },
        'biomarkers': {
            'microsatellite-instability': {
                '@status': 'MSS'
            },
            'tumor-mutation-burden': {
                '@status': 'low'
            }
        }
    }
}


class ConvertTest(TestCase):

    def test_convert(self):
        cnv_resources = extract_copy_numbers(foundation_source_dict)
        self.maxDiff = None
        self.assertDictEqual(expected_results, cnv_resources)

    def test_calculate_status(self):
        self.assertEqual('gain', calculate_status('true', 'amplification'))
        self.assertEqual('amplification', calculate_status('false', 'amplification'))
        self.assertEqual('partial_loss', calculate_status('true', 'loss'))
        self.assertEqual('loss', calculate_status('false', 'loss'))
        self.assertEqual('gain', calculate_status('false', 'partial amplification'))
        self.assertEqual('gain', calculate_status('true', 'partial amplification'))
        self.assertEqual('', calculate_status('true', 'fred'))

    def test_gather_attributes_with_partial_amplification(self):
        copy_number = {
            '@gene': 'RAD21',
            '@position': 'chr8:117859738-117878968',
            '@copy-number': '7',
            '@equivocal': 'true',
            '@ratio': 2.69,
            '@status': 'known',
            '@type': 'partial amplification',
            '@number-of-exons': '13 of 13',
            'dna-evidence': {
                '@sample': 'SA-1612348'
            }
        }
        self.assertEqual({'number-of-exons': '13 of 13', 'partial amplification': True}, gather_attributes(copy_number))

    def test_gather_attributes_with_partial_amplification(self):
        copy_number = {
            '@gene': 'RAD21',
            '@position': 'chr8:117859738-117878968',
            '@copy-number': '7',
            '@equivocal': 'true',
            '@ratio': 2.69,
            '@status': 'known',
            '@type': 'loss',
            '@number-of-exons': '13 of 13',
            'dna-evidence': {
                '@sample': 'SA-1612348'
            }
        }
        self.assertEqual({'number-of-exons': '13 of 13'}, gather_attributes(copy_number))

    def test_gather_attributes_with_no_exons(self):
        copy_number = {
            '@gene': 'RAD21',
            '@position': 'chr8:117859738-117878968',
            '@copy-number': '7',
            '@equivocal': 'true',
            '@ratio': 2.69,
            '@status': 'known',
            '@type': 'loss',
            'dna-evidence': {
                '@sample': 'SA-1612348'
            }
        }
        self.assertEqual({}, gather_attributes(copy_number))

