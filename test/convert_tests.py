from unittest import TestCase
from src.convert import extract_copy_numbers
from src.convert import calculate_status
from src.convert import gather_attributes
from src.convert import calculate_interpretation
from src.convert import extract_sample

expected_results = {
    'CopyNumbers': [{
        'status': 'amplification',
        'sample_id': 'SA-1612348',
        'start_position': '58093932',
        'end_position': '58188144',
        'copy_number': 44.0,
        'gene': 'CDK4',
        'chromosome': 'chr12',
        'interpretation': 'Pathogenic',
        'attributes': {'number-of-exons': '7 of 7', 'ratio': 11.63, 'interpretation': 'known',
                       'status': 'amplification'}
    }, {
        'status': 'gain',
        'sample_id': 'SA-1612348',
        'start_position': '41853880',
        'end_position': '41956362',
        'copy_number': 6.0,
        'gene': 'CCND3',
        'chromosome': 'chr6',
        'interpretation': 'Likely pathogenic',
        'attributes': {'number-of-exons': '5 of 5', 'ratio': 2.17, 'interpretation': 'likely',
                       'status': 'amplification'}
    }, {
        'status': 'loss',
        'sample_id': 'SA-1612348',
        'start_position': '128706589',
        'end_position': '128801451',
        'copy_number': 41.0,
        'gene': 'MYC',
        'chromosome': 'chr8',
        'interpretation': 'Uncertain significance',
        'attributes': {'number-of-exons': '5 of 5', 'ratio': 10.34, 'interpretation': 'unknown',
                       'status': 'loss'}
    }, {
        'status': 'partial_loss',
        'sample_id': 'SA-1612348',
        'start_position': '37138078',
        'end_position': '37141867',
        'copy_number': 6.0,
        'gene': 'PIM1',
        'chromosome': 'chr6',
        'interpretation': 'other',
        'attributes': {'number-of-exons': '7 of 7', 'ratio': 2.14, 'interpretation': 'ambiguous',
                       'status': 'loss'}
    }, {
        'status': 'gain',
        'sample_id': 'SA-1612348',
        'start_position': '117859738',
        'end_position': '117878968',
        'copy_number': 7.0,
        'gene': 'RAD21',
        'chromosome': 'chr8',
        'interpretation': 'Pathogenic',
        'attributes': {'number-of-exons': '13 of 13', 'ratio': 2.69,
                       'interpretation': 'known', 'status': 'partial amplification'}
    }]
}

expected_results_multiple_samples = {
    'CopyNumbers': [{
        'status': 'amplification',
        'sample_id': 'SA-1612348',
        'start_position': '58093932',
        'end_position': '58188144',
        'copy_number': 44.0,
        'gene': 'CDK4',
        'chromosome': 'chr12',
        'interpretation': 'Pathogenic',
        'attributes': {'number-of-exons': '7 of 7', 'ratio': 11.63, 'interpretation': 'known',
                       'status': 'amplification'}
    }, {
        'status': 'gain',
        'sample_id': 'SA-1612348',
        'start_position': '41853880',
        'end_position': '41956362',
        'copy_number': 6.0,
        'gene': 'CCND3',
        'chromosome': 'chr6',
        'interpretation': 'Likely pathogenic',
        'attributes': {'number-of-exons': '5 of 5', 'ratio': 2.17, 'interpretation': 'likely',
                       'status': 'amplification'}
    }, {
        'status': 'loss',
        'sample_id': 'SA-1612348',
        'start_position': '128706589',
        'end_position': '128801451',
        'copy_number': 41.0,
        'gene': 'MYC',
        'chromosome': 'chr8',
        'interpretation': 'Uncertain significance',
        'attributes': {'number-of-exons': '5 of 5', 'ratio': 10.34, 'interpretation': 'unknown',
                       'status': 'loss'}
    }, {
        'status': 'partial_loss',
        'sample_id': 'SA-1612348',
        'start_position': '37138078',
        'end_position': '37141867',
        'copy_number': 6.0,
        'gene': 'PIM1',
        'chromosome': 'chr6',
        'interpretation': 'other',
        'attributes': {'number-of-exons': '7 of 7', 'ratio': 2.14, 'interpretation': 'ambiguous',
                       'status': 'loss'}
    }, {
        'status': 'gain',
        'sample_id': 'sample1',
        'start_position': '117859738',
        'end_position': '117878968',
        'copy_number': 7.0,
        'gene': 'RAD21',
        'chromosome': 'chr8',
        'interpretation': 'Pathogenic',
        'attributes': {'number-of-exons': '13 of 13', 'ratio': 2.69,
                       'interpretation': 'known', 'status': 'partial amplification'}
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
                '@name': 'SA-1612348'
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
                    '@status': 'likely',
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
                    '@status': 'unknown',
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
                    '@status': 'ambiguous',
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
                    '@number-of-exons': '13 of 13'
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

foundation_source_dict_multiple_samples = {
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
            'sample': [
                {
                    '@bait-set': 'D2',
                    '@mean-exon-depth': '811.98',
                    '@name': 'sample1',
                    '@nucleic-acid-type': 'DNA'
                },
                {
                    '@bait-set': 'R2',
                    '@mean-exon-depth': '413.25',
                    '@name': 'sample2',
                    '@nucleic-acid-type': 'RNA'
                }
            ]
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
                    '@status': 'likely',
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
                    '@status': 'unknown',
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
                    '@status': 'ambiguous',
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
                    '@number-of-exons': '13 of 13'
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

    def test_convert_with_multiple_samples(self):
        cnv_resources = extract_copy_numbers(foundation_source_dict_multiple_samples)
        self.maxDiff = None
        self.assertDictEqual(expected_results_multiple_samples, cnv_resources)

    def test_calculate_status(self):
        self.assertEqual('gain', calculate_status('true', 'amplification'))
        self.assertEqual('amplification', calculate_status('false', 'amplification'))
        self.assertEqual('partial_loss', calculate_status('true', 'loss'))
        self.assertEqual('loss', calculate_status('false', 'loss'))
        self.assertEqual('gain', calculate_status('false', 'partial amplification'))
        self.assertEqual('gain', calculate_status('true', 'partial amplification'))
        self.assertEqual('', calculate_status('true', 'fred'))

    def test_calculate_interpretation(self):
        self.assertEqual('Pathogenic', calculate_interpretation('known'))
        self.assertEqual('Likely pathogenic', calculate_interpretation('likely'))
        self.assertEqual('Uncertain significance', calculate_interpretation('unknown'))
        self.assertEqual('other', calculate_interpretation('ambiguous'))
        self.assertEqual('', calculate_interpretation('fred'))

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
        self.assertEqual(
            {'number-of-exons': '13 of 13', 'status': 'partial amplification', 'ratio': 2.69,
             'interpretation': 'known'}, gather_attributes(copy_number))

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
        self.assertEqual({'ratio': 2.69, 'interpretation': 'known', 'status': 'loss'}, gather_attributes(copy_number))

    def test_gather_attributes_with_no_ratio(self):
        copy_number = {
            '@gene': 'RAD21',
            '@position': 'chr8:117859738-117878968',
            '@copy-number': '7',
            '@equivocal': 'true',
            '@status': 'known',
            '@type': 'loss',
            'dna-evidence': {
                '@sample': 'SA-1612348'
            }
        }
        self.assertEqual({'interpretation': 'known', 'status': 'loss'}, gather_attributes(copy_number))

    def test_gather_attributes_with_no_status(self):
        copy_number = {
            '@gene': 'RAD21',
            '@position': 'chr8:117859738-117878968',
            '@copy-number': '7',
            '@equivocal': 'true',
            '@type': 'loss',
            'dna-evidence': {
                '@sample': 'SA-1612348'
            }
        }
        self.assertEqual({'status': 'loss'}, gather_attributes(copy_number))

    def test_gather_attributes_with_no_type(self):
        copy_number = {
            '@gene': 'RAD21',
            '@position': 'chr8:117859738-117878968',
            '@copy-number': '7',
            '@equivocal': 'true',
            'dna-evidence': {
                '@sample': 'SA-1612348'
            }
        }
        self.assertEqual({}, gather_attributes(copy_number))

    def test_extract_samples_empty_sample(self):
        samples = {
            'samples': {
                'sample': [
                ]
            }
        }
        self.assertIsNone(extract_sample(samples))

    def test_extract_samples_missing_dna_sample(self):
        samples = {
            'samples': {
                'sample': [
                    {
                        '@bait-set': 'D2',
                        '@mean-exon-depth': '811.98',
                        '@name': 'sample1',
                        '@nucleic-acid-type': 'NOT-DNA'
                    },
                    {
                        '@bait-set': 'R2',
                        '@mean-exon-depth': '413.25',
                        '@name': 'sample2',
                    }
                ]
            }
        }
        self.assertIsNone(extract_sample(samples))

    def test_extract_samples_multiple(self):
        samples = {
            'samples': {
                'sample': [
                    {
                        '@bait-set': 'D2',
                        '@mean-exon-depth': '811.98',
                        '@name': 'sample1',
                        '@nucleic-acid-type': 'DNA'
                    },
                    {
                        '@bait-set': 'R2',
                        '@mean-exon-depth': '413.25',
                        '@name': 'sample2',
                        '@nucleic-acid-type': 'RNA'
                    }
                ]
            }
        }
        self.assertEquals('sample1', extract_sample(samples))

    def test_extract_samples(self):
        samples = {
            'samples': {
                'sample': {
                    '@name': 'SA-1612348'
                }
            }
        }
        self.assertEquals('SA-1612348', extract_sample(samples))
