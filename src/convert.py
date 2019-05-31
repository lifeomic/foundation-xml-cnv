#!/usr/bin/env python
import argparse
import csv
import json, ast
import logging

import xmltodict

logging.basicConfig(level=logging.INFO,
                    format='[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s')
logger = logging.getLogger(__name__)


def read_xml(xml_file):
    with open(xml_file) as fd:
        return xmltodict.parse(fd.read())


def calculate_status(equivocal, copy_type):
    if copy_type == 'amplification':
        if equivocal == 'true':
            return 'gain'
        return 'amplification'
    if copy_type == 'loss':
        if equivocal == 'true':
            return 'partial_loss'
        return 'loss'
    if copy_type == 'partial amplification':
        return 'gain'

    logger.error('Failed to resolve copy type: %s, equivocal: %s', copy_type, equivocal)
    return ''


def calculate_interpretation(status):
    if status == 'known':
        return 'Pathogenic'
    if status == 'likely':
        return 'Likely pathogenic'
    if status == 'unknown':
        return 'Uncertain significance'
    if status == 'ambiguous':
        return 'other'

    logger.error('Failed to resolve interpretation: %s', status)
    return ''


def gather_attributes(copy_number):
    attributes = {}
    if '@number-of-exons' in copy_number.keys():
        attributes['number-of-exons'] = copy_number['@number-of-exons']
    if '@type' in copy_number.keys():
        attributes['status'] = copy_number['@type']
    if '@ratio' in copy_number.keys():
        attributes['ratio'] = copy_number['@ratio']
    if '@status' in copy_number.keys():
        attributes['interpretation'] = copy_number['@status']

    return attributes


def extract_sample(samples):
    if not samples:
        return None

    sample = samples.get('samples', {}).get('sample', {})
    if isinstance(sample, list):
        # Multiple sample (rna/dna) provided, find dna:
        dna_samples = [s for s in sample if s.get('@nucleic-acid-type', None) == 'DNA']
        if dna_samples and len(dna_samples) >= 1:
            return dna_samples[0].get('@name', None)
        return None

    return sample.get('@name', None)


def extract_copy_numbers(results_payload_dict):
    logger.info('Extracting copy numbers from xml')
    copy_number_list = {'CopyNumbers': []}

    if 'copy-number-alterations' in results_payload_dict['variant-report'].keys():
        if (results_payload_dict['variant-report']['copy-number-alterations'] is not None and
                'copy-number-alteration' in results_payload_dict['variant-report']['copy-number-alterations'].keys()):

            sample_id = extract_sample(results_payload_dict['variant-report'])
            variants_dict = results_payload_dict['variant-report']['copy-number-alterations']['copy-number-alteration']
            copy_numbers = variants_dict if isinstance(variants_dict, list) else [variants_dict]

            for copy_number in copy_numbers:
                copy_number_value = {'sample_id': copy_number.get('dna-evidence', {}).get('@sample', sample_id),
                                     'gene': copy_number['@gene'],
                                     'copy_number': float(format(copy_number['@copy-number'])),
                                     'status': calculate_status(copy_number['@equivocal'], copy_number['@type']),
                                     'chromosome': copy_number['@position'].split(":")[0],
                                     'start_position': copy_number['@position'].split(":")[1].split("-")[0],
                                     'end_position': copy_number['@position'].split(":")[1].split("-")[1],
                                     'attributes': gather_attributes(copy_number),
                                     'interpretation': calculate_interpretation(copy_number['@status'])}

                copy_number_list['CopyNumbers'].append(ast.literal_eval(json.dumps(copy_number_value)))

    return copy_number_list


def write_copy_numbers_to_cnv(cnv_dict, args):
    logger.info('Saving copy numbers to cnv file')

    with open(args.out_file, 'wb') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=',')
        csv_writer.writerow(['sample_id', 'gene', 'copy_number', 'status', 'attributes',
                             'chromosome', 'start_position', 'end_position', 'interpretation'])
        for cnv in cnv_dict['CopyNumbers']:
            csv_writer.writerow([cnv['sample_id'], cnv['gene'], cnv['copy_number'],
                                 cnv['status'], cnv['attributes'], cnv['chromosome'],
                                 cnv['start_position'], cnv['end_position'], cnv['interpretation']])


def main():
    parser = argparse.ArgumentParser(
        prog='foundation-xml-cnv',
        description='Extracts copy number information from FoundationOne XML reports into CSV resources.')
    parser.add_argument('-x, --xml', dest='xml_file',
                        required=True, help='Path to the XML file')
    parser.add_argument('-o, --output', dest='out_file',
                        required=True, help='Path to write the CNV file')
    args = parser.parse_args()

    logger.info('Extracting copy numbers from XML and into CNV with args: %s', json.dumps(args.__dict__))
    xml_dict = read_xml(args.xml_file)

    cnv_dict = extract_copy_numbers(
        xml_dict['rr:ResultsReport']['rr:ResultsPayload'])

    write_copy_numbers_to_cnv(cnv_dict, args)


if __name__ == '__main__':
    main()
