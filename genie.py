from urllib.request import urlopen
import logging
import time
import json
import sys
import argparse
from pathlib import Path
import xml.etree.ElementTree as ET


debug_handler = logging.FileHandler("./debug.log")
debug_handler.setLevel(logging.DEBUG)

console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s %(levelname)-0.11s %(message)s',
    handlers=[
        debug_handler,
        console_handler
    ]
)
logging.addLevelName(logging.CRITICAL, ' [ FATAL ]')
logging.addLevelName(logging.ERROR,    ' [ ERROR ]')
logging.addLevelName(logging.WARNING,  ' [ WARN  ]')
logging.addLevelName(logging.INFO, '')
logging.addLevelName(logging.DEBUG,    '   [ DEBUG ]')


class ClinVarVariation:
    id: int = -1  # 125749 etc
    title = ""  # NM_007294.4(BRCA1):c.4987-92A>G / NM_007294.4(BRCA1):c.4987-93_4987-92inv / NM_007294.4(BRCA1):c.4987-94_4987-92delinsGCG / etc # nopep8
    clinical_significance = ""  # Benign / Pathogenic
    accession = ""  # e.g VCV000125749 / VCV000491093
    obj_type = ""  # single nucleotide variant / Inversion / Indel etc
    canonical_spdi = ""  # NC_000017.11:43067786:T:C"
    accession_version = ""  # 1/2/../23 etc
    cdna_change = ""  # NM_007294.4(BRCA1):c.4987-92A>G / NM_007294.4(BRCA1):c.4987-93_4987-92inv / NM_007294.4(BRCA1):c.4987-94_4987-92delinsGCG / etc # nopep8
    grch38_position: int = 0


class Mutation:
    def __init__(self, line: str):
        tokens = line.split(sep='\t')
        self.chromosome = tokens[0].replace('chr', '')  # chr13 -> 13
        self.base_position = tokens[1]
        self.reference_base = tokens[2]
        self.variant_base = tokens[3]
        self.frequency = tokens[6]
        self.type = tokens[9]
        logging.debug(
            f'\nself.chromosome: {self.chromosome}\nself.base_position: {self.base_position}\nself.reference_base: {self.reference_base}\nself.variant_base: {self.variant_base}')
        self.validate(line)
        self.translate_type()

        self.variation_ids: list[str] = []
        self.variations: list[ClinVarVariation] = []

    def translate_type(self):
        iontorrent_to_clinvar = {
            'SNP': 'single nucleotide variant'}
        if self.type in iontorrent_to_clinvar:
            self.type = iontorrent_to_clinvar[self.type]

    def validate(self, line: str):
        if not self.chromosome.isnumeric():
            logging.error(
                f'parsing failed: first column (Chrom) expected to be chromosome number, got "{self.chromosome}"')
            logging.debug(f'line was: "{line}"')
            exit(1)

        if not self.base_position.isnumeric():
            logging.error(
                f'parsing failed: second column (Position) expected to be position, got "{self.base_position}"')
            logging.debug(f'line was: "{line}"')
            exit(1)

        if not self.reference_base.isalpha() and self.reference_base != '-':
            logging.error(
                f'parsing failed: third column (Ref) expected to be a base or "-", got "{self.reference_base}"')
            logging.debug(f'line was: "{line}"')
            exit(1)

        if not self.variant_base.isalpha() and self.variant_base != '-':
            logging.error(
                f'parsing failed: fourth column (Variant) expected to be a base or "-", got "{self.variant_base}"')
            logging.debug(f'line was: "{line}"')
            exit(1)

        if not self.frequency.replace('.', '').isnumeric():
            logging.error(
                f'parsing failed: seventh column (Frequency) expected to be a floating point number", got "{self.frequency}"')
            logging.debug(f'line was: "{line}"')
            exit(1)

        if int(self.frequency.split('.')[0]) > 100:
            logging.error(
                f'parsing failed: seventh column (Frequency) expected to be in [0-100]", got "{self.frequency}"')
            logging.debug(f'line was: "{line}"')
            exit(1)

        if self.type not in ['SNP', 'INS', 'DEL', 'COMPLEX', 'MNP']:
            logging.error(
                f'parsing failed: tenth column (Type) expected to be SNP/DEL/INS/COMPLEX, got "{self.type}"')
            logging.debug(f'line was: "{line}"')
            exit(1)

    def search_clinvar(self) -> ClinVarVariation | None:
        # ClinVar records do not store the mutation (deletion, SNV, insertion etc) in a serachable way.
        # E.g a search like this is not possible: "chr13:32890572 AND SNV AND G>A"
        # The esearch API returns all records identified by chr13:32890572 and those
        # records need to be studied one by one to find the e.g G>A SNV mutation.
        self.variation_ids = ClinVarAPI.find_all_variations(self)
        if len(self.variation_ids) == 0:
            logging.error(
                f' > cannot find "chr{self.chromosome}:{self.base_position}" in ClinVar')
            return None

        if match := self.select_best_match():
            return match

        # If the mutation is COMPLEX or we can't match, bail out and just give the probable candidates for human review
        logging.error(
            f' > cannot find matching ClinVar recrod for "chr{m.chromosome}:{m.base_position} {m.reference_base}>{m.variant_base}"')
        logging.info(f'           > candidates are:')
        for variation in self.variations:
            logging.info(
                f'             > {variation.accession} - {variation.title}')
        return None

    def select_best_match(self) -> ClinVarVariation | None:
        for variation_id in self.variation_ids:
            variation = ClinVarAPI.find_specific_variation(variation_id)
            self.variations.append(variation)

            # cdna_change ofthen contains the mutation, this is the most straightforward way to match
            if(f'{self.reference_base}>{self.variant_base}' in variation.cdna_change):
                self.find_accession_version(variation)
                return variation

            # Try again as sometimes an A>G SNP is stored as T>C, but the canonical_spdi will have A>G
            if variation.obj_type == self.type and f'{self.base_position}:{self.reference_base}:{self.variant_base}' in variation.canonical_spdi:
                self.find_accession_version(variation)
                logging.warning(
                    f' > chr{self.chromosome}:{self.base_position} matched using canonical_spdi, review needed')
                return variation

            # Try again as sometimes an A>G SNP is stored as T>C, but the canonical_spdi will have A>G with GRCh38_base_position-1
            if variation.obj_type == self.type and f'{int(variation.grch38_position)-1}:{self.reference_base}:{self.variant_base}' in variation.canonical_spdi:
                self.find_accession_version(variation)
                logging.warning(
                    ' > matched using canonical_spdi and the corresponding GRCh38 base position, review needed')
                return variation
        return None

    def find_accession_version(self, variation: ClinVarVariation):
        # This takes another ~600ms so we only do this when there's a match
        var_with_accession = ClinVarAPI.get_variation_accession_version(
            variation.id)
        if variation.accession == var_with_accession.accession:
            variation.accession_version = var_with_accession.accession_version
        else:
            logging.warning(
                f' cannot find accession version for {self.chromosome} {self.base_position} (VCV{variation.accession})')
        logging.debug(
            f'found match for chr{self.chromosome}:{self.base_position} {self.reference_base}>{self.variant_base}: {variation.title}, {variation.accession}.{variation.accession_version}, {self.frequency}%, {variation.clinical_significance}')


class ClinVarAPI:
    BASE_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
    BASE_ESEARCH_URL = f"{BASE_URL}/esearch.fcgi?db=clinvar"
    BASE_ESUMMARY_URL = f"{BASE_URL}/esummary.fcgi?db=clinvar"
    BASE_EFETCH_URL = f"{BASE_URL}/efetch.fcgi?db=clinvar"

    @staticmethod
    def find_all_variations(mutation: Mutation) -> list[str]:
        # e.g: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=13[Chromosome]+AND+32906729[Base+Position+for+Assembly+GRCh37]&retmode=json
        esearch_request = f'{ClinVarAPI.BASE_ESEARCH_URL}&term={mutation.chromosome}[Chromosome]+AND+{mutation.base_position}[Base+Position+for+Assembly+GRCh37]&retmode=json'
        logging.debug(esearch_request)
        start_time = time.time()
        with urlopen(esearch_request) as esearch_response:
            if not esearch_response.status == 200:
                logging.error(
                    'esearch request failed: {esearch_response.status}')
                exit(1)
            end_time = time.time()
            logging.debug(
                f'esearch request took: {(end_time - start_time) * 1000:.0f}ms')
            variation_ids = json.loads(esearch_response.read())[
                'esearchresult']['idlist']
        return variation_ids

    @staticmethod
    def find_specific_variation(variation_id: str | int) -> ClinVarVariation:
        # e.g: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=822206&retmode=json
        esummary_request = f'{ClinVarAPI.BASE_ESUMMARY_URL}&id={variation_id}&retmode=json'
        logging.debug(esummary_request)
        start_time = time.time()
        with urlopen(esummary_request) as esummary_response:
            if not esummary_response.status == 200:
                logging.error(
                    'efetch request failed: {esummary_response.status}')
                exit(1)
            end_time = time.time()
            logging.debug(
                f'esummary request took: {(end_time - start_time) * 1000:.0f}ms')
            json_data = json.loads(esummary_response.read())[
                'result'][variation_id]

        var = ClinVarVariation()
        var.id = int(variation_id)
        var.title = json_data['title']
        var.clinical_significance = json_data['clinical_significance']['description']
        var.accession = json_data['accession']
        var.obj_type = json_data['obj_type']
        var.cdna_change = json_data['variation_set'][0]['cdna_change']

        for variation_loc in json_data['variation_set'][0]['variation_loc']:
            if 'GRCh38' in variation_loc['assembly_name']:
                var.grch38_position = variation_loc['start']
                break

        var.canonical_spdi = json_data['variation_set'][0]['canonical_spdi']
        return var

    @staticmethod
    def get_variation_accession_version(variation_id: str | int) -> ClinVarVariation:
        # e.g: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id=125889
        efetch_request = f'{ClinVarAPI.BASE_EFETCH_URL}&rettype=vcv&is_variationid&id={variation_id}'
        logging.debug(efetch_request)
        start_time = time.time()
        with urlopen(efetch_request) as efetch_response:
            if not efetch_response.status == 200:
                logging.error(
                    'efetch request failed: {efetch_response.status}')
                exit(1)
            end_time = time.time()
            logging.debug(
                f'efetch request took: {(end_time - start_time) * 1000:.0f}ms')
            root = ET.fromstring(efetch_response.read())

        var = ClinVarVariation()
        var.accession = root[0].attrib['Accession']
        var.accession_version = root[0].attrib['Version']
        return var


def write_csv(results_file_name: str, search_and_results: list[tuple[Mutation, ClinVarVariation | None]]):
    with open(f'./results/{results_file_name}', 'w') as output:
        output.write(
            f'Chromosome;Base position;ClinVar title;ClinVar accession;Frequency;Clinical significance')
        for mutation, clinvar_record in search_and_results:
            if clinvar_record:
                output.write(
                    f'\nchr{mutation.chromosome};{mutation.base_position};{clinvar_record.title};{clinvar_record.accession}.{clinvar_record.accession_version};{mutation.frequency}%;{clinvar_record.clinical_significance}')
            else:
                output.write(
                    f'\nchr{mutation.chromosome};{mutation.base_position};;;;')


def get_linecount(input_file):
    line_count = sum(1 for line in input_file if line.strip())
    input_file.seek(0)
    return line_count - 1  # first line is a header


if __name__ == '__main__':

    arg_parser = argparse.ArgumentParser(
        description='Queries the ClinVar database for mutations identified by Ion Torrent')
    arg_parser.add_argument('--input', '-i',
                            metavar='INPUT',
                            dest='input_files',
                            nargs="+",
                            type=str,
                            help='the path to the Ion Torrent raw output file',
                            required=True)

    args = arg_parser.parse_args()

    for file_idx, input_file_path in enumerate(args.input_files):
        with open(input_file_path, 'r') as input_file:
            logging.info(f'[{file_idx+1:>3}/{len(args.input_files):<3}]{input_file.name}')
    line_count = get_linecount(input_file)

            mutations_and_records: list[tuple[Mutation,
                                              ClinVarVariation | None]] = []

            for line_idx, line in enumerate(input_file.readlines()[1:]):
        if len(line) == 0:
            continue
        m = Mutation(line)
        logging.info(
                    f'[{line_idx+1:>3}/{line_count:<3}] chr{m.chromosome}:{m.base_position}')
        clinvar_record = m.search_clinvar()
                mutations_and_records.append((m, clinvar_record))

    timestamp = time.strftime('%Y%m%d-%H%M%S')
            results_file_name = Path(
                input_file.name).stem + f'_{timestamp}_results.csv'
            write_csv(results_file_name, mutations_and_records)
