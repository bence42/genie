from urllib.request import urlopen
import logging
import typing
import time
import json
import sys
import xml.etree.ElementTree as ET


debug_handler = logging.FileHandler("./debug.log")
debug_handler.setLevel(logging.DEBUG)

console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s [%(levelname)-5.5s] %(message)s',
    handlers=[
        debug_handler,
        console_handler
    ]
)


class ClinVarVariation:
    title = ""
    clinical_significance = ""
    accession = ""
    accession_version = ""
    cdna_change = ""


class Mutation:
    def __init__(self, line: str):
        tokens = line.split(sep='\t')
        self.chromosome = tokens[0].replace('chr', '')  # chr13 -> 13
        self.base_position = tokens[1]
        self.reference_base = tokens[2]
        self.variant_base = tokens[3]
        logging.debug(
            f'\nself.chromosome: {self.chromosome}\nself.base_position: {self.base_position}\nself.reference_base: {self.reference_base}\nself.variant_base: {self.variant_base}')
        self.validate(line)

        self.variation_ids: list[str] = []
        self.variations: list[ClinVarVariation] = []
        self.best_match: typing.Optional[ClinVarVariation] = None

    def validate(self, line: str):
        if not self.chromosome.isnumeric():
            logging.error(
                f'parsing failed: first column assumed to be chromosome number, got "{self.chromosome}"')
            logging.debug(f'line was: "{line}"')
            exit(1)

        if not self.base_position.isnumeric():
            logging.error(
                f'parsing failed: second column assumed to be position, got "{self.base_position}"')
            logging.debug(f'line was: "{line}"')
            exit(1)

        if not self.reference_base.isalpha() or len(self.reference_base) != 1:
            logging.error(
                f'parsing failed: third column assumed to be a 1 letter base, got "{self.reference_base}"')
            logging.debug(f'line was: "{line}"')
            exit(1)

        if not self.variant_base.isalpha() or len(self.variant_base) != 1:
            logging.error(
                f'parsing failed: third column assumed to be a 1 letter base, got "{self.variant_base}"')
            logging.debug(f'line was: "{line}"')
            exit(1)

    def search_clinvar(self):
        # ClinVar records do not store the mutation (deletion, SNV, insertion etc) in a serachable way.
        # Such a search is not possible: chr13:32890572 AND SNV AND G>A
        # The esearch API returns all records identified by chr13:32890572 and those
        # records need to be studied one by one to find the e.g G>A SNV mutation.
        self.variation_ids = ClinVarAPI.find_all_variations(m)
        if len(self.variation_ids) == 0:
            logging.warning(
                f'variation list for "{m.chromosome} {m.base_position}" is empty')
            return

        for variation_id in self.variation_ids:
            variation = ClinVarAPI.find_specific_variation(variation_id)
            self.variations.append(variation)

            if(f'{self.reference_base}>{self.variant_base}' in variation.cdna_change):
                self.best_match = variation
                var_with_accession = ClinVarAPI.get_variation_accession_version(
                    variation_id)
                if variation.accession == var_with_accession.accession:
                    variation.accession_version = var_with_accession.accession_version
                else:
                    logging.warning(
                        f'cannot find accession version for {self.chromosome} {self.base_position} (VCV{variation.accession})')
                logging.info(
                    f'found match for chr{self.chromosome}:{self.base_position} {self.reference_base}>{self.variant_base}: {self.best_match.title}, {self.best_match.accession}.{self.best_match.accession_version}, {self.best_match.clinical_significance}')

        if not self.best_match:
            print('\n', end='')
            logging.error(
                f'cannot find matching ClinVar recrod for "chr{m.chromosome}:{m.base_position} {m.reference_base}>{m.variant_base}"')
            logging.info(f'candidates are:')
            for variation in self.variations:
                print(
                    f'{variation.accession} - {variation.title} {variation.cdna_change}')
            print('\n', end='')


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
        esearch_response = urlopen(esearch_request)
        if not esearch_response.status == 200:
            logging.error('esearch request failed: {esearch_response.status}')
            exit(1)
        end_time = time.time()
        logging.debug(
            f'esearch request took: {(end_time - start_time) * 1000:.0f}ms')
        variation_ids = json.loads(esearch_response.read())[
            'esearchresult']['idlist']
        return variation_ids

    @staticmethod
    def find_specific_variation(variation_id: str | int):
        # e.g: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=822206&retmode=json
        esummary_request = f'{ClinVarAPI.BASE_ESUMMARY_URL}&id={variation_id}&retmode=json'
        logging.debug(esummary_request)
        start_time = time.time()
        esummary_response = urlopen(esummary_request)
        if not esummary_response.status == 200:
            logging.error('efetch request failed: {esummary_response.status}')
            exit(1)
        end_time = time.time()
        logging.debug(
            f'esummary request took: {(end_time - start_time) * 1000:.0f}ms')
        json_data = json.loads(esummary_response.read())
        var = ClinVarVariation()
        var.title = json_data['result'][variation_id]['title']
        var.clinical_significance = json_data['result'][variation_id]['clinical_significance']['description']
        var.accession = json_data['result'][variation_id]['accession']
        var.cdna_change = json_data['result'][variation_id]['variation_set'][0]['cdna_change']
        return var

    @staticmethod
    def get_variation_accession_version(variation_id: str | int):
        # e.g: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id=125889
        efetch_request = f'{ClinVarAPI.BASE_EFETCH_URL}&rettype=vcv&is_variationid&id={variation_id}'
        logging.debug(efetch_request)
        start_time = time.time()
        efetch_response = urlopen(efetch_request)
        if not efetch_response.status == 200:
            logging.error('efetch request failed: {efetch_response.status}')
            exit(1)
        end_time = time.time()
        logging.debug(
            f'efetch request took: {(end_time - start_time) * 1000:.0f}ms')
        root = ET.fromstring(efetch_response.read())
        var = ClinVarVariation()
        var.accession = root[0].attrib['Accession']
        var.accession_version = root[0].attrib['Version']
        return var


with open('./BRCA/BRCA_2022/1_22.txt', 'r') as input_file:
    logging.debug(f'file: {input_file.name}')
    for line in input_file.readlines()[1:]:
        if len(line) == 0:
            continue
        m = Mutation(line)
        m.search_clinvar()
