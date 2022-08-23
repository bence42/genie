

import enum
import os
import re
import json
from pathlib import Path


def get_raw_file_list(dir_to_parse):
    raw_results = re.compile(r'\d{1,3}_?\d{1,3}?\.txt')
    raw_files: list[str] = []
    for root, dirs, files in os.walk(dir_to_parse):
        for file in files:
            if raw_results.match(file):
                raw_files.append(os.path.join(root, file))
    return raw_files


def find_matching_result_files(dir_to_parse, raw_files):
    raw_to_result: list[tuple[str, str]] = []
    processed = re.compile(r'BRCA_.*\.txt')
    for root, dirs, files in os.walk(dir_to_parse):
        for file in files:
            for raw_file in raw_files:
                processed_with_raw = re.compile(
                    f'BRCA_?{Path(raw_file).stem}.*\\.txt')
                if processed_with_raw.match(file):
                    raw_to_result.append((raw_file, os.path.join(root, file)))
    return raw_to_result


def file_len_match(path1, path2):
    with open(f'{path1}') as file1:
        file1_line_count = len(file1.readlines()) - 1  # first line is header
    with open(f'{path2}') as file2:
        file2_line_count = len(file2.readlines())
    return file1_line_count == file2_line_count


def create_record(result_line, raw_line) -> dict | None:
    raw_tokens = raw_line.split('\t')
    # "NM_007294.4(BRCA1):c.4308T>C (p.Ser1436=),VCV000125703.29	45.8%"
    # "00125703.29	45.8%"
    result_tokens = result_line.split('CV0')
    try:
        # "00125703.29	45.8%"
        # "VCV000125703.29"
        VCV_id = "VCV0" + result_tokens[1].split('\t')[0]
        # "VCV000125703.29"
        # "VCV000125703"
        VCV_id = VCV_id.split('.')[0]
    except:
        # if a result line cannot be parsed, it has no vcv_id,  so it's unusable
        return None
    record = {}
    record['chromosome'] = raw_tokens[0]
    record['position'] = raw_tokens[1]
    record['ref'] = raw_tokens[2]
    record['var'] = raw_tokens[3]
    record['vcv'] = VCV_id

    return record


def dump_records(raw_lines, result_lines, matched_file):
    records = []
    for idx, raw_line in enumerate(raw_lines):
        result_line = result_lines[idx]
        if record := create_record(result_line, raw_line):
            records.append(record)

    json.dump(records, matched_file)


if __name__ == '__main__':

    dir_to_parse = './test/BRCA_2022'
    result_dir = os.path.join(dir_to_parse, "matched")
    if not os.path.exists(result_dir):
        os.mkdir(os.path.join(dir_to_parse, "matched"))

    raw_files = get_raw_file_list(dir_to_parse)
    raw_to_result = find_matching_result_files(dir_to_parse, raw_files)

    for raw, result in raw_to_result:
        if file_len_match(raw, result):
            raw_lines = open(f'{raw}').readlines()[1:]
            result_lines = open(f'{result}').readlines()
            result_file = os.path.join(
                result_dir, f'matched_{Path(raw).stem}.json')
                
            with open(result_file, 'w') as matched_file:
                dump_records(raw_lines, result_lines, matched_file)
