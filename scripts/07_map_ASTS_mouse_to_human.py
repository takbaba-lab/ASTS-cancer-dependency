#!/usr/bin/env python3
"""ASTRA-Drug Step 07: frozen mouse ASTS → human ortholog mapping.

"""
from pathlib import Path

# ==================== Path / Parameter ====================
PROJECT_DIR = Path(__file__).resolve().parents[1]
INPUT_FILE = PROJECT_DIR / 'results/06b_ASTS_v1/ASTS_v1.0_primary_100genes.tsv'
REFERENCE_DIR = PROJECT_DIR / 'data/reference/step07_ensembl114'
OUTPUT_DIR = PROJECT_DIR / 'results/07_ASTS_human_v1'
BIOMART_URL = 'https://may2025.archive.ensembl.org/biomart/martservice'
ENSEMBL_RELEASE = '114'
DATASET = 'mmusculus_gene_ensembl'
ID_COLUMN = 'gene_id'  # Mouse gene symbol or ENSMUSG identifier; preserve case.
COMPONENT_COLUMN = 'score_component'
WEIGHT_COLUMN = 'score_weight'
EXPECTED_GENES = 100
EXPECTED_PER_COMPONENT = 50
COMPONENT_WEIGHTS = {'DHT_relative_up': 1, 'DHT_relative_down': -1}
TIMEOUT_SECONDS = 180
ATTEMPTS = 3
RETRY_SECONDS = 5
ATTRIBUTES = [
    'ensembl_gene_id', 'external_gene_name',
    'hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_associated_gene_name',
    'hsapiens_homolog_orthology_type', 'hsapiens_homolog_orthology_confidence',
]
REFERENCE_COLUMNS = ['mouse_ensembl_id', 'mouse_symbol', 'human_ensembl_id',
                     'human_symbol', 'orthology_type', 'orthology_confidence']
import argparse
import csv
import hashlib
import io
import json
import platform
import re
import shutil
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from datetime import datetime, timezone


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def write_tsv(path, fields, rows):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames
        if not fields or len(fields) != len(set(fields)):
            raise ValueError(f'Missing or duplicated column name: {path}')
        rows = list(reader)
    if any(None in row or None in row.values() for row in rows):
        raise ValueError(f'Inconsistent number of TSV columns: {path}')
    return fields, rows


def new_directory(path):
    path.mkdir(parents=True, exist_ok=False)


def request(url, data=None):
    for attempt in range(ATTEMPTS):
        try:
            req = urllib.request.Request(url, data=data, headers={'User-Agent': 'ASTRA-Drug-Step07/1.0'})
            with urllib.request.urlopen(req, timeout=TIMEOUT_SECONDS) as response:
                final_url = response.geturl()
                if urllib.parse.urlparse(final_url).hostname != urllib.parse.urlparse(BIOMART_URL).hostname:
                    raise ValueError(f'Redirect to a different host detected: {final_url}')
                return response.read().decode('utf-8'), final_url
        except (OSError, ValueError):
            if attempt + 1 == ATTEMPTS:
                raise
            time.sleep(RETRY_SECONDS)


def query_xml(attributes):
    root = ET.Element('Query', virtualSchemaName='default', formatter='TSV',
                      header='0', uniqueRows='1', datasetConfigVersion='0.6', completionStamp='1')
    dataset = ET.SubElement(root, 'Dataset', name=DATASET, interface='default')
    for attribute in attributes:
        ET.SubElement(dataset, 'Attribute', name=attribute)
    return ET.tostring(root, encoding='unicode')


def parse_biomart(text, width):
    lines = text.splitlines()
    if not lines or lines[-1].strip() != '[success]':
        raise ValueError('Incomplete BioMart response: [success] marker not found. Stopping instead of treating entries as no ortholog.')
    rows = list(csv.reader(lines[:-1], delimiter='\t'))
    if not rows or any(len(row) != width for row in rows):
        raise ValueError('BioMart response is empty or has an inconsistent number of columns')
    if any(not re.fullmatch(r'ENSMUSG\d+', row[0]) for row in rows):
        raise ValueError('Invalid mouse Ensembl ID detected in the BioMart response')
    return rows


def fetch_reference(directory):
    new_directory(directory)
    registry, _ = request(BIOMART_URL + '?type=registry')
    (directory / 'registry.xml').write_text(registry, encoding='utf-8')
    marts = ET.fromstring(registry).iter('MartURLLocation')
    if not any(m.attrib.get('database', '') == f'ensembl_mart_{ENSEMBL_RELEASE}' for m in marts):
        raise ValueError('Ensembl release in the registry does not match the configured release')
    records = []
    for name, attributes in [('inventory', ATTRIBUTES[:2]), ('homology', ATTRIBUTES)]:
        xml = query_xml(attributes)
        (directory / f'{name}.query.xml').write_text(xml, encoding='utf-8')
        raw, final_url = request(BIOMART_URL, urllib.parse.urlencode({'query': xml}).encode())
        (directory / f'{name}.raw.tsv').write_text(raw, encoding='utf-8')
        rows = parse_biomart(raw, len(attributes))
        write_tsv(directory / f'{name}.tsv', REFERENCE_COLUMNS[:len(attributes)],
                  [dict(zip(REFERENCE_COLUMNS, row)) for row in rows])
        records.append({'query': name, 'retrieved_at_utc': now(), 'url': final_url, 'rows': len(rows)})
    metadata = {'source': 'Ensembl BioMart', 'url': BIOMART_URL,
                'ensembl_release': ENSEMBL_RELEASE, 'dataset': DATASET,
                'queries': records, 'sha256': {p.name: sha256(p) for p in sorted(directory.iterdir())}}
    save_json(directory / 'reference_metadata.json', metadata)
    print(f'Reference retrieval completed: {directory}')


def load_reference(directory):
    metadata = json.loads((directory / 'reference_metadata.json').read_text(encoding='utf-8'))
    if (metadata['ensembl_release'], metadata['url'], metadata['dataset']) != (ENSEMBL_RELEASE, BIOMART_URL, DATASET):
        raise ValueError('Reference metadata do not match the configured parameters')
    required = {'registry.xml', 'inventory.tsv', 'homology.tsv', 'inventory.raw.tsv',
                'homology.raw.tsv', 'inventory.query.xml', 'homology.query.xml'}
    if set(metadata['sha256']) != required:
        raise ValueError('Reference checksum list is incomplete')
    for filename, expected in metadata['sha256'].items():
        if sha256(directory / filename) != expected:
            raise ValueError(f'Reference checksum mismatch: {filename}')
    inventory_fields, inventory = read_tsv(directory / 'inventory.tsv')
    homology_fields, homology = read_tsv(directory / 'homology.tsv')
    if inventory_fields != REFERENCE_COLUMNS[:2] or homology_fields != REFERENCE_COLUMNS:
        raise ValueError('Unexpected column names in the reference TSV')
    ids = {r['mouse_ensembl_id'] for r in inventory}
    if not ids or any(r['mouse_ensembl_id'] not in ids for r in homology):
        raise ValueError('Mouse IDs are inconsistent between inventory and homology tables')
    return metadata, inventory, homology


def load_input(path):
    fields, rows = read_tsv(path)
    required = {ID_COLUMN, COMPONENT_COLUMN, WEIGHT_COLUMN}
    if not required.issubset(fields):
        raise ValueError(f'Required columns: {sorted(required)} / input columns: {fields}')
    if len(rows) != EXPECTED_GENES or len({r[ID_COLUMN] for r in rows}) != EXPECTED_GENES:
        raise ValueError('Mouse input must contain exactly 100 unique genes')
    if any(not r[ID_COLUMN] or r[ID_COLUMN] != r[ID_COLUMN].strip() for r in rows):
        raise ValueError('gene_id contains an empty value or leading/trailing whitespace')
    counts = Counter(r[COMPONENT_COLUMN] for r in rows)
    if counts != Counter({c: EXPECTED_PER_COMPONENT for c in COMPONENT_WEIGHTS}):
        raise ValueError(f'Expected a 50/50 up/down composition: {counts}')
    for row in rows:
        if float(row[WEIGHT_COLUMN]) != COMPONENT_WEIGHTS[row[COMPONENT_COLUMN]]:
            raise ValueError('score_weight is inconsistent with component')
    return fields, rows


def build_mapping(inputs, inventory, homology):
    symbols, by_mouse, by_human = defaultdict(set), defaultdict(list), defaultdict(set)
    valid_ids = {r['mouse_ensembl_id'] for r in inventory}
    for row in inventory:
        if row['mouse_symbol']:
            symbols[row['mouse_symbol']].add(row['mouse_ensembl_id'])
    seen = set()
    for row in homology:
        key = tuple(row[c] for c in REFERENCE_COLUMNS)
        if key in seen:
            continue
        seen.add(key)
        by_mouse[row['mouse_ensembl_id']].append(row)
        if row['human_ensembl_id']:
            if not re.fullmatch(r'ENSG\d+', row['human_ensembl_id']):
                raise ValueError('Invalid human Ensembl ID')
            by_human[row['human_ensembl_id']].add(row['mouse_ensembl_id'])
    audit, gene_summary = [], []
    for index, source in enumerate(inputs, 1):
        identifier = source[ID_COLUMN]
        if re.fullmatch(r'ENSMUSG\d+(?:\.\d+)?', identifier):
            stable = identifier.split('.')[0]  # Preserve the original versioned identifier in the source column.
            candidates = {stable} & valid_ids
        else:
            candidates = symbols.get(identifier, set())
        pairs = [r for mid in sorted(candidates) for r in by_mouse[mid] if r['human_ensembl_id']]
        humans = {r['human_ensembl_id'] for r in pairs}
        types = {r['orthology_type'] for r in pairs}
        if not candidates:
            category, reason = 'unresolved_mouse_id', 'not_found_in_reference'
        elif len(candidates) > 1:
            category, reason = 'ambiguous_mouse_id', 'symbol_matches_multiple_mouse_ids'
        elif not pairs:
            category, reason = 'no_ortholog', 'recognized_mouse_without_human_ortholog_in_this_release'
        elif types == {'ortholog_one2one'} and len(humans) == 1:
            human = next(iter(humans))
            if len(by_human[human]) == 1:
                category, reason = '1:1', 'included'
            else:
                category, reason = 'inconsistent_mapping', 'human_has_multiple_mouse_sources'
        elif types == {'ortholog_one2many'}:
            category, reason = '1:many', 'ensembl_one2many_excluded'
        elif types == {'ortholog_many2many'}:
            category, reason = 'many:many', 'ensembl_many2many_excluded'
        else:
            category, reason = 'inconsistent_mapping', 'unexpected_or_conflicting_orthology_type'
        base = {'input_row': index, 'input_gene_id': identifier,
                'score_component': source[COMPONENT_COLUMN], 'score_weight': source[WEIGHT_COLUMN],
                'mapping_class': category, 'primary_included': str(category == '1:1').lower(),
                'decision_reason': reason, 'n_mouse_candidates': len(candidates),
                'n_human_targets': len(humans)}
        gene_summary.append(dict(base))
        for mid in sorted(candidates) or ['']:
            details = by_mouse[mid] or [dict.fromkeys(REFERENCE_COLUMNS, '')]
            for detail in details:
                entry = {**base, **detail, 'mouse_ensembl_id': mid}
                entry['n_mouse_sources_for_human'] = len(by_human[detail['human_ensembl_id']])
                entry.update({'source_' + k: v for k, v in source.items()})
                audit.append(entry)
    human_inputs = defaultdict(set)
    for row in audit:
        if row['primary_included'] == 'true':
            human_inputs[row['human_ensembl_id']].add(row['input_row'])
    duplicates = {i for indices in human_inputs.values() if len(indices) > 1 for i in indices}
    for row in audit + gene_summary:
        if row['input_row'] in duplicates:
            row.update(primary_included='false', decision_reason='duplicate_human_across_input_rows')
    return audit, gene_summary


def map_signature(input_path, reference_dir, output_dir):
    original_hash = sha256(input_path)
    fields, inputs = load_input(input_path)
    metadata, inventory, homology = load_reference(reference_dir)
    audit, per_gene = build_mapping(inputs, inventory, homology)
    primary = []
    seen = set()
    for row in audit:
        if row['primary_included'] == 'true' and row['human_ensembl_id'] not in seen:
            primary.append(row)
            seen.add(row['human_ensembl_id'])
    if not primary:
        raise ValueError('The primary human signature contains zero genes. Check the reference and input identifiers.')
    new_directory(output_dir)
    shutil.copytree(reference_dir, output_dir / 'reference_snapshot')
    shutil.copy2(input_path, output_dir / 'mouse_primary_input.tsv')
    shutil.copy2(Path(__file__), output_dir / '07_script_snapshot.py')
    write_tsv(output_dir / 'ASTS_v1.0_mouse_human_mapping.tsv', list(audit[0]), audit)
    write_tsv(output_dir / 'ASTS_v1.0_mouse_mapping_status.tsv', list(per_gene[0]), per_gene)
    write_tsv(output_dir / 'ASTS_human_v1.0_primary.tsv', list(audit[0]), primary)
    counts = Counter(r['mapping_class'] for r in per_gene)
    included = Counter(r['score_component'] for r in primary)
    summary = ['ASTRA-Drug Step 07 summary', f'Mouse input: {len(inputs)} (up=50, down=50)',
               f'Ensembl release: {ENSEMBL_RELEASE}', f'Primary human genes: {len(primary)}']
    for category in ['1:1', '1:many', 'many:many', 'no_ortholog', 'unresolved_mouse_id',
                     'ambiguous_mouse_id', 'inconsistent_mapping']:
        summary.append(f'Mouse classification {category}: {counts[category]}')
    for component in COMPONENT_WEIGHTS:
        n = included[component]
        summary.append(f'Primary {component}: {n}/50 retained ({n / 50:.1%})')
    summary += [f'Exclusion reasons: {dict(Counter(r["decision_reason"] for r in per_gene if r["primary_included"] == "false"))}',
                'No replacement, re-ranking, direction change, or rebalancing.',
                'no_ortholog means absent in this reference release, not proven biological absence.',
                'Human response is projected from mouse, not experimentally validated in human.',
                'ASTS represents a DHT-induced relative transcriptional state; not androgen-specific.']
    (output_dir / '07_summary.txt').write_text('\n'.join(summary) + '\n', encoding='utf-8')
    if sha256(input_path) != original_hash or sha256(output_dir / 'mouse_primary_input.tsv') != original_hash:
        raise ValueError('Input changed during processing. The freeze will not be created.')
    manifest = {'status': 'FROZEN', 'created_at_utc': now(), 'input_path': str(input_path.resolve()),
                'mouse_input_sha256': original_hash, 'mouse_primary_modified': False,
                'reference': metadata, 'python': platform.python_version(),
                'platform': platform.platform(), 'command': sys.argv,
                'primary_human_genes': len(primary), 'primary_by_component': dict(included),
                'rules': ['Ensembl ortholog_one2one only; unique mouse ID and human ID',
                          'Exclude ambiguous IDs, conflicting types and duplicate human targets',
                          'Orthology confidence recorded, not used as extra selection filter',
                          'Exact case-sensitive symbol match; no alias guessing',
                          'Keep original input order, rank, component and weight',
                          'No replacement or up/down rebalancing; mouse primary stays frozen',
                          'No DepMap/PRISM/GDSC data used in this script'],
                'sha256': {str(p.relative_to(output_dir)): sha256(p)
                           for p in sorted(output_dir.rglob('*')) if p.is_file()}}
    save_json(output_dir / 'ASTS_human_v1.0_FREEZE_MANIFEST.json', manifest)
    files = [p for p in sorted(output_dir.rglob('*')) if p.is_file()]
    (output_dir / 'ASTS_human_v1.0_SHA256SUMS.txt').write_text(
        ''.join(f'{sha256(p)}  {p.relative_to(output_dir)}\n' for p in files), encoding='utf-8')
    print('\n'.join(summary))
    print(f'Completed successfully; freeze created: {output_dir}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['fetch', 'map'])
    parser.add_argument('--input', type=Path, default=INPUT_FILE)
    parser.add_argument('--reference-dir', type=Path, default=REFERENCE_DIR)
    parser.add_argument('--output-dir', type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    if args.mode == 'fetch':
        fetch_reference(args.reference_dir)
    else:
        map_signature(args.input, args.reference_dir, args.output_dir)


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
