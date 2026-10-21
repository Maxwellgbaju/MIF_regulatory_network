import xenaPython as xena
import json

host = 'https://toil.xenahubs.net'
dataset = 'TcgaTargetGTEX_phenotype.txt'

print("Fetching sample list...")
samples = xena.dataset_samples(host, dataset, None)
print(f"Total samples: {len(samples)}")

fields = ['_study', '_sample_type', 'primary disease or tissue']

print("Fetching codebooks...")
codebooks_raw = xena.field_codes(host, dataset, fields)
codebooks = {c['name']: c['code'].split('\t') for c in codebooks_raw}

print("Fetching raw coded values for all samples...")
values = xena.dataset_probe_values(host, dataset, samples, fields)
raw_codes = values[1]  # values[0] is None (coding handled separately via field_codes)

result = {}
for i, s in enumerate(samples):
    entry = {}
    for j, field in enumerate(fields):
        code_idx = raw_codes[j][i]
        try:
            code_idx = int(code_idx)
            label = codebooks[field][code_idx] if 0 <= code_idx < len(codebooks[field]) else None
        except (ValueError, TypeError):
            label = None
        entry[field] = label
    result[s] = entry

with open("xena_phenotype_decoded.json", "w") as f:
    json.dump(result, f)

print(f"Saved decoded metadata for {len(result)} samples")

# Sanity check
from collections import Counter
study_counts = Counter(v['_study'] for v in result.values())
sample_type_counts = Counter(v['_sample_type'] for v in result.values())
print(f"\n_study distribution: {dict(study_counts)}")
print(f"\n_sample_type distribution (top 10): {sample_type_counts.most_common(10)}")
