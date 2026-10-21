import xenaPython as xena
import json

host = 'https://toil.xenahubs.net'
dataset = 'TcgaTargetGTEX_phenotype.txt'

print("Fetching sample list...")
samples = xena.dataset_samples(host, dataset, None)
print(f"Total samples in cohort: {len(samples)}")

fields = ['_study', '_sample_type', '_primary_site', 'primary disease or tissue']
print(f"Fetching fields: {fields}")

values = xena.dataset_probe_values(host, dataset, samples, fields)

# values is a list of lists (one list per field, same order as `fields`), each aligned with `samples`
result = {}
for i, s in enumerate(samples):
    result[s] = {fields[j]: values[j][i] for j in range(len(fields))}

with open("xena_phenotype.json", "w") as f:
    json.dump(result, f)

print(f"Saved metadata for {len(result)} samples to xena_phenotype.json")

# Quick sanity check: show unique values for _study and _sample_type
studies = set(v['_study'] for v in result.values())
sample_types = set(v['_sample_type'] for v in result.values())
print(f"\nUnique _study values: {studies}")
print(f"\nUnique _sample_type values: {sample_types}")
