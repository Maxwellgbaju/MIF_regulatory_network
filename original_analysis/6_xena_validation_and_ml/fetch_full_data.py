import xenaPython as xena
import json

host = 'https://toil.xenahubs.net'
pheno_dataset = 'TcgaTargetGTEX_phenotype.txt'
expr_dataset = 'TcgaTargetGtex_rsem_gene_tpm'

with open("xena_phenotype_decoded.json") as f:
    metadata = json.load(f)
samples = list(metadata.keys())
print(f"Loaded metadata for {len(samples)} samples")

# Fetch _primary_site and merge in
print("Fetching _primary_site...")
codebooks_raw = xena.field_codes(host, pheno_dataset, ['_primary_site'])
site_codebook = codebooks_raw[0]['code'].split('\t')
values = xena.dataset_probe_values(host, pheno_dataset, samples, ['_primary_site'])
raw_codes = values[1][0]

for i, s in enumerate(samples):
    try:
        idx = int(raw_codes[i])
        metadata[s]['_primary_site'] = site_codebook[idx] if 0 <= idx < len(site_codebook) else None
    except (ValueError, TypeError):
        metadata[s]['_primary_site'] = None

with open("xena_phenotype_full.json", "w") as f:
    json.dump(metadata, f)
print("Saved xena_phenotype_full.json with _primary_site added")

# Fetch MIF expression for all samples
print("\nFetching MIF expression for all samples (this may take a moment)...")
mif_values = xena.dataset_gene_values(host, expr_dataset, samples, ['MIF'])
mif_expr = {samples[i]: mif_values[1][0][i] for i in range(len(samples))}

with open("mif_expression.json", "w") as f:
    json.dump(mif_expr, f)
print(f"Saved MIF expression for {len(mif_expr)} samples")

# Quick sanity check
valid_vals = [v for v in mif_expr.values() if v is not None]
print(f"\nValid (non-null) expression values: {len(valid_vals)} / {len(mif_expr)}")
print(f"Range: min={min(valid_vals):.3f}, max={max(valid_vals):.3f}")
