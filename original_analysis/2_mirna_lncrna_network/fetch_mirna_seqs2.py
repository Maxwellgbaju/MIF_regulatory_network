import requests
import time

targets = ["hsa-miR-1228-5p", "hsa-miR-23a-5p", "hsa-miR-363-5p",
           "hsa-miR-23b-5p", "hsa-miR-1537-3p", "hsa-miR-92a-1-5p"]

records = []
for mir in targets:
    query = {
        "query": f'"{mir}" AND expert_db:"mirbase"',
        "fields": "description",
        "format": "json",
    }

    data = None
    for attempt in range(3):
        try:
            r = requests.get("https://www.ebi.ac.uk/ebisearch/ws/rest/rnacentral", params=query, timeout=15)
            data = r.json()
            break
        except Exception as e:
            print(f"{mir}: attempt {attempt+1} failed ({e}), retrying...")
            time.sleep(2)

    if data is None:
        print(f"{mir}: FAILED after 3 attempts, skipping")
        continue

    entries = data.get("entries", [])
    if not entries:
        print(f"{mir}: NO MATCH FOUND")
        continue

    urs = entries[0]["id"].split("_")[0]
    print(f"{mir}: matched {len(entries)} entr(y/ies), using {urs}")

    try:
        seq_r = requests.get(f"https://rnacentral.org/api/v1/rna/{urs}.fasta", timeout=15)
        fasta = seq_r.text.strip()
        print(fasta)
        records.append((mir, urs, fasta))
    except Exception as e:
        print(f"{mir}: sequence fetch failed ({e})")

    time.sleep(1)  # slightly longer pause between queries

with open("candidate_miRNA_sequences_verified.fa", "w") as f:
    for mir, urs, fasta in records:
        lines = fasta.splitlines()
        seq = "".join(lines[1:])
        f.write(f">{mir}|{urs}\n{seq}\n")

print(f"\nSaved {len(records)}/{len(targets)} sequences to candidate_miRNA_sequences_verified.fa")
