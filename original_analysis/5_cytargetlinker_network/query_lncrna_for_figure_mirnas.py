import requests
import time
import json

# The 8 miRNAs exactly as shown in the current CyTargetLinker/Cytoscape figure
figure_mirnas = [
    "hsa-miR-320a", "hsa-miR-744-5p", "hsa-miR-146a-5p", "hsa-miR-1246",
    "hsa-miR-769-5p", "hsa-miR-654-5p", "hsa-miR-1228-5p", "hsa-miR-451a"
]

results = {}
for mirna in figure_mirnas:
    print(f"=== {mirna} ===")
    url = f"https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38&geneType=lncRNA&miRNA={mirna}&clipExpNum=1&source=all"
    try:
        r = requests.get(url, timeout=20)
        text = r.text.strip()
        if text and text not in ("[]", "") and "error" not in text.lower():
            lines = [l for l in text.split("\n") if l.strip()]
            print(f"  {len(lines)} lncRNA target record(s) found")
            for l in lines[:5]:
                print(f"    {l}")
            results[mirna] = lines
        else:
            print("  No lncRNA targets returned")
            results[mirna] = []
    except Exception as e:
        print(f"  Query failed: {e}")
        results[mirna] = None
    time.sleep(1)

with open("figure_mirna_lncrna_results.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nSaved figure_mirna_lncrna_results.json")
