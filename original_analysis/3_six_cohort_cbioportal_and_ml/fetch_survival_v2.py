import requests, csv, time

studies = ["kirc_tcga_pan_can_atlas_2018","lihc_tcga_pan_can_atlas_2018",
           "stad_tcga_pan_can_atlas_2018","brca_tcga_pan_can_atlas_2018",
           "lusc_tcga_pan_can_atlas_2018","skcm_tcga_pan_can_atlas_2018"]

rows = []
for study in studies:
    url = f"https://www.cbioportal.org/api/studies/{study}/clinical-data?clinicalDataType=PATIENT&projection=SUMMARY"
    data = None
    for attempt in range(4):
        try:
            r = requests.get(url, timeout=30)
            if r.status_code != 200:
                print(f"  {study}: attempt {attempt+1} got HTTP {r.status_code}, retrying...")
                time.sleep(3)
                continue
            data = r.json()
            break
        except Exception as e:
            print(f"  {study}: attempt {attempt+1} failed ({e}), retrying...")
            time.sleep(3)

    if data is None:
        print(f"  {study}: FAILED after 4 attempts, skipping")
        continue

    surv = [d for d in data if d["clinicalAttributeId"] in ("OS_MONTHS", "OS_STATUS")]
    print(f"{study}: {len(surv)} survival-field records")
    for d in surv:
        rows.append({"studyId": study, "patientId": d["patientId"],
                     "attribute": d["clinicalAttributeId"], "value": d["value"]})
    time.sleep(1)

with open("survival_6cancers_raw.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["studyId","patientId","attribute","value"])
    w.writeheader()
    w.writerows(rows)
print(f"\nSaved {len(rows)} rows")
