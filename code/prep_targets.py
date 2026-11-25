import os, subprocess, sys
DATA = os.path.expanduser("~/invfold_survey/data")
IDS = """1CRN 1PGB 1SHG 2IGD 1UBQ 2CI2 1LMB 1TIT 1TEN 1RIS 1FKB 1BNI 1MJC 3CHY 2RN2
2LZM 1AKE 1BTA 1CSP 1PHT 1RUN 2PTL 1QYS 1SRL 1TCA 1WHI 2AYH 1AYE 1CEX 3PGK
1BDD 2GB1 1MBN 1HRC 1UTG 1SN3 1PGX 1OPD 1FYN 1PLQ""".split()
def ca_count(path, chain="A"):
    n = 0
    for line in open(path, errors="ignore"):
        if line.startswith("ATOM") and line[21] == chain and line[12:16].strip() == "CA":
            n += 1
    return n
ok = []
for pid in IDS:
    p = os.path.join(DATA, f"{pid}.pdb")
    if not os.path.exists(p) or os.path.getsize(p) < 2000:
        r = subprocess.run(["curl", "-sL", "-m", "60", "-o", p,
                            f"https://files.rcsb.org/download/{pid}.pdb"], capture_output=True)
    n = ca_count(p) if os.path.exists(p) else 0
    if 50 <= n <= 450:
        ok.append(pid)
        print(f"{pid} L={n} ok", flush=True)
    else:
        print(f"{pid} L={n} SKIP", flush=True)
open(os.path.join(DATA, "targets30.txt"), "w").write("\n".join(ok) + "\n")
print("VALID", len(ok), flush=True)
