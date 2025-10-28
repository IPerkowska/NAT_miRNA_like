#!/usr/bin/env python3
"""
generate_RNA_structures_VARNA.py
----------------------------------------------------------
Generates publication-quality RNA secondary structure images
for selected candidates using VARNA (v3.93).

Input:
    NAT_noncano_candidates.xlsx  – table with sequences and group_id column

Output:
    ./RNA_structures_VARNA/  – PNG files (per candidate × temperature)
"""

import os
import pandas as pd
import subprocess

# === CONFIGURATION ===========================================================
EXCEL_FILE = "NAT_noncano_candidates.xlsx"   # input file
SELECTED_IDS = [229, 253, 267, 275, 283, 289]      # candidate IDs to visualize
TEMPS = [4, 22, 37]                                # temperatures to fold
OUT_DIR = "RNA_structures_VARNA"                   # output folder
VARNA_JAR = "VARNAv3-93.jar"                       # path to VARNA JAR file
# ============================================================================

os.makedirs(OUT_DIR, exist_ok=True)

print(f"📂 Loading candidate sequences from {EXCEL_FILE} ...")
df = pd.read_excel(EXCEL_FILE)

# Filter only selected candidates
df = df[df["group_id"].isin(SELECTED_IDS)]
if df.empty:
    raise ValueError("❌ None of the selected group_ids were found in the Excel file.")

print(f"🧬 Generating VARNA structures for candidates: {', '.join(map(str, SELECTED_IDS))}")

# === Helper: Run RNAfold and generate .dbn files =============================
def run_rnafold(sequence, temperature):
    """Run RNAfold and return (structure, MFE)."""
    sequence = sequence.upper().replace("T", "U")
    cmd = ["RNAfold", "--noPS", f"--temp={temperature}"]
    result = subprocess.run(cmd, input=sequence, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"RNAfold failed at {temperature}°C: {result.stderr}")
    lines = result.stdout.strip().split("\n")
    struct_line = lines[-1]
    struct, mfe = struct_line.split(" (")
    mfe = mfe.replace(")", "").strip()
    struct = struct.strip()
    return struct, mfe, sequence

# === Helper: Generate PNG with VARNA ========================================
def run_varna(sequence, structure, out_path):
    """Use VARNA to render RNA structure to PNG."""
    cmd = [
        "java", "-cp", VARNA_JAR,
        "fr.orsay.lri.varna.applications.VARNAcmd",
        "-sequenceDBN", sequence,
        "-structureDBN", structure,
        "-algorithm", "naview",
        "-resolution", "8.0",
        "-o", out_path,
        "-bpStyle", "line",
        "-bpColor", "#0000FF",
        "-unpairedColor", "#B0B0B0",
        "-noTitle", "true"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0:
        print(f"✅ {os.path.basename(out_path)} generated")
    else:
        print(f"⚠️ VARNA failed for {os.path.basename(out_path)}:\n{result.stderr}")

# === Main loop ==============================================================
for _, row in df.iterrows():
    gid = int(row["group_id"])
    seq = row["seq"]

    for temp in TEMPS:
        struct, mfe, rna_seq = run_rnafold(seq, temp)
        out_png = os.path.join(OUT_DIR, f"Candidate_{gid}_{temp}C.png")

        run_varna(rna_seq, struct, out_png)

print(f"\n🎉 All PNG structures saved in: {OUT_DIR}/")
print("🧾 Each candidate folded with RNAfold (4°C, 22°C, 37°C) and visualized via VARNA (naview layout).")
