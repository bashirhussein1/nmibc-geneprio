# Validation Guide for nmibc-geneprio Pipeline

This guide describes how to perform script-level validation for each module of the nmibc-geneprio pipeline.

## Required tools
- Apptainer/Singularity installed
- Container built: `apptainer build containers/pipeline.sif containers/singularity.def`
- LDREF data available in `LDREF/` directory

---

## Module 1: eQTL Lookup
```bash
apptainer exec containers/pipeline.sif Rscript scripts/validation/eqtl_lookup_simulate.R
```
Expected output: `Number of overlapping variants found: 6`

---

## Module 2: Colocalization
```bash
apptainer exec containers/pipeline.sif Rscript scripts/validation/coloc_simulate.R
```
Expected output: PP.H4 ≥ 0.8

---

## Module 3: TWAS

### Step 1 — Download required utility
```bash
wget https://raw.githubusercontent.com/gusevlab/fusion_twas/master/utils/plink_utils.R \
  -O scripts/utils/plink_utils.R
```

### Step 2 — Prepare LDREF without chr prefix
```bash
mkdir -p LDREF_sim
for chr in 21 22; do
    cp LDREF/1000G.EUR.$chr.bed LDREF_sim/1000G.EUR.$chr.bed
    cp LDREF/1000G.EUR.$chr.fam LDREF_sim/1000G.EUR.$chr.fam
    awk '{sub(/^chr/, "", $1); sub(/^chr/, "", $2); print}' OFS='\t' \
        LDREF/1000G.EUR.$chr.bim > LDREF_sim/1000G.EUR.$chr.bim
done
```

### Step 3 — Run validation
```bash
apptainer exec containers/pipeline.sif Rscript scripts/validation/fusion_twas_validation.R
```

Expected output:
- chr22: 18/20 causal genes significant (sensitivity = 90%)
- chr21: 0/10 null genes significant (false positive rate = 0%)
