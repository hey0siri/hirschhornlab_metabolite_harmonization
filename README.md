# Metabolite Harmonization for metabolite GWAS

> **Status:** ⚠️ Work in progress
> This repository contains an evolving metabolite harmonization workflow. Scripts, parameters, and outputs are under active development and should be validated before use in production analyses.

---

## Repository Purpose

This repository implements a multi-stage metabolite harmonization pipeline designed to:

* Reconcile metabolite identifiers and names across heterogeneous sources
* Construct a reusable reference harmonization database
* Apply that database consistently to  cohort-specific metabolite lists

---

## High-Level Workflow

The pipeline is organized into these stages:

```
[Raw source lists]
        ↓
[Reference harmonization database]
        ↓
[Cohort-specific harmonization]
```

Each stage corresponds to one or more scripts in this repository.

---

## Script Overview & Execution Order (Subject to change)

### 1. Cleaned Identifier Construction (Run Once)

**Purpose:**

* Resolve internal duplicates and harmonize within Moore's list

**Typical outputs:**

* Moore's metabolite tables

This step establishes the metabolite lists provided by Dr. Steven Moore to be used in downstream harmonization.

---

### 2. Reference Harmonization Database Construction (Run Once)

**Purpose:**

* Merge metabolite annotations from:

  * COMETS
  * Dr. Steven Moore’s metabolite list
  * MetLinkR (external harmonization)

**Key outputs:**

* `final_merged_metlinkr_filt.Rdata` (reference database)

⚠️ Due to current MetLinkR limitations, this stage may rely on intermediate files (e.g. `mapping_library_COMETS_Moore.xlsx`).

---

### 3. Cohort / Sample Harmonization (Run Repeatedly)

**Purpose:**

* Map a new cohort-specific metabolite list onto the reference database
* Harmonize using a prioritized identifier strategy
* Retain unmatched metabolites for QC

**Typical outputs: (subject to change)**

* `harmonized_output_<cohort>.Rdata`

  * `$harmonized`: matched metabolites
  * `$remaining`: unmatched metabolites

This stage is designed to be reusable across studies.

---

## Harmonization Philosophy

Harmonization proceeds from highest-confidence identifiers to lowest-confidence matches:

1. HMDB IDs
2. KEGG IDs
3. PubChem CIDs
4. InChIKeys
5. ChemSpider IDs
6. Biochemical / metabolite names

## Key Inputs

* COMETS metabolite lists for reference harmonization database (`.RData`)
* Moore metabolite lists for reference harmonization database (`.csv`)
* Cohort-specific metabolite lists (`.xlsx`, `.csv`)

Input schemas may differ; scripts rely on explicit column mapping rather than fixed column names.

---

## Key Outputs

* Reference harmonization database (`final_merged_metlinkr_filt.Rdata`)
* Cohort-specific harmonization results
* Intermediate CSVs for inspection and debugging

---

## Known Limitations

* MetLinkR integration is currently semi-manual
* Many-to-many joins are permitted and require interpretation
* Name-based harmonization is inherently ambiguous

Results should be reviewed in the context of biological plausibility.

---

*Last updated: Draft*
