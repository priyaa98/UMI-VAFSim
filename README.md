# UMIVAFSim (UMI-Aware Variant Allele Fraction Simulator)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![pysam](https://img.shields.io/badge/dependency-pysam-green.svg)](https://github.com/pysam-developers/pysam)

**UMIVAFSim** is a data-driven, empirical in-silico variant spike-in simulator designed for high-depth Next-Generation Sequencing (NGS) and liquid biopsy (ctDNA/cfDNA) assays. 

Unlike traditional *de novo* read simulators that generate artificial reads from reference FASTA files, **UMIVAFSim** modifies real sequencing alignments directly. It incorporates **empirical Phred error profiling**, **coherent UMI molecular family mutations**, **Negative Binomial amplification dispersion**, and **asymmetric strand bias** to generate ultra-realistic synthetic variants at precise Variant Allele Fractions (VAFs).

---

## Highlights & Key Features

* **Preserves Native Backgrounds:** Modifies real BAM alignments in place, retaining original sequencing artifacts, real coverage fluctuations, and background error profiles.
* **Coherent UMI Family Mutations:** Groups reads by Unique Molecular Identifier (`MI` or `RX` tags) to mutate entire molecular families simultaneously, preventing single-read artificial variants.
* **Empirical Error Profiling:** Scans flanking non-variant background reads within a configurable window (e.g., $\pm 5000\text{ bp}$) to recalibrate Phred quality scores and capture cycle-dependent error decay.
* **Negative Binomial Family Boosting:** Simulates non-uniform PCR duplication bias and family size dispersion by scaling molecular family depths via Negative Binomial sampling ($\mu, \phi$).
* **Optional BED Panel Annotation & Validation:** Queries custom target BED files to ensure spike-in sites fall within defined target regions and reports distance-to-edge metrics.
* **Local COSMIC SQLite Engine:** Automatically builds and indexes an offline SQLite database from raw COSMIC mutation archives for fast locus and variant lookup.

---

## Architecture Overview

```mermaid
graph TD
    subgraph Pass0 ["Pass 0: Empirical Error Profiling Engine"]
        P0A["Scan local genomic window"] --> P0B["Build Phred recalibration & cycle error profiles"]
    end

    subgraph Pass1 ["Pass 1: UMI Family Scan & Selection"]
        P1A["Group alignments via MI/RX tags"] --> P1B["Stochastically select target families (Binomial)"]
        P1B --> P1C["Calculate family depth multipliers (Negative Binomial)"]
    end

    subgraph Pass2 ["Pass 2: Global FASTQ Stream"]
        P2A["Inject allele & recalibrate Phred scores"] --> P2B["Stream synchronized paired-end FASTQs"]
    end

    Pass0 --> Pass1
    Pass1 --> Pass2
```


---

## Installation

### Prerequisites
* Python $\ge 3.8$
* `pysam`
* `numpy`

### Setup
```bash
# Clone repository
git clone [https://github.com/your-username/UMIVAFSim.git](https://github.com/your-username/UMIVAFSim.git)
cd UMIVAFSim

# Install requirements
pip install -r requirements.txt
