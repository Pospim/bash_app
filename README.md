# GoLizard

GoLizard is a robust bioinformatics pipeline designed to predict the function and localization of unknown or mutated protein sequences. It leverages the Gene Ontology (GO) framework and integrates multiple inference strategies, including sequence similarity, structural homology, and motif detection, to deliver comprehensive functional annotation.

---

## 🧬 Introduction

The increasing discovery of novel and mutated proteins has created a critical gap in functional annotation. Traditional approaches often rely solely on sequence or structure-based inference, which limits their effectiveness. GoLizard addresses this challenge by combining multiple complementary strategies to enhance GO-based functional interpretation.

---

## 🚀 What is GoLizard?

GoLizard is a multi-strategy bioinformatics pipeline that:

- Predicts functional roles and subcellular localization of proteins.
- Utilizes Gene Ontology (GO) as the central annotation framework.
- Integrates diverse inference methods: sequence, structural, and motif-based.

---

## 🛠️ Core Tools and Methodology

- **BLAST**: Identifies sequence similarity to known proteins for GO transfer.
- **FoldSeek**: Detects structural homology for function prediction when sequence data is insufficient.
- **ELM (Eukaryotic Linear Motif)**: Searches for functional motifs in disordered regions, used as fallback when BLAST results are missing.
- **REVIGO**: Reduces redundancy in GO terms and visualizes key biological insights.

---

## 📥 Input

- FASTA formatted protein sequence.

---

## 📤 Output

- Ranked GO terms with associated evidence scores.
  - Prioritized sources: BLAST, FoldSeek > ELM.
- Link to REVIGO.
- List of possibly similar proteins and their GO annotations.

---

## 💡 Applications

- Functional annotation of proteins in non-model organisms.
- Mutational impact analysis in proteomics.
- Characterization of hypothetical or novel proteins.
- Predictive insights into pathogen protein function.

---

## 🛣️ Roadmap

- Web-based interface for user-friendly access.
- Integration with pathway and interaction databases for expanded context.

---

## 📦 Installation

- Bash required / WSL for Windows
- Install neccessary libraries from requirements.txt, environment.yml

---

## ▶️ Usage

```bash
# Basic usage with default settings (remote BLAST against SwissProt)
./golizard.sh --file tst_input/protein.fasta

# High stringency search with custom output
./golizard.sh --file tst_input/protein.fa --outputdir high_stringency_results \
  --blast-dbs swissprot pdb --blast-min_ident 90.0 --blast-max_eval 1e-10 \
  --fs-min_ident 90.0 --fs-max_eval 1.0 --revigo-cutoff 0.6

# Local-only BLAST with custom parameters and specific REVIGO output
./golizard.sh --file tst_input/protein.fasta --blast-local /path/to/blastdb/swissprot \
  --blast-min_ident 50.0 --blast-kMax 15 --revigo-result jTable jCytoscape
```

---

## 🤝 Contributing

Contributions are welcome! Please open issues or submit pull requests to help improve GoLizard.

---

## 📄 License

_Include your preferred license here (e.g., MIT, GPL-3.0)._

---

## 📫 Contact

_For questions or support, please contact <marek.pospisil@uochb.cas.cz> <jiri.vondrasek@uochb.cas.cz>._  
