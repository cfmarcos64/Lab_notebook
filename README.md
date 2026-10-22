# 🧪 Lab Reaction Suite & Scheme Exporter

![Python](https://img.shields.io/badge/python-3.10+-blue)
![License](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-green)
![Streamlit](https://img.shields.io/badge/Streamlit-Cloud-orange)

A comprehensive **Streamlit** suite designed for organic chemists, researchers, and students to seamlessly draw chemical reactions, calculate stoichiometry, query PubChem data, generate laboratory notebooks, and export high-resolution publication-ready reaction schemes.

All tools share a unified design philosophy: draw reactions interactively using **Ketcher**, leverage **RDKit** and **PubChem API** for chemical intelligence, and export formatted Word documents or high-DPI images for reports, papers, and theses.

---

## 🚀 Included Modules

1. **Reaction Scheme & Molecule Image Exporter (`Reaction_Scheme_Exporter.py`)**  
   * Generates publication-ready chemical reaction schemes and molecule grid figures with precise physical measurements (cm, mm, pt, DPI).
   * Supports automatic or manual molecule labeling, R-group generic labels ($R_1, R_2$, etc.), reaction conditions display, and high-resolution PNG export (300, 600, 1200 DPI).

2. **Reaction Stoichiometry & Lab Notebook (`Reaction_Lab_Notebook.py`)**  
   * Automatically parses drawn reactions into individual reactants and products.
   * Fetches official IUPAC names via InChIKey querying against PubChem REST API.
   * Calculates stoichiometry dynamically based on user-defined reference limiting reagents, densities, and solvent volumes.
   * Exports fully formatted laboratory reports as Word documents (`.docx`) including automated pagination and formatted headers.

---

## ✨ Key Features

- **Interactive Drawing Canvas**: Integrated with `streamlit-ketcher` for intuitive drawing of single molecules or complete multi-step reaction schemes.
- **Automated Chemical Intelligence**: InChIKey generation, SMILES parsing, and IUPAC name resolution powered by RDKit and PubChem.
- **Publication-Ready Figures**: Physical-unit rendering engine ensuring exact bond lengths, font scaling, and line weights matching journal requirements (e.g., ACS style).
- **Dynamic Lab Notebook Reports**: Live stoichiometric table generation with customizable reactant equivalence, density, and physical state settings.
- **Seamless Document Export**: One-click generation of `.docx` lab reports and publication-grade `.png` image schemes.

---

## 🛠️ How to Run the Suite Locally

1. **Clone the Repository & Install Python Dependencies**:

```bash
git clone https://github.com/cfmarcos64/Lab_notebook.git
cd Lab_notebook

pip install -r requirements.txt
```

2. **Launch the Suite with Streamlit**:

```bash
streamlit run Chem_tool.py
```

---

## 📦 Project Dependencies

### `requirements.txt`
```text
streamlit
pandas
rdkit
streamlit-ketcher
python-docx
Pillow
```

### `packages.txt` (For Linux / Streamlit Cloud RDKit & Cairo Stability)
```text
libxrender1
libxext6
libfontconfig1
```

---

## 📁 Repository Structure

The repository is organized into a clean, modular structure:

```text
Reaction_Lab_Suite/
├── Chem_tool.py                   # Master navigation script & sidebar theme controller
├── Reaction_Scheme_Exporter.py    # Reaction scheme & molecule figure renderer module
├── Reaction_Lab_Notebook.py      # Reaction stoichiometry & Word report generator module
├── requirements.txt               # Python package dependencies
├── packages.txt                   # Linux system dependencies
└── README.md                      # Project documentation
```

---

## 📄 Licence

This project is licensed under the CC BY-NC-SA 4.0 Licence. See the [LICENSE](LICENSE) file for further details.

---

## 👤 Author & Contact
* **Created by**: Carlos Fernández Marcos
* **Institution**: Universidad de Extremadura (University of Extremadura, Spain)
* **GitHub**: [@cfmarcos64](https://github.com/cfmarcos64)
* **Email**: cfmarcos@unex.es