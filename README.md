

# 🐍 Vritra

**A Scalable Pipeline for Functional Gene Detection and Species Attribution from Metagenomic and Metatranscriptomic Data**

---

## 📖 Overview

**Vritra** is a lightweight and scalable pipeline designed for:

* Functional gene detection
* Species-level attribution
* Integration of metagenomic and metatranscriptomic data

This repository contains the implementation used in the paper:

> *“A Scalable Pipeline for Functional Gene Detection and Species Attribution from Metagenomic and Metatranscriptomic Data”*

---

## ⚙️ Requirements

### 🔧 External Tools

Please install the following tools before running Vritra:

* **DIAMOND2** (tested with v2.0.15)
  [https://github.com/bbuchfink/diamond](https://github.com/bbuchfink/diamond)
* **Kraken2**
* **Bowtie2**
* **samtools**

**Note:** Make sure all tools are added to your PATH environment variable.

---

### 🐍 Python Environment

Vritra was developed and tested with:

* Python **3.6.5**

Required Python packages:

```
pip install numpy pandas biopython ete3
```

| Package   | Tested Version            |
| --------- | ------------------------- |
| numpy     | 1.19.5                    |
| pandas    | 1.1.5                     |
| biopython | latest stable recommended |
| ete3      | latest stable recommended |

**Tip:** If you encounter compatibility issues, consider recreating this environment.

---

## 🚀 Installation

```
git clone https://github.com/BoyanZhou/Vritra.git
cd Vritra
```

No additional compilation is required.

---

## ▶️ Usage

Detailed documentation and examples are available in the Wiki:

[https://github.com/BoyanZhou/Vritra/wiki](https://github.com/BoyanZhou/Vritra/wiki)

---

## 📬 Support & Contact

For questions, bug reports, or suggestions:

[boyanzhou1992@gmail.com](mailto:boyanzhou1992@gmail.com)

---

## 🧩 Notes

* Designed for scalability on large sequencing datasets
* Compatible with both metagenomic and metatranscriptomic workflows
* Minimal dependencies for easier deployment

---

## ⭐ Citation

If you find this tool useful, please consider citing:

[Add your paper citation here]
