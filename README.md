# RBPscan: RNA-Binding Protein Motif Scanner

**RBPscan** is a Python-based tool designed to identify and analyze RNA-binding protein (RBP) motifs from sequence data. This tool enables researchers to filter, count, and analyze sequences containing specific motifs, with a focus on motifs involved in RNA editing. The tool supports the analysis of sequences containing recorder motifs, identifies insertions, and calculates editing efficiencies.


## Purpose

The primary purpose of RBPscan is to facilitate the exploration of RNA editing events and their associated motifs, providing insights into the editing efficiency and expression of RNA molecules. The tool processes sequencing data, identifies relevant motifs, and generates detailed statistics for further analysis.


## Files in This Repository

### 1. **RBPscan\_editing\_counting\_Algorithm.txt**

- **Description**: An annotated guide to the algorithm implemented in RBPscan. It provides a comprehensive walkthrough of the steps involved in filtering, motif identification, and calculation of editing statistics.

- **Key Features**:

  - Filtering sequences containing recorder motifs.
  - Counting and summarizing editing events.
  - Generating and saving data frames with editing statistics.


### 2. **recorder\_7N\_and\_EMPTY\_Pumilio.py**

- **Description**: A Python script focused on analyzing sequences with 7N or EMPTY motifs. It identifies motifs with insertions and calculates various statistics like editing efficiencies and occurrence frequencies.

- **Key Features**:

  - Filtering sequences with specific 5' ends and 7N or EMPTY motifs.
  - Generating density plots and exporting motif statistics to CSV files.
  - Visualization support using Seaborn and Matplotlib.


### 3. **RBPscan\_editing\_counting.py**

- **Description**: The main script for RBPscan, enabling motif filtering and analysis. This script is highly customizable for different insertion criteria.

- **Key Features**:

  - Identification of recorder motifs with and without insertions.
  - Calculation of expression percentages and editing efficiencies.
  - Export of results as a CSV file for downstream analysis.


## Usage

1. Install the required Python libraries:

       bash
       Copy code
       pip install seaborn matplotlib numpy pandas

2. Prepare your sequence data in a text file (`YOUR_FILE.txt`).

3. Run the desired script:

       bash
       Copy code
       python3 RBPscan_editing_counting.py

   Replace `YOUR_FILE.txt` with the path to your input file.

4. Modify the regular expressions and cutoff parameters as necessary to tailor the analysis to your data.


## Output

- **CSV files**: Contain statistics such as motif occurrences, editing efficiencies, and percentage of edited reads.
- **Visualizations**: Generate density plots and other charts to illustrate motif distributions and editing metrics.


## License

This project is released under the MIT License. See `LICENSE` for details.
