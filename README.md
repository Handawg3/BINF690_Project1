# Scaffold Filtering and Trimming Pipeline

## Overview
This Python script processes scaffold sequences by:
1. **Filtering** out scaffolds listed in an exclusion file.
2. **Trimming** scaffolds based on a trimming list while retaining sequences above a set threshold.
3. **Generating** processed FASTA files for further analysis.

## Features
- **Exclusion Filtering:** Removes scaffold sequences present in `exclude_list.tab`.
- **Sequence Trimming:** Trims sequences based on specified spans in `trim_list.tab`, with additional handling for mitochondrial sequences.
- **Threshold Enforcement:** Ensures only sequences longer than 200 bases are included in the output.
- **Automated File Handling:** Reads and writes FASTA files efficiently.

## File Structure

## Installation & Requirements
### Prerequisites:
- **Python 3.x**
- No additional libraries required beyond standard Python.

### Running the Script:
1. Place the required files (`Scaffolds.fasta`, `exclude_list.tab`, `trim_list.tab`) in the same directory as the script.
2. Run the script using:
   ```bash
   python project1new.py
