
# FSPBWT

## Overview
FSPBWT is a C++ program designed for processing and querying genomic data stored in VCF (Variant Call Format) files. It supports fuzzy panel creation and querying with configurable parameters for in-panel and out-panel queries. The program allows saving and loading preprocessed panels for efficient querying and supports different syllable sizes (64 or 128 bits) and fuzzy levels.

## Features
- **Input/Output**: Reads VCF files for panel and query data, outputs results to text files.
- **Fuzzy Panel Creation**: 
  - **Global (Count-based)**: Creates fuzzy panels based on variant counts across the dataset.
  - **Even (Position-based)**: Creates fuzzy panels based on evenly distributed positions in the dataset.
- **Query Modes**: 
  - **In-panel queries**: Matches within the panel data.
  - **Out-panel queries**: Matches query data against the panel.
- **Configurable Parameters**:
  - `B`: Syllable size (64 or 128 bits).
  - `F`: Fuzzy level (1 to 4).
  - `L`: Minimum match length for queries.
  - `even`: Toggle for position-based fuzzy panel creation.
- **Save/Load**: Save processed panels to binary files and load them for faster queries.
- **Static Linking**: Ensures portability with static binary compilation.

## Prerequisites
- **CMake**: Version 3.28 or higher.
- **C++ Compiler**: Supporting C++14 standard (e.g., g++).
- **Operating System**: Compatible with Linux/Unix systems (tested on Ubuntu).

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/BioLab310/FSPBWT.git
   cd FSPBWT
   ```
2. Create a build directory and generate build files:
   ```bash
   mkdir build && cd build
   cmake ..
   ```
3. Compile the program:
   ```bash
   make
   ```

## Usage
Run the program with the following command-line options:
```bash
./FSPBWT [options]
```

### Options
- `-h, -H`: Print help message and exit.
- `-B, -b <value>`: Set syllable size (default: 64, options: 64 or 128).
- `-F, -f <value>`: Set fuzzy level (default: 1, range: 1–4).
- `-L, -l <value>`: Set minimum match length .
- `-i, -I <file>`: Specify input panel file .
- `-o, -O <file>`: Specify output file (default: auto-generated based on parameters).
- `-m, -M <mode>`: Query mode (`in` for in-panel, `out` for out-panel).
- `-q, -Q <file>`: Specify query file for out-panel mode .
- `-e, -even`: Enable position-based fuzzy panel creation (default: false, uses count-based mode).

### Example Commands
1. Run an in-panel query with custom parameters:
   ```bash
   ./FSPBWT -m in -i panel.vcf -q query.vcf -b 64 -f 1 -l 1000
   ```
2. Run an out-panel query with custom parameters:
   ```bash
   ./FSPBWT -m out -B 128 -F 2 -L 2000 -i panel.vcf -q query.vcf 
   ```

## Output
- **Query Results**: Written to the specified output file (or auto-generated name).
- **Information File**: Additional information about the query is saved to a file prefixed with `information_` (e.g., `information_panel_inPanelQuery_B_F_L.txt`).

## Notes
- Ensure that the input VCF files (`panel.vcf`, `query.vcf`) are correctly formatted.
- The program supports syllable sizes of 64 or 128 bits; other values will result in an error.
- Fuzzy level (`F`) must be between 1 and 4.
- The program is statically linked for portability but may require significant memory for large VCF files.

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contributing
Contributions are welcome! Please submit a pull request or open an issue for bug reports or feature requests.

## Contact
For questions or support, contact [your email] or open an issue on GitHub.
