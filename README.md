# SAM Parser

An Elixir library for parsing and manipulating SAM (Sequence Alignment/Map) and BAM (Binary Alignment/Map) format files according to the SAM format specification v1.6.


## Overview

SAM Parser provides a convenient and efficient way to work with genomic alignment data in both SAM and BAM formats. This library offers comprehensive functionality for reading, parsing, manipulating, and writing alignment data, making it suitable for bioinformatics workflows and genomic data analysis.

## Features

- Fast parsing of SAM and BAM files
- Complete header and alignment record parsing
- Support for all standard and optional fields
- Simple filtering capabilities for alignments based on reference and position
- CIGAR string parsing for alignment operations
- SAM file writing capability
- Comprehensive FLAG field interpretation (is_paired?, is_mapped?, etc.)

## Installation

You can install the package directly from GitHub:

```elixir
def deps do
  [
    {:sam_parser, git: "https://github.com//tripitakit/sam_parser.git"}
  ]
end
```

Then run:

```
$ mix deps.get
```


## License

This project is licensed under the GNU General Public License v3.0 - see the [COPYING](COPYING) file for details.


## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the project
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Contact

Author: Patrick De Marta  
Email: patrick.demarta@gmail.com
