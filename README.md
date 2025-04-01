<img src="logo.png" alt="SAM Parser Logo" width="200"/>

# SAM Parser

A robust Elixir library for parsing and manipulating SAM (Sequence Alignment/Map) and BAM (Binary Alignment/Map) format files according to the SAM format specification v1.6.

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

Add `sam_parser` to your list of dependencies in `mix.exs`:

```elixir
def deps do
  [
    {:sam_parser, "~> 0.1.0"}
  ]
end
```

Then run:

```
$ mix deps.get
```

## Usage

### Parsing SAM/BAM Files

```elixir
# Parse any SAM/BAM file based on extension
sam_file = SamParser.parse_file("path/to/alignment.sam")

# Or parse specific format
sam_file = SamParser.parse_sam("path/to/alignment.sam")
bam_file = SamParser.parse_bam("path/to/alignment.bam")
```

### Accessing Data

```elixir
# Get header information
header = sam_file.header
reference_sequences = SamParser.reference_sequences(sam_file)

# Access all alignments
alignments = sam_file.alignments

# Access a single alignment
first_alignment = Enum.at(alignments, 0)
IO.puts("Read name: #{first_alignment.qname}")
IO.puts("Mapped to: #{first_alignment.rname} at position: #{first_alignment.pos}")
IO.puts("Sequence: #{first_alignment.seq}")
```

### Filtering Alignments

```elixir
# Filter alignments by reference sequence
chr1_alignments = SamParser.filter_by_reference(sam_file, "chr1")

# Filter alignments by position
region_alignments = SamParser.filter_by_position(sam_file, 1000, 2000)
```

### Flag Operations

```elixir
# Check alignment flags
alignment = Enum.at(sam_file.alignments, 0)
if SamParser.is_mapped?(alignment) do
  IO.puts("Read is mapped")
end

if SamParser.is_paired?(alignment) do
  IO.puts("Read is paired")
end

if SamParser.is_properly_paired?(alignment) do
  IO.puts("Read is properly paired")
end

if SamParser.is_reverse?(alignment) do
  IO.puts("Read is reverse complemented")
end
```

### CIGAR Operations

```elixir
# Parse CIGAR string
cigar_operations = SamParser.parse_cigar(alignment.cigar)

# Process CIGAR operations
Enum.each(cigar_operations, fn {length, op} ->
  case op do
    "M" -> IO.puts("#{length} bases match/mismatch")
    "I" -> IO.puts("#{length} bases insertion")
    "D" -> IO.puts("#{length} bases deletion")
    _ -> IO.puts("#{length} bases #{op}")
  end
end)
```

### Writing SAM Files

```elixir
# Write a SAM file
SamParser.write_sam(sam_file, "path/to/output.sam")

# Modify and then write
modified = %{sam_file | alignments: Enum.take(sam_file.alignments, 100)}
SamParser.write_sam(modified, "path/to/subset.sam")
```

## Working with Optional Tags

```elixir
alignment = Enum.at(sam_file.alignments, 0)

# Access optional tags
if Map.has_key?(alignment.tags, "NM") do
  {_type, edit_distance} = alignment.tags["NM"]
  IO.puts("Edit distance: #{edit_distance}")
end

if Map.has_key?(alignment.tags, "MD") do
  {_type, md_string} = alignment.tags["MD"]
  IO.puts("MD string: #{md_string}")
end
```

## Example Workflow

```elixir
# Complete example workflow
defmodule MyGenomicAnalysis do
  def analyze_coverage(sam_path, target_chromosome, region_start, region_end) do
    # Parse the SAM file
    sam_file = SamParser.parse_file(sam_path)
    
    # Filter to relevant region
    region_data = sam_file
      |> SamParser.filter_by_reference(target_chromosome)
      |> SamParser.filter_by_position(region_start, region_end)
      
    # Count reads in region
    read_count = length(region_data.alignments)
    
    # Calculate basic coverage
    region_size = region_end - region_start + 1
    coverage = read_count / region_size
    
    IO.puts("Region #{target_chromosome}:#{region_start}-#{region_end}")
    IO.puts("Read count: #{read_count}")
    IO.puts("Average coverage: #{coverage}")
    
    # Return the filtered data for further processing
    region_data
  end
end
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
