<img src="logo.png" alt="SAM Parser Logo" width="200"/>

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

## Data Structures

The library represents SAM/BAM data using the following key structures:

### `SamParser.SamFile`

The top-level structure containing both the header and alignment records:

```elixir
%SamParser.SamFile{
  header: %SamParser.Header{},  # Header information
  alignments: []                # List of alignment records
}
```

### `SamParser.Header`

Stores the SAM file header components:

```elixir
%SamParser.Header{
  hd: nil,   # @HD line fields as a map
  sq: [],    # List of @SQ entries (reference sequences)
  rg: [],    # List of @RG entries (read groups)
  pg: [],    # List of @PG entries (programs)
  co: []     # List of @CO entries (comments)
}
```

### `SamParser.Alignment`

Represents a single alignment record with all SAM fields:

```elixir
%SamParser.Alignment{
  # Mandatory fields
  qname: nil,   # Query template NAME
  flag: 0,      # Bitwise FLAG
  rname: "*",   # Reference sequence NAME
  pos: 0,       # 1-based leftmost mapping POSition
  mapq: 0,      # MAPping Quality
  cigar: "*",   # CIGAR string
  rnext: "*",   # Ref. name of the mate/next read
  pnext: 0,     # Position of the mate/next read
  tlen: 0,      # Observed Template LENgth
  seq: "*",     # Segment SEQuence
  qual: "*",    # ASCII of Phred-scaled base QUALity+33
  # Optional field storage
  tags: %{}     # Map to store TAG:TYPE:VALUE optional fields
}
```

Optional tags are stored in the `tags` map as `{tag => {type, value}}` entries.

## API Reference

### File Parsing

| Function | Description |
|----------|-------------|
| `parse_file(path)` | Parses either SAM or BAM file based on file extension |
| `parse_sam(path)` | Parses a SAM format file |
| `parse_bam(path)` | Parses a BAM format file |

### Data Access and Processing

| Function | Description |
|----------|-------------|
| `parse_header(header_lines)` | Parses SAM header lines into a Header struct |
| `parse_alignment(line)` | Parses a single SAM alignment line |
| `parse_optional_fields(alignment, fields)` | Parses optional TAG:TYPE:VALUE fields |
| `parse_tag_value(type, value)` | Parses tag values based on their declared type |
| `parse_cigar(cigar)` | Parses CIGAR string into a list of operations |
| `reference_sequences(sam_file)` | Returns a list of reference sequence names from the header |

### Filtering

| Function | Description |
|----------|-------------|
| `filter_by_reference(sam_file, reference)` | Filters alignments by reference sequence name |
| `filter_by_position(sam_file, start_pos, end_pos)` | Filters alignments by position range |

### Flag Operations

| Function | Description |
|----------|-------------|
| `is_mapped?(alignment)` | Checks if read is mapped (0x4 flag bit) |
| `is_paired?(alignment)` | Checks if read is paired (0x1 flag bit) |
| `is_properly_paired?(alignment)` | Checks if read is properly paired (0x2 flag bit) |
| `is_reverse?(alignment)` | Checks if read is reverse complemented (0x10 flag bit) |
| `is_secondary?(alignment)` | Checks if read is a secondary alignment (0x100 flag bit) |
| `is_supplementary?(alignment)` | Checks if read is a supplementary alignment (0x800 flag bit) |

### File Writing

| Function | Description |
|----------|-------------|
| `write_sam(sam_file, path)` | Writes a SamFile struct to a SAM format file |
| `format_header(header)` | Formats header section for writing |
| `format_header_section(type, fields)` | Formats a header section for writing |
| `format_alignment(alignment)` | Formats an alignment record for writing |
| `format_tag_value(type, value)` | Formats tag values for writing based on their type |

### Utility Functions

| Function | Description |
|----------|-------------|
| `infer_array_type(array)` | Infers the array type for "B" tag values |
| `parse_header_line(line)` | Parses a header line into a map of tag:value pairs |

## Examples of Data Structure Usage

```elixir
# Examine the structure of a parsed SAM file
sam_file = SamParser.parse_file("sample.sam")
IO.inspect(sam_file.header.sq)  # List all reference sequences

# Create a new alignment manually
new_alignment = %SamParser.Alignment{
  qname: "read123",
  flag: 0,
  rname: "chr1",
  pos: 100,
  mapq: 60,
  cigar: "100M",
  seq: "ACGT...",
  qual: "FFFF..."
}

# Add a tag to an alignment
alignment_with_tag = %{new_alignment | 
  tags: Map.put(new_alignment.tags, "NM", {"i", 0})
}

# Add this alignment to the file
updated_sam_file = %{sam_file | 
  alignments: [alignment_with_tag | sam_file.alignments]
}

# Write the updated file
SamParser.write_sam(updated_sam_file, "updated.sam")
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
