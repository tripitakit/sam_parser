# SAM Parser

An Elixir library for parsing and manipulating SAM (Sequence Alignment/Map) and BAM (Binary Alignment/Map) format files according to the [SAM format specification v1.6](https://samtools.github.io/hts-specs/SAMv1.pdf).


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

## Usage

### Basic File Handling

Load and parse SAM/BAM files with a simple API:

```elixir
# Parse a file based on its extension (.sam or .bam)
sam_file = SamParser.parse_file("path/to/alignment.sam")

# Explicitly parse SAM or BAM
sam_data = SamParser.parse_sam("path/to/alignment.sam")
bam_data = SamParser.parse_bam("path/to/alignment.bam")

# Write modified data back to a SAM file
SamParser.write_sam(sam_file, "path/to/output.sam")
```

### Working with Headers

SAM headers contain important metadata about the reference sequences, read groups, and programs:

```elixir
# Access header information
header = sam_file.header

# Get all reference sequences from the header
ref_seqs = SamParser.reference_sequences(sam_file)
IO.inspect(ref_seqs)  # ["chr1", "chr2", ...]

# Access specific header sections
hd_line = header.hd        # @HD line (header)
sq_lines = header.sq       # @SQ lines (reference sequences)
rg_lines = header.rg       # @RG lines (read groups)
pg_lines = header.pg       # @PG lines (programs)
comments = header.co       # @CO lines (comments)

# Example: Print all reference sequence lengths
Enum.each(sq_lines, fn sq ->
  IO.puts("#{sq["SN"]} length: #{sq["LN"]}")
end)
```

### Manipulating Alignments

Access and process individual alignment records:

```elixir
# Get all alignments
alignments = sam_file.alignments

# Access a specific alignment
first_aln = Enum.at(alignments, 0)

# Basic alignment properties
IO.puts("Read name: #{first_aln.qname}")
IO.puts("Mapped to: #{first_aln.rname} at position: #{first_aln.pos}")
IO.puts("Mapping quality: #{first_aln.mapq}")
IO.puts("CIGAR: #{first_aln.cigar}")
IO.puts("Sequence: #{first_aln.seq}")
IO.puts("Quality: #{first_aln.qual}")

# Check flag properties using helper functions
if SamParser.is_mapped?(first_aln) do
  IO.puts("Read is mapped")
else
  IO.puts("Read is unmapped")
end

# Access optional tags
if Map.has_key?(first_aln.tags, "NM") do
  {_type, edit_distance} = first_aln.tags["NM"]
  IO.puts("Edit distance: #{edit_distance}")
end

# Add or modify a tag
updated_aln = %{first_aln | 
  tags: Map.put(first_aln.tags, "AS", {"i", 42})
}
```

### Filtering and Searching

Filter alignments by various criteria:

```elixir
# Get alignments from a specific chromosome
chr1_alignments = SamParser.filter_by_reference(sam_file, "chr1")

# Get alignments within a specific region
region_alignments = SamParser.filter_by_position(sam_file, 1000, 2000)

# Combine filters (alignments on chr1 between positions 1000-2000)
filtered = sam_file
  |> SamParser.filter_by_reference("chr1") 
  |> SamParser.filter_by_position(1000, 2000)

# Custom filtering with Enum functions
paired_reads = Enum.filter(alignments, &SamParser.is_paired?/1)
high_quality = Enum.filter(alignments, fn aln -> aln.mapq > 30 end)
```

### Working with CIGAR Strings

Parse and analyze CIGAR strings to understand alignment details:

```elixir
# Parse a CIGAR string into a list of operations
cigar_ops = SamParser.parse_cigar("3M2I4M1D2M")
# Result: [{3, "M"}, {2, "I"}, {4, "M"}, {1, "D"}, {2, "M"}]

# Get a detailed analysis of a CIGAR string
alignment = Enum.at(sam_file.alignments, 0)
cigar_analysis = SamParser.Alignment.Parser.analyze_cigar(alignment.cigar)

# Example CIGAR statistics
IO.puts("Total aligned reference bases: #{cigar_analysis.aligned_ref_bases}")
IO.puts("Total aligned read bases: #{cigar_analysis.aligned_read_bases}")
IO.puts("Insertions: #{cigar_analysis.insertions}")
IO.puts("Deletions: #{cigar_analysis.deletions}")
IO.puts("Clipped bases: #{cigar_analysis.clipped_bases}")

# Extract reference sequence covered by an alignment
ref_seq = SamParser.Alignment.Parser.extract_reference_sequence(alignment, reference_sequence)

# Create a visual alignment representation
alignment_view = SamParser.Alignment.Parser.create_alignment_view(alignment, reference_sequence)
IO.puts(alignment_view)
# Output example:
# Ref:  ACGTACGTA
#       |||  |||
# Read: ACG--GTA
```

### Flag Handling

SAM flags contain important information about alignment properties:

```elixir
# Get a map of all flag properties for an alignment
flag_map = SamParser.Alignment.Parser.interpret_flags(alignment)

# Check specific flags
if flag_map.paired do
  IO.puts("Read is paired")
end

if flag_map.reversed do
  IO.puts("Read is reverse complemented")
end

# Create a flag from a map of properties
new_flag = SamParser.Alignment.Parser.build_flag(%{
  paired: true,
  proper_pair: true,
  unmapped: false,
  reversed: false,
  first: true,
  last: false
})
# new_flag value: 67 (0x43)
```

### Processing Large Files

For large SAM/BAM files, consider using Stream-based processing:

```elixir
# Process alignments one at a time
File.stream!("large_alignment.sam")
|> Stream.drop_while(&String.starts_with?(&1, "@"))  # Skip header lines
|> Stream.map(&SamParser.parse_alignment/1)
|> Stream.filter(&SamParser.is_mapped?/1)
|> Stream.each(fn alignment ->
    # Process each alignment here
    IO.puts("Processing #{alignment.qname}")
end)
|> Stream.run()
```

## Advanced Usage

### Working with Region Queries

Extract and analyze data from specific genomic regions:

```elixir
defmodule RegionAnalysis do
  def analyze_region(sam_path, chrom, start_pos, end_pos) do
    # Load SAM file
    sam_file = SamParser.parse_file(sam_path)
    
    # Filter alignments in region
    region_data = sam_file
    |> SamParser.filter_by_reference(chrom)
    |> SamParser.filter_by_position(start_pos, end_pos)
    
    # Calculate coverage
    read_count = length(region_data.alignments)
    region_size = end_pos - start_pos + 1
    coverage = read_count / region_size
    
    # Calculate average mapping quality
    avg_mapq = region_data.alignments
    |> Enum.map(& &1.mapq)
    |> Enum.sum()
    |> Kernel./(max(1, read_count))
    
    %{
      region: "#{chrom}:#{start_pos}-#{end_pos}",
      read_count: read_count,
      coverage: coverage,
      avg_mapq: avg_mapq
    }
  end
end

# Usage
result = RegionAnalysis.analyze_region("sample.bam", "chr1", 1000, 2000)
IO.inspect(result)
```

### Creating and Modifying Alignments

Create new alignments or modify existing ones:

```elixir
# Create a new alignment record
new_alignment = %SamParser.Alignment{
  qname: "read_1",
  flag: 0,
  rname: "chr1",
  pos: 100,
  mapq: 60,
  cigar: "100M",
  rnext: "*",
  pnext: 0,
  tlen: 0,
  seq: String.duplicate("A", 100),
  qual: String.duplicate("I", 100),  # Quality score 40 in Phred+33
  tags: %{"NM" => {"i", 0}, "MD" => {"Z", "100"}}
}

# Add the new alignment to a file
updated_file = %{sam_file | alignments: [new_alignment | sam_file.alignments]}

# Modify an alignment property
modified = %{new_alignment | mapq: 50}

# Write the updated file
SamParser.write_sam(updated_file, "updated.sam")
```

### Alignment Statistics

Gather statistics across alignments:

```elixir
defmodule AlignmentStats do
  def calculate(alignments) do
    # Basic counts
    total = length(alignments)
    mapped = Enum.count(alignments, &SamParser.is_mapped?/1)
    unmapped = total - mapped
    paired = Enum.count(alignments, &SamParser.is_paired?/1)
    
    # Quality distribution
    qualities = alignments 
    |> Enum.filter(&SamParser.is_mapped?/1)
    |> Enum.map(& &1.mapq)
    
    avg_quality = Enum.sum(qualities) / max(1, length(qualities))
    
    # Mapping by chromosome
    by_chrom = Enum.reduce(alignments, %{}, fn aln, acc ->
      if aln.rname != "*" do
        Map.update(acc, aln.rname, 1, &(&1 + 1))
      else
        acc
      end
    end)
    
    %{
      total: total,
      mapped: mapped,
      unmapped: unmapped,
      paired: paired,
      mapping_rate: mapped / total,
      avg_quality: avg_quality,
      by_chromosome: by_chrom
    }
  end
end

# Usage
stats = AlignmentStats.calculate(sam_file.alignments)
IO.inspect(stats)
```

## API Reference

### Core Modules

#### `SamParser`

Main module with functions for parsing and manipulating SAM/BAM files.

| Function | Description |
|----------|-------------|
| `parse_file(path)` | Parse a SAM or BAM file based on extension |
| `parse_sam(path)` | Parse a SAM format file |
| `parse_bam(path)` | Parse a BAM format file |
| `write_sam(sam_file, path)` | Write a SamFile struct to a SAM file |
| `parse_header(header_lines)` | Parse SAM header lines |
| `parse_alignment(line)` | Parse a single SAM alignment line |
| `parse_optional_fields(alignment, fields)` | Parse optional fields |
| `parse_cigar(cigar)` | Parse CIGAR string into operations list |
| `reference_sequences(sam_file)` | Get reference sequence names |
| `filter_by_reference(sam_file, reference)` | Filter alignments by reference |
| `filter_by_position(sam_file, start_pos, end_pos)` | Filter alignments by position |
| `is_mapped?(alignment)` | Check if alignment is mapped |
| `is_paired?(alignment)` | Check if alignment is paired |
| `is_properly_paired?(alignment)` | Check if alignment is properly paired |
| `is_reverse?(alignment)` | Check if alignment is reverse strand |
| `is_secondary?(alignment)` | Check if alignment is secondary |
| `is_supplementary?(alignment)` | Check if alignment is supplementary |

#### `SamParser.Alignment.Parser`

Specialized module for parsing and analyzing alignment records.

| Function | Description |
|----------|-------------|
| `interpret_flags(alignment)` | Get a map of flags from an alignment |
| `build_flag(flag_map)` | Build a flag integer from a map |
| `analyze_cigar(cigar)` | Detailed analysis of a CIGAR string |
| `get_end_position(alignment)` | Get the rightmost position of alignment |
| `overlaps_region?(alignment, start_pos, end_pos)` | Check if alignment overlaps a region |
| `extract_reference_sequence(alignment, reference)` | Extract reference sequence |
| `create_alignment_view(alignment, reference)` | Create visual alignment representation |
| `extract_quality_scores(alignment)` | Extract quality scores as integers |

### Data Structures

#### `SamParser.SamFile`

```elixir
%SamParser.SamFile{
  header: %SamParser.Header{},  # Header information
  alignments: []                # List of alignment records
}
```

#### `SamParser.Header`

```elixir
%SamParser.Header{
  hd: %{},  # @HD line fields (header)
  sq: [],   # @SQ lines (reference sequences)
  rg: [],   # @RG lines (read groups)
  pg: [],   # @PG lines (programs)
  co: []    # @CO lines (comments)
}
```

#### `SamParser.Alignment`

```elixir
%SamParser.Alignment{
  qname: "",    # Query name
  flag: 0,      # Bitwise FLAG
  rname: "*",   # Reference name
  pos: 0,       # 1-based position
  mapq: 0,      # Mapping quality
  cigar: "*",   # CIGAR string
  rnext: "*",   # Reference of next segment
  pnext: 0,     # Position of next segment
  tlen: 0,      # Template length
  seq: "*",     # Sequence
  qual: "*",    # Quality string
  tags: %{}     # Optional tags
}
```

## Error Handling

SAM Parser includes comprehensive error handling:

```elixir
try do
  sam_file = SamParser.parse_file("non_existent.sam")
rescue
  e in File.Error -> 
    IO.puts("File error: #{e.reason}")
  e in RuntimeError ->
    IO.puts("Runtime error: #{e.message}")
end
```

## Performance Tips

- For large files, use Stream-based processing instead of loading the entire file
- Filter alignments as early as possible to reduce memory usage
- When writing custom filters, combine them to minimize passes through the data
- Consider using `Flow` or `GenStage` for parallel processing of alignments

## License

This project is licensed under the GNU General Public License v3.0 - see the [COPYING](COPYING) file for details.


## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the project
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request


## Citation

If you use ReadsMap in your research, please cite:

```
De Marta, P. (2025). SAM Parser: An Elixir library for parsing and manipulating SAM and BAM files.
https://github.com/tripitakit/sam_parser
```


## Contact

Author: Patrick De Marta  
Email: patrick.demarta@gmail.com
