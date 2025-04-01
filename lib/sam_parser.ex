defmodule SamParser do
  @moduledoc """
  Module for parsing SAM/BAM files according to the SAM format specification v1.6.

  This module allows loading, parsing and accessing data from Sequence Alignment/Map format files,
  supporting both text-based SAM and binary BAM formats.
  """

  import Bitwise
  alias SamParser.BamParser

  # Define structs for storing SAM file components
  defmodule Header do
    @moduledoc "Structure for storing SAM header information"
    defstruct hd: nil, sq: [], rg: [], pg: [], co: []
  end

  defmodule Alignment do
    @moduledoc "Structure for storing SAM alignment records"
    defstruct [
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
    ]
  end

  defmodule SamFile do
    @moduledoc "Structure for storing a complete SAM file with header and alignments"
    defstruct header: %SamParser.Header{}, alignments: []
  end

  @doc """
  Parses a SAM or BAM file based on the file extension.
  """
  def parse_file(path) do
    cond do
      String.ends_with?(path, ".sam") -> parse_sam(path)
      String.ends_with?(path, ".bam") -> parse_bam(path)
      true -> raise "Unsupported file format. Expected .sam or .bam file"
    end
  end

  @doc """
  Parses a SAM file from the given path and returns a SamFile struct.
  """
  def parse_sam(path) do
    file_content = File.read!(path)
    lines = String.split(file_content, ~r/\r?\n/, trim: true)

    {header_lines, alignment_lines} = Enum.split_with(lines, &String.starts_with?(&1, "@"))

    %SamFile{
      header: parse_header(header_lines),
      alignments: Enum.map(alignment_lines, &parse_alignment/1)
    }
  end

  @doc """
  Parses a BAM file from the given path and returns a SamFile struct.
  """
  def parse_bam(path) do
    BamParser.parse_bam(path)
  end

  @doc """
  Parses the header section from a list of header lines.
  """
  def parse_header(header_lines) do
    Enum.reduce(header_lines, %Header{}, fn line, header ->
      cond do
        String.starts_with?(line, "@HD") -> %{header | hd: parse_header_line(line)}
        String.starts_with?(line, "@SQ") -> %{header | sq: [parse_header_line(line) | header.sq]}
        String.starts_with?(line, "@RG") -> %{header | rg: [parse_header_line(line) | header.rg]}
        String.starts_with?(line, "@PG") -> %{header | pg: [parse_header_line(line) | header.pg]}
        String.starts_with?(line, "@CO") ->
          comment = String.trim_leading(line, "@CO\t")
          %{header | co: [comment | header.co]}
        true -> header
      end
    end)
    |> then(fn header ->
      # Reverse the lists since we've been prepending
      %{header |
        sq: Enum.reverse(header.sq),
        rg: Enum.reverse(header.rg),
        pg: Enum.reverse(header.pg),
        co: Enum.reverse(header.co)
      }
    end)
  end

  @doc """
  Parses a header line into a map of tag:value pairs.
  """
  def parse_header_line(line) do
    [type | fields] = String.split(line, "\t")
    type = String.trim_leading(type, "@")

    Enum.reduce(fields, %{type: type}, fn field, acc ->
      case String.split(field, ":", parts: 2) do
        [tag, value] -> Map.put(acc, tag, value)
        _ -> acc
      end
    end)
  end

  @doc """
  Parses an alignment line from SAM into an Alignment struct.
  """
  def parse_alignment(line) do
    fields = String.split(line, "\t")

    if length(fields) < 11 do
      raise "Invalid SAM alignment: fewer than 11 mandatory fields"
    end

    [qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual | optional] = fields

    alignment = %Alignment{
      qname: qname,
      flag: String.to_integer(flag),
      rname: rname,
      pos: String.to_integer(pos),
      mapq: String.to_integer(mapq),
      cigar: cigar,
      rnext: rnext,
      pnext: String.to_integer(pnext),
      tlen: String.to_integer(tlen),
      seq: seq,
      qual: qual
    }

    parse_optional_fields(alignment, optional)
  end

  @doc """
  Parses optional fields in TAG:TYPE:VALUE format and adds them to the alignment.
  """
  def parse_optional_fields(alignment, []), do: alignment
  def parse_optional_fields(alignment, optional_fields) do
    tags = Enum.reduce(optional_fields, %{}, fn field, acc ->
      case String.split(field, ":", parts: 3) do
        [tag, type, value] -> Map.put(acc, tag, {type, parse_tag_value(type, value)})
        _ -> acc
      end
    end)

    %{alignment | tags: tags}
  end

  @doc """
  Parse the value of an optional field based on its type.
  """
  def parse_tag_value("A", value), do: value
  def parse_tag_value("i", value), do: String.to_integer(value)
  def parse_tag_value("f", value), do: String.to_float(value)
  def parse_tag_value("Z", value), do: value
  def parse_tag_value("H", value), do: value  # Hex string, could parse to bytes if needed
  def parse_tag_value("B", value) do
    [array_type | elements] = String.split(value, ",")

    case array_type do
      "c" -> Enum.map(elements, &String.to_integer/1) # int8 array
      "C" -> Enum.map(elements, &String.to_integer/1) # uint8 array
      "s" -> Enum.map(elements, &String.to_integer/1) # int16 array
      "S" -> Enum.map(elements, &String.to_integer/1) # uint16 array
      "i" -> Enum.map(elements, &String.to_integer/1) # int32 array
      "I" -> Enum.map(elements, &String.to_integer/1) # uint32 array
      "f" -> Enum.map(elements, &String.to_float/1)   # float array
      _ -> elements  # Default as strings if type is unknown
    end
  end

  @doc """
  Parse CIGAR string into a list of operations.
  Returns a list of {op_length, op_type} tuples.
  """
  def parse_cigar("*"), do: []
  def parse_cigar(cigar) do
    Regex.scan(~r/(\d+)([MIDNSHP=X])/, cigar)
    |> Enum.map(fn [_, length, op] ->
      {String.to_integer(length), op}
    end)
  end

  @doc """
  Checks if a read is properly mapped based on its FLAG field.
  """
  def is_mapped?(alignment) do
    (alignment.flag &&& 0x4) == 0
  end

  @doc """
  Checks if a read is paired based on its FLAG field.
  """
  def is_paired?(alignment) do
    (alignment.flag &&& 0x1) != 0
  end

  @doc """
  Checks if a read is properly paired based on its FLAG field.
  """
  def is_properly_paired?(alignment) do
    (alignment.flag &&& 0x2) != 0
  end

  @doc """
  Checks if a read is reverse complemented based on its FLAG field.
  """
  def is_reverse?(alignment) do
    (alignment.flag &&& 0x10) != 0
  end

  @doc """
  Checks if a read is a secondary alignment based on its FLAG field.
  """
  def is_secondary?(alignment) do
    (alignment.flag &&& 0x100) != 0
  end

  @doc """
  Checks if a read is a supplementary alignment based on its FLAG field.
  """
  def is_supplementary?(alignment) do
    (alignment.flag &&& 0x800) != 0
  end

  @doc """
  Returns a list of reference sequences from the header.
  """
  def reference_sequences(sam_file) do
    Enum.map(sam_file.header.sq, fn sq -> sq["SN"] end)
  end

  @doc """
  Filters alignments by reference name.
  """
  def filter_by_reference(sam_file, reference) do
    filtered = Enum.filter(sam_file.alignments, fn aln ->
      aln.rname == reference
    end)

    %{sam_file | alignments: filtered}
  end

  @doc """
  Filters alignments by position range.
  """
  def filter_by_position(sam_file, start_pos, end_pos) do
    filtered = Enum.filter(sam_file.alignments, fn aln ->
      aln.pos >= start_pos && aln.pos <= end_pos
    end)

    %{sam_file | alignments: filtered}
  end

  @doc """
  Writes a SamFile struct back to a SAM format file.
  """
  def write_sam(sam_file, path) do
    header_lines = format_header(sam_file.header)
    alignment_lines = Enum.map(sam_file.alignments, &format_alignment/1)

    content = Enum.join(header_lines ++ alignment_lines, "\n")
    File.write!(path, content)
  end

  @doc """
  Formats the header section for writing.
  """
  def format_header(header) do
    hd_line = if header.hd, do: [format_header_section("@HD", header.hd)], else: []

    sq_lines = Enum.map(header.sq, &format_header_section("@SQ", &1))
    rg_lines = Enum.map(header.rg, &format_header_section("@RG", &1))
    pg_lines = Enum.map(header.pg, &format_header_section("@PG", &1))
    co_lines = Enum.map(header.co, &"@CO\t#{&1}")

    hd_line ++ sq_lines ++ rg_lines ++ pg_lines ++ co_lines
  end

  @doc """
  Formats a header section for writing.
  """
  def format_header_section(type, fields) do
    field_strings = Enum.map(fields, fn
      {:type, _} -> nil  # Skip the type field as we're using the type parameter
      {tag, value} -> "#{tag}:#{value}"
    end)
    |> Enum.reject(&is_nil/1)

    [type | field_strings] |> Enum.join("\t")
  end

  @doc """
  Formats an alignment record for writing.
  """
  def format_alignment(alignment) do
    mandatory = [
      alignment.qname,
      Integer.to_string(alignment.flag),
      alignment.rname,
      Integer.to_string(alignment.pos),
      Integer.to_string(alignment.mapq),
      alignment.cigar,
      alignment.rnext,
      Integer.to_string(alignment.pnext),
      Integer.to_string(alignment.tlen),
      alignment.seq,
      alignment.qual
    ]

    optional = Enum.map(alignment.tags, fn {tag, {type, value}} ->
      formatted_value = format_tag_value(type, value)
      "#{tag}:#{type}:#{formatted_value}"
    end)

    Enum.join(mandatory ++ optional, "\t")
  end

  @doc """
  Formats a tag value for writing.
  """
  def format_tag_value("B", value) when is_list(value) do
    # Assume first element determines the array type
    {array_type, _} = infer_array_type(value)

    array_values = Enum.map(value, &to_string/1) |> Enum.join(",")
    "#{array_type},#{array_values}"
  end

  def format_tag_value(_type, value), do: to_string(value)

  @doc """
  Infers the array type for a B tag.
  """
  def infer_array_type([first | _]) when is_integer(first) and first >= -128 and first <= 127, do: {"c", "int8"}
  def infer_array_type([first | _]) when is_integer(first) and first >= 0 and first <= 255, do: {"C", "uint8"}
  def infer_array_type([first | _]) when is_integer(first) and first >= -32768 and first <= 32767, do: {"s", "int16"}
  def infer_array_type([first | _]) when is_integer(first) and first >= 0 and first <= 65535, do: {"S", "uint16"}
  def infer_array_type([first | _]) when is_integer(first), do: {"i", "int32"}
  def infer_array_type([first | _]) when is_float(first), do: {"f", "float"}
  def infer_array_type(_), do: {"i", "int32"}  # Default to int32 if empty or unknown
end
