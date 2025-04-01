defmodule SamParser.BamParser do
  @moduledoc """
  Module for parsing BAM files according to the SAM format specification v1.6.

  BAM is the binary version of SAM with BGZF compression. This module handles
  the decompression and binary parsing of BAM files.
  """

  import Bitwise
  alias SamParser.{Alignment, Header, SamFile}
  require Logger

  @typedoc "Reference sequence entry in format {name, length}"
  @type ref_seq :: {String.t(), non_neg_integer()}

  @bam_magic "BAM\1"
  @bgzf_eof_marker <<31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 27, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0>>

  @doc """
  Parses a BAM file from the given path and returns a SamFile struct.

  ## Examples

      iex> SamParser.BamParser.parse_bam("test.bam")
      %SamParser.SamFile{header: %SamParser.Header{...}, alignments: [...]}

  """
  @spec parse_bam(String.t()) :: SamFile.t() | no_return()
  def parse_bam(path) do
    with true <- File.exists?(path) || raise("BAM file not found: #{path}"),
         {:ok, compressed_data} <- File.read(path),
         decompressed_data when byte_size(decompressed_data) > 0 <- decompress_bgzf(compressed_data),
         {:ok, sam_file} <- do_parse_bam_content(decompressed_data) do
      sam_file
    else
      false -> raise "BAM file not found: #{path}"
      {:error, reason} -> raise "Failed to parse BAM file: #{reason}"
      <<>> -> raise "BAM decompression failed: empty result"
      error -> raise "Failed to parse BAM file: #{inspect(error)}"
    end
  end

  # Private function to handle parsing of the BAM content with proper magic byte checking
  defp do_parse_bam_content(data) do
    case data do
      <<@bam_magic, _::binary>> ->
        parse_bam_content(data)

      <<"BAM", 1, _::binary>> ->
        parse_bam_content(data)

      _ ->
        header_bytes = binary_part(data, 0, min(4, byte_size(data)))
        {:error, "Invalid BAM format: magic bytes mismatch. Got: #{inspect(header_bytes, binaries: :as_binaries)}"}
    end
  end

  @doc """
  Decompress BGZF-compressed BAM file contents.

  BAM files use Block GZip Format (BGZF) compression, which consists of
  multiple concatenated gzip blocks.

  ## Parameters
    * `compressed_data` - The binary compressed data from a BAM file

  ## Returns
    Binary data containing the decompressed BAM file content
  """
  @spec decompress_bgzf(binary()) :: binary()
  def decompress_bgzf(compressed_data) do
    Logger.debug("Input compressed data size: #{byte_size(compressed_data)} bytes")

    # First try standard gzip decompression for simpler files
    try_standard_decompression(compressed_data)
  end

  defp try_standard_decompression(compressed_data) do
    try do
      case :zlib.gunzip(compressed_data) do
        decompressed when is_binary(decompressed) and byte_size(decompressed) > 0 ->
          Logger.debug("Successfully decompressed with standard gzip")
          decompressed
        _ ->
          Logger.debug("Standard gzip decompression produced empty result, trying BGZF blocks")
          decompress_bgzf_blocks(compressed_data)
      end
    rescue
      error ->
        Logger.debug("Standard gzip decompression failed: #{inspect(error)}, trying BGZF blocks")
        decompress_bgzf_blocks(compressed_data)
    end
  end

  defp decompress_bgzf_blocks(compressed_data) do
    # Extract and decompress all BGZF blocks
    {blocks, _} = extract_bgzf_blocks(compressed_data, 0, [])

    block_count = length(blocks)
    Logger.debug("Found #{block_count} BGZF blocks")

    # Debug first block if any
    if block_count > 0 do
      [first_block | _] = blocks
      block_size = byte_size(first_block)
      sample_size = min(10, block_size)

      Logger.debug(fn ->
        block_start = binary_part(first_block, 0, sample_size)
        "First block size: #{block_size} bytes, starts with: #{inspect(block_start)}"
      end)
    end

    # Decompress and concatenate all blocks
    blocks
    |> Enum.map(&decompress_block/1)
    |> Enum.join()
  end

  @doc """
  Extracts individual BGZF blocks from the compressed data.

  A BGZF block starts with the gzip magic bytes (1F 8B) and includes
  the block size information.

  ## Parameters
    * `data` - The binary data to search through
    * `pos` - The current position in the binary
    * `blocks` - Accumulated blocks found so far

  ## Returns
    A tuple of {blocks, position} where blocks is a list of binary blocks
  """
  @spec extract_bgzf_blocks(binary(), non_neg_integer(), [binary()]) :: {[binary()], non_neg_integer()}
  def extract_bgzf_blocks(data, pos, blocks) when byte_size(data) <= pos do
    {Enum.reverse(blocks), pos}
  end

  def extract_bgzf_blocks(data, pos, blocks) do
    case find_next_bgzf_block(data, pos) do
      {:ok, block_start, block_size} ->
        block = binary_part(data, block_start, block_size)
        next_pos = block_start + block_size
        extract_bgzf_blocks(data, next_pos, [block | blocks])

      :eof ->
        {Enum.reverse(blocks), pos}

      {:error, :invalid_block} ->
        # Skip this byte and continue searching
        extract_bgzf_blocks(data, pos + 1, blocks)
    end
  end

  # Finds the next BGZF block in the data starting at the given position
  # Returns {:ok, start_position, block_size} or {:error, reason}
  defp find_next_bgzf_block(data, pos) do
    case data do
      # Match BGZF block header (gzip magic + flags)
      <<_::binary-size(pos), 31, 139, 8, 4, rest::binary>> ->
        with true <- byte_size(rest) >= 8,
             <<_::binary-size(2), xlen::little-16, rest2::binary>> <- rest,
             true <- byte_size(rest2) >= xlen,
             {:ok, block_size} <- find_bc_block_size(binary_part(rest2, 0, xlen)),
             total_size = block_size + 1, # As per BAM spec
             true <- pos + total_size <= byte_size(data) do
          {:ok, pos, total_size}
        else
          false -> {:error, :invalid_block}
          _ -> {:error, :invalid_block}
        end

      <<_::binary-size(pos), _::binary>> ->
        {:error, :invalid_block}

      _ ->
        :eof
    end
  end

  # Finds the BC (block content) subfield in the extra field and extracts the block size
  defp find_bc_block_size(extra_field) do
    case :binary.match(extra_field, "BC") do
      {pos, _} ->
        if byte_size(extra_field) >= pos + 4 do
          <<_::binary-size(pos + 2), bsize::little-16, _::binary>> = extra_field
          {:ok, bsize}
        else
          :error
        end

      :nomatch ->
        :error
    end
  end

  # Decompresses a single BGZF block
  defp decompress_block(block) do
    try do
      Logger.debug("Decompressing block of size: #{byte_size(block)} bytes")

      case :zlib.gunzip(block) do
        decompressed when is_binary(decompressed) and byte_size(decompressed) > 0 ->
          Logger.debug("Successfully decompressed to size: #{byte_size(decompressed)} bytes")
          decompressed
        _ ->
          Logger.debug("Block decompressed to empty binary")
          <<>>
      end
    rescue
      error ->
        Logger.debug("Error decompressing block: #{inspect(error)}")
        # Check if it's the EOF marker
        if block == @bgzf_eof_marker do
          Logger.debug("Found EOF marker")
        end
        <<>>
    end
  end

  @doc """
  Parses BAM binary content into a SamFile struct.

  ## Parameters
    * `data` - The binary content of a BAM file after decompression

  ## Returns
    * `{:ok, sam_file}` - Successfully parsed BAM file
    * `{:error, reason}` - Failed to parse BAM file
  """
  @spec parse_bam_content(binary()) :: {:ok, SamFile.t()} | {:error, atom() | String.t()}
  def parse_bam_content(<<@bam_magic, l_text::little-32, header_text::binary-size(l_text), rest::binary>>) do
    # Parse reference sequences
    case rest do
      <<n_ref::little-32, rest_after_n_ref::binary>> ->
        with {ref_seqs, alignment_data} <- parse_reference_sequences(rest_after_n_ref, n_ref, []),
             header <- parse_bam_header(header_text, ref_seqs),
             alignments <- parse_bam_alignments(alignment_data, header, []) do
          {:ok, %SamFile{header: header, alignments: alignments}}
        end

      _ ->
        {:error, "Invalid BAM format: missing reference sequence count"}
    end
  end

  # Handle the case when binary starts with "BAM\1" but doesn't match the constant
  def parse_bam_content(<<"BAM", 1, l_text::little-32, header_text::binary-size(l_text), rest::binary>>) do
    # Parse reference sequences
    case rest do
      <<n_ref::little-32, rest_after_n_ref::binary>> ->
        with {ref_seqs, alignment_data} <- parse_reference_sequences(rest_after_n_ref, n_ref, []),
             header <- parse_bam_header(header_text, ref_seqs),
             alignments <- parse_bam_alignments(alignment_data, header, []) do
          {:ok, %SamFile{header: header, alignments: alignments}}
        end

      _ ->
        {:error, "Invalid BAM format: missing reference sequence count"}
    end
  end

  def parse_bam_content(_) do
    {:error, "Invalid BAM format: magic bytes not found"}
  end

  @doc """
  Parses the BAM header text into a Header struct.

  ## Parameters
    * `header_text` - The text part of the BAM header
    * `ref_seqs` - List of reference sequences from the BAM binary section

  ## Returns
    A SamParser.Header struct
  """
  @spec parse_bam_header(binary(), [{String.t(), non_neg_integer()}]) :: Header.t()
  def parse_bam_header(header_text, ref_seqs) do
    # Split header text into lines and parse as regular SAM header
    lines = String.split(header_text, "\n", trim: true)
    header = SamParser.parse_header(lines)

    # Add reference sequences to header if not present in header text
    sq_names = MapSet.new(Enum.map(header.sq, & &1["SN"]))

    # Convert ref_seqs to SQ entries and add them if not already present
    sq_entries =
      Enum.reduce(ref_seqs, header.sq, fn {name, length}, acc ->
        if MapSet.member?(sq_names, name) do
          acc
        else
          [%{"SN" => name, "LN" => Integer.to_string(length)} | acc]
        end
      end)

    %{header | sq: sq_entries}
  end

  @doc """
  Parses reference sequences from BAM file.

  ## Parameters
    * `data` - Binary data containing reference sequences
    * `n_ref` - Number of reference sequences to parse
    * `acc` - Accumulator for parsed reference sequences

  ## Returns
    A tuple of {reference_sequences, remaining_data}
  """
  @spec parse_reference_sequences(binary(), non_neg_integer(), [ref_seq()]) ::
        {[ref_seq()], binary()}
  def parse_reference_sequences(data, 0, acc) do
    {Enum.reverse(acc), data}
  end

  def parse_reference_sequences(data, n_ref, acc) do
    case data do
      <<l_name::little-32, name_with_null::binary-size(l_name), l_ref::little-32, rest::binary>> ->
        # Remove null terminator
        name = binary_part(name_with_null, 0, l_name - 1)
        parse_reference_sequences(rest, n_ref - 1, [{name, l_ref} | acc])

      _ ->
        Logger.error("Invalid reference sequence data")
        {Enum.reverse(acc), data}
    end
  end

  @doc """
  Parses BAM alignments from binary data.

  ## Parameters
    * `data` - Binary data containing alignment records
    * `header` - The parsed BAM header
    * `acc` - Accumulator for parsed alignments

  ## Returns
    List of parsed alignments
  """
  @spec parse_bam_alignments(binary(), Header.t(), [Alignment.t()]) :: [Alignment.t()]
  def parse_bam_alignments(<<>>, _header, acc) do
    Enum.reverse(acc)
  end

  def parse_bam_alignments(data, header, acc) do
    case parse_single_alignment(data, header) do
      {alignment, rest} ->
        parse_bam_alignments(rest, header, [alignment | acc])
      :eof ->
        Enum.reverse(acc)
    end
  end

  @doc """
  Parses a single BAM alignment record.

  ## Parameters
    * `data` - Binary data containing an alignment record
    * `header` - The parsed BAM header

  ## Returns
    * `{alignment, remaining_data}` - Successfully parsed alignment and remaining data
    * `:eof` - End of file or incomplete data
  """
  @spec parse_single_alignment(binary(), Header.t()) ::
        {Alignment.t(), binary()} | :eof
  def parse_single_alignment(data, header) do
    case data do
      <<block_size::little-32, rest::binary>> when byte_size(rest) >= block_size ->
        <<aln_block::binary-size(block_size), remaining::binary>> = rest
        alignment = parse_alignment_block(aln_block, header)
        {alignment, remaining}
      _ ->
        :eof
    end
  end

  @doc """
  Parses a BAM alignment block into an Alignment struct.

  ## Parameters
    * `block` - Binary data for a single alignment record
    * `header` - The parsed BAM header with reference sequences

  ## Returns
    A SamParser.Alignment struct
  """
  @spec parse_alignment_block(binary(), Header.t()) :: Alignment.t()
  def parse_alignment_block(block, header) do
    # Extract mandatory fields with pattern matching
    <<
      ref_id::little-signed-32,
      pos::little-signed-32,
      l_read_name::unsigned-8,
      mapq::unsigned-8,
      _bin::little-16,  # Unused bin field
      n_cigar_op::little-16,
      flag::little-16,
      l_seq::little-32,
      next_ref_id::little-signed-32,
      next_pos::little-signed-32,
      tlen::little-signed-32,
      rest::binary
    >> = block

    # Get reference name from reference ID
    rname = get_reference_name(ref_id, header)

    # Get next reference name or special markers
    rnext = get_next_reference_name(next_ref_id, ref_id, header)

    # Parse read name (null-terminated)
    <<read_name::binary-size(l_read_name - 1), 0, rest_after_name::binary>> = rest

    # Extract and process CIGAR, sequence, and quality data
    {cigar_string, rest_after_cigar} = parse_cigar_data(rest_after_name, n_cigar_op)
    {seq, qual, rest_after_seq} = parse_seq_and_qual(rest_after_cigar, l_seq)

    # Parse optional fields
    tags = parse_auxiliary_fields(rest_after_seq, %{})

    # Create alignment struct with 0-based to 1-based position conversion
    %Alignment{
      qname: read_name,
      flag: flag,
      rname: rname,
      pos: pos + 1,  # Convert from 0-based to 1-based
      mapq: mapq,
      cigar: cigar_string,
      rnext: rnext,
      pnext: next_pos + 1,  # Convert from 0-based to 1-based
      tlen: tlen,
      seq: seq,
      qual: qual,
      tags: tags
    }
  end

  # Gets the reference name from the reference ID
  @spec get_reference_name(integer(), Header.t()) :: String.t()
  defp get_reference_name(ref_id, header) do
    cond do
      ref_id == -1 -> "*"  # Unmapped read
      ref_id >= 0 && ref_id < length(header.sq) ->
        Enum.at(header.sq, ref_id)["SN"]
      true -> "*"  # Out of range or invalid
    end
  end

  # Gets the next reference name from the next reference ID
  @spec get_next_reference_name(integer(), integer(), Header.t()) :: String.t()
  defp get_next_reference_name(next_ref_id, ref_id, header) do
    cond do
      next_ref_id == -1 -> "*"  # Mate unmapped
      next_ref_id == ref_id -> "="  # Same reference as query
      next_ref_id >= 0 && next_ref_id < length(header.sq) ->
        Enum.at(header.sq, next_ref_id)["SN"]
      true -> "*"  # Out of range or invalid
    end
  end

  @doc """
  Parses CIGAR operations from binary data.

  ## Parameters
    * `data` - Binary data containing CIGAR operations
    * `n_cigar_op` - Number of CIGAR operations to parse

  ## Returns
    A tuple of {cigar_string, remaining_data}
  """
  @spec parse_cigar_data(binary(), non_neg_integer()) :: {String.t(), binary()}
  def parse_cigar_data(data, n_cigar_op) do
    {cigar_ops, rest} = parse_cigar_ops(data, n_cigar_op, [])
    cigar_string = cigar_ops_to_string(cigar_ops)
    {cigar_string, rest}
  end

  # Recursive function to parse CIGAR operations
  @spec parse_cigar_ops(binary(), non_neg_integer(), [{non_neg_integer(), String.t()}]) ::
        {[{non_neg_integer(), String.t()}], binary()}
  defp parse_cigar_ops(data, 0, acc) do
    {Enum.reverse(acc), data}
  end

  defp parse_cigar_ops(data, n, acc) do
    <<cigar::little-32, rest::binary>> = data
    op_len = cigar >>> 4
    op_code = cigar &&& 0xF

    op_char = cigar_op_to_char(op_code)
    parse_cigar_ops(rest, n - 1, [{op_len, op_char} | acc])
  end

  # Map CIGAR operation codes to characters
  @spec cigar_op_to_char(integer()) :: String.t()
  defp cigar_op_to_char(op_code) do
    case op_code do
      0 -> "M" # alignment match
      1 -> "I" # insertion
      2 -> "D" # deletion
      3 -> "N" # skipped region (e.g., intron)
      4 -> "S" # soft clipping
      5 -> "H" # hard clipping
      6 -> "P" # padding
      7 -> "=" # sequence match
      8 -> "X" # sequence mismatch
      _ -> "?" # unknown (shouldn't occur)
    end
  end

  # Convert CIGAR operations to string representation
  @spec cigar_ops_to_string([{non_neg_integer(), String.t()}]) :: String.t()
  defp cigar_ops_to_string([]) do
    "*"  # Empty CIGAR
  end

  defp cigar_ops_to_string(ops) do
    Enum.map_join(ops, "", fn {len, op} -> "#{len}#{op}" end)
  end

  @doc """
  Parses sequence and quality scores from binary data.

  ## Parameters
    * `data` - Binary data containing sequence and quality
    * `l_seq` - Length of the sequence

  ## Returns
    Tuple of {sequence, quality, remaining_data}
  """
  @spec parse_seq_and_qual(binary(), non_neg_integer()) :: {String.t(), String.t(), binary()}
  def parse_seq_and_qual(data, l_seq) do
    case l_seq do
      0 ->
        # No sequence data
        {"*", "*", data}
      _ ->
        # Calculate bytes needed for sequence (4-bit packed)
        seq_bytes = div(l_seq + 1, 2)

        <<seq_data::binary-size(seq_bytes), qual_data::binary-size(l_seq), rest::binary>> = data

        seq = decode_sequence(seq_data, l_seq)
        qual = decode_quality(qual_data)

        {seq, qual, rest}
    end
  end

  # Decode the 4-bit encoded sequence
  @spec decode_sequence(binary(), non_neg_integer()) :: String.t()
  defp decode_sequence(_data, 0) do
    "*"
  end

  defp decode_sequence(data, length) do
    decode_sequence_bytes(data, length, [])
    |> Enum.reverse()
    |> Enum.join()
  end

  # Recursive function to decode sequence bytes
  @spec decode_sequence_bytes(binary(), non_neg_integer(), [String.t()]) :: [String.t()]
  defp decode_sequence_bytes(<<>>, 0, acc) do
    acc
  end

  defp decode_sequence_bytes(<<byte, rest::binary>>, length, acc) when length >= 2 do
    # Extract two bases from each byte
    base1 = decode_base(byte >>> 4)
    base2 = decode_base(byte &&& 0xF)
    decode_sequence_bytes(rest, length - 2, [base2, base1 | acc])
  end

  defp decode_sequence_bytes(<<byte, _rest::binary>>, 1, acc) do
    # Only one base left
    base = decode_base(byte >>> 4)
    [base | acc]
  end

  # Map 4-bit encoded values to nucleotides
  @spec decode_base(integer()) :: String.t()
  defp decode_base(value) do
    case value do
      0 -> "="  # Equal to reference (not typically used in SEQ)
      1 -> "A"  # Adenine
      2 -> "C"  # Cytosine
      3 -> "M"  # aMino (A or C)
      4 -> "G"  # Guanine
      5 -> "R"  # puRine (A or G)
      6 -> "S"  # Strong (C or G)
      7 -> "V"  # not T (A, C or G)
      8 -> "T"  # Thymine
      9 -> "W"  # Weak (A or T)
      10 -> "Y" # pYrimidine (C or T)
      11 -> "H" # not G (A, C or T)
      12 -> "K" # Keto (G or T)
      13 -> "D" # not C (A, G or T)
      14 -> "B" # not A (C, G or T)
      15 -> "N" # aNy (A, C, G or T)
      _ -> "?"  # Invalid value
    end
  end

  # Decode quality scores
  @spec decode_quality(binary()) :: String.t()
  defp decode_quality(<<0xFF, _::binary>>) do
    "*"  # Quality is not available (0xFF bytes)
  end

  defp decode_quality(data) do
    # Convert to ASCII Phred+33 encoding
    data
    |> :binary.bin_to_list()
    |> Enum.map(&(&1 + 33))
    |> List.to_string()
  end

  @doc """
  Parses auxiliary fields (tags) from binary data.

  ## Parameters
    * `data` - Binary data containing auxiliary fields
    * `acc` - Accumulator for parsed tags

  ## Returns
    Map of tag keys to {type, value} pairs
  """
  @spec parse_auxiliary_fields(binary(), map()) :: map()
  def parse_auxiliary_fields(<<>>, acc) do
    acc
  end

  def parse_auxiliary_fields(data, acc) when byte_size(data) < 3 do
    # Not enough data for a valid tag
    acc
  end

  def parse_auxiliary_fields(data, acc) do
    <<tag1, tag2, value_type, rest::binary>> = data
    tag = <<tag1, tag2>>

    case parse_tag_value(value_type, rest) do
      {value, remaining} ->
        # Add tag to accumulator and continue parsing
        updated_acc = Map.put(acc, tag, {<<value_type>>, value})
        parse_auxiliary_fields(remaining, updated_acc)
      :error ->
        # Error parsing tag, return accumulated tags so far
        Logger.warning("Error parsing auxiliary field with tag #{tag} and type #{<<value_type>>}")
        acc
    end
  end

  # Parse tag values based on their type
  @spec parse_tag_value(byte(), binary()) :: {term(), binary()} | :error

  # Character
  defp parse_tag_value(?A, <<value, rest::binary>>) do
    {<<value>>, rest}
  end

  # Int8
  defp parse_tag_value(?c, <<value::signed-8, rest::binary>>) do
    {value, rest}
  end

  # UInt8
  defp parse_tag_value(?C, <<value::unsigned-8, rest::binary>>) do
    {value, rest}
  end

  # Int16
  defp parse_tag_value(?s, <<value::little-signed-16, rest::binary>>) do
    {value, rest}
  end

  # UInt16
  defp parse_tag_value(?S, <<value::little-unsigned-16, rest::binary>>) do
    {value, rest}
  end

  # Int32
  defp parse_tag_value(?i, <<value::little-signed-32, rest::binary>>) do
    {value, rest}
  end

  # UInt32
  defp parse_tag_value(?I, <<value::little-unsigned-32, rest::binary>>) do
    {value, rest}
  end

  # Float
  defp parse_tag_value(?f, <<value::little-float-32, rest::binary>>) do
    {value, rest}
  end

  # String
  defp parse_tag_value(?Z, data) do
    case :binary.split(data, <<0>>) do
      [string, rest] -> {string, rest}
      _ -> :error
    end
  end

  # Hex string
  defp parse_tag_value(?H, data) do
    case :binary.split(data, <<0>>) do
      [hex_string, rest] -> {hex_string, rest}
      _ -> :error
    end
  end

  # Array
  defp parse_tag_value(?B, <<sub_type, count::little-32, rest::binary>>) do
    parse_array_values(sub_type, count, rest, [])
  end

  # Invalid type
  defp parse_tag_value(_, _) do
    :error
  end

  # Parse array values based on their type
  @spec parse_array_values(byte(), non_neg_integer(), binary(), list()) ::
        {list(), binary()} | :error
  defp parse_array_values(_sub_type, 0, rest, acc) do
    {Enum.reverse(acc), rest}
  end

  # Parse array of Int8 values
  defp parse_array_values(?c, count, <<value::signed-8, rest::binary>>, acc) do
    parse_array_values(?c, count - 1, rest, [value | acc])
  end

  # Parse array of UInt8 values
  defp parse_array_values(?C, count, <<value::unsigned-8, rest::binary>>, acc) do
    parse_array_values(?C, count - 1, rest, [value | acc])
  end

  # Parse array of Int16 values
  defp parse_array_values(?s, count, <<value::little-signed-16, rest::binary>>, acc) do
    parse_array_values(?s, count - 1, rest, [value | acc])
  end

  # Parse array of UInt16 values
  defp parse_array_values(?S, count, <<value::little-unsigned-16, rest::binary>>, acc) do
    parse_array_values(?S, count - 1, rest, [value | acc])
  end

  # Parse array of Int32 values
  defp parse_array_values(?i, count, <<value::little-signed-32, rest::binary>>, acc) do
    parse_array_values(?i, count - 1, rest, [value | acc])
  end

  # Parse array of UInt32 values
  defp parse_array_values(?I, count, <<value::little-unsigned-32, rest::binary>>, acc) do
    parse_array_values(?I, count - 1, rest, [value | acc])
  end

  # Parse array of Float values
  defp parse_array_values(?f, count, <<value::little-float-32, rest::binary>>, acc) do
    parse_array_values(?f, count - 1, rest, [value | acc])
  end

  # Error case: not enough data or invalid sub-type
  defp parse_array_values(_, _, _, _) do
    :error
  end
end
