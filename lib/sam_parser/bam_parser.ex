defmodule SamParser.BamParser do
  @moduledoc """
  Module for parsing BAM files according to the SAM format specification v1.6.

  BAM is the binary version of SAM with BGZF compression. This module handles
  the decompression and binary parsing of BAM files.
  """

  import Bitwise
  alias SamParser.{Alignment, SamFile}

  @bam_magic "BAM\1"

  @bgzf_eof_marker <<31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 27, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0>>

  @doc """
  Parses a BAM file from the given path and returns a SamFile struct.
  """
  def parse_bam(path) do
    try do
      # Check if file exists
      unless File.exists?(path) do
        raise "BAM file not found: #{path}"
      end

      # Read file as binary
      {:ok, compressed_data} = File.read(path)

      # Decompress the BAM content
      decompressed_data = decompress_bgzf(compressed_data)

      # Check if decompression was successful
      if byte_size(decompressed_data) == 0 do
        raise "BAM decompression failed: empty result"
      end

      # Check for BAM magic bytes
      first_bytes = if byte_size(decompressed_data) >= 4 do
        binary_part(decompressed_data, 0, 4)
      else
        decompressed_data
      end

      # Comparing using inspected values (for debugging)
      if first_bytes != @bam_magic do
        # Try to parse anyway - the magic bytes might be correct but string comparison fails
        case parse_bam_content(decompressed_data) do
          {:ok, sam_file} -> sam_file
          {:error, reason} ->
            raise "Invalid BAM format: magic bytes mismatch. Got: #{inspect(first_bytes, binaries: :as_binaries)}, Expected: #{inspect(@bam_magic, binaries: :as_binaries)}. Error: #{reason}"
        end
      else
        # Parse the decompressed content
        case parse_bam_content(decompressed_data) do
          {:ok, sam_file} -> sam_file
          {:error, reason} -> raise "Failed to parse BAM content: #{reason}"
        end
      end
    rescue
      e ->
        raise "Failed to parse BAM file: #{inspect(e)}"
    end
  end

  @doc """
  Decompress BGZF-compressed BAM file contents.

  BAM files use Block GZip Format (BGZF) compression, which consists of
  multiple concatenated gzip blocks.
  """
  def decompress_bgzf(compressed_data) do
    # Log the size of input data
    IO.puts("Input compressed data size: #{byte_size(compressed_data)}")

    # First try standard gzip decompression
    try do
      case :zlib.gunzip(compressed_data) do
        decompressed when is_binary(decompressed) and byte_size(decompressed) > 0 ->
          IO.puts("Successfully decompressed with standard gzip")
          decompressed
        _ ->
          IO.puts("Standard gzip decompression failed, trying BGZF blocks")
          decompress_bgzf_blocks(compressed_data)
      end
    rescue
      _ ->
        IO.puts("Standard gzip decompression failed with error, trying BGZF blocks")
        decompress_bgzf_blocks(compressed_data)
    end
  end

  defp decompress_bgzf_blocks(compressed_data) do
    # Extract and decompress all BGZF blocks
    {blocks, _} = extract_bgzf_blocks(compressed_data, 0, [])

    # Log number of blocks found
    IO.puts("Found #{length(blocks)} BGZF blocks")

    # Debug first block if any
    if length(blocks) > 0 do
      [first_block | _] = blocks
      IO.puts("First block size: #{byte_size(first_block)}")
      IO.puts("First block starts with: #{inspect(binary_part(first_block, 0, min(10, byte_size(first_block))))}")
    end

    # Decompress and concatenate all blocks
    blocks
    |> Enum.map(&decompress_block/1)
    |> Enum.join()
  end

  @doc """
  Extract individual BGZF blocks from the compressed data.
  A BGZF block starts with the gzip magic bytes (1F 8B) and includes
  the block size information.
  """
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

  defp find_next_bgzf_block(data, pos) do
    case data do
      # Match BGZF block header
      <<_::binary-size(pos), 31, 139, 8, 4, rest::binary>> ->
        if byte_size(rest) >= 8 do
          <<_::binary-size(2), xlen::little-16, rest2::binary>> = rest

          if byte_size(rest2) >= xlen do
            case find_bc_block_size(binary_part(rest2, 0, xlen)) do
              {:ok, block_size} ->
                total_size = block_size + 1  # As per BAM spec
                if pos + total_size <= byte_size(data) do
                  {:ok, pos, total_size}
                else
                  {:error, :invalid_block}
                end

              _ ->
                {:error, :invalid_block}
            end
          else
            {:error, :invalid_block}
          end
        else
          :eof
        end

      <<_::binary-size(pos), _::binary>> ->
        {:error, :invalid_block}

      _ ->
        :eof
    end
  end

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

  defp decompress_block(block) do
    try do
      IO.puts("Attempting to decompress block of size: #{byte_size(block)}")
      case :zlib.gunzip(block) do
        decompressed when is_binary(decompressed) and byte_size(decompressed) > 0 ->
          IO.puts("Successfully decompressed to size: #{byte_size(decompressed)}")
          decompressed
        _ ->
          IO.puts("Block decompressed to empty binary")
          # If decompression results in empty data, return empty binary
          <<>>
      end
    rescue
      e ->
        IO.puts("Error decompressing block: #{inspect(e)}")
        # If it's the EOF marker, return empty binary
        if block == @bgzf_eof_marker do
          IO.puts("Found EOF marker")
          <<>>
        else
          # For other decompression errors, return empty binary
          <<>>
        end
    end
  end

  @doc """
  Parses BAM binary content into a SamFile struct.
  """
  def parse_bam_content(<<@bam_magic, data::binary>>) do
    # Parse header text length
    <<l_text::little-32, rest::binary>> = data

    # Extract header text and parse it
    <<header_text::binary-size(l_text), alignment_data::binary>> = rest

    # Parse reference sequences
    <<n_ref::little-32, rest::binary>> = alignment_data
    {ref_seqs, alignment_data} = parse_reference_sequences(rest, n_ref, [])

    # Parse header
    header = parse_bam_header(header_text, ref_seqs)

    # Parse alignments
    alignments = parse_bam_alignments(alignment_data, header, [])

    {:ok, %SamFile{header: header, alignments: alignments}}
  end

  # Handle the case when binary starts with "BAM\1" but doesn't match the constant
  def parse_bam_content(<<"BAM", 1, data::binary>>) do
    # Parse header text length
    <<l_text::little-32, rest::binary>> = data

    # Extract header text and parse it
    <<header_text::binary-size(l_text), alignment_data::binary>> = rest

    # Parse reference sequences
    <<n_ref::little-32, rest::binary>> = alignment_data
    {ref_seqs, alignment_data} = parse_reference_sequences(rest, n_ref, [])

    # Parse header
    header = parse_bam_header(header_text, ref_seqs)

    # Parse alignments
    alignments = parse_bam_alignments(alignment_data, header, [])

    {:ok, %SamFile{header: header, alignments: alignments}}
  end

  def parse_bam_content(_) do
    {:error, :invalid_bam_format}
  end

  @doc """
  Parses the BAM header text into a Header struct.
  """
  def parse_bam_header(header_text, ref_seqs) do
    # Split header text into lines and parse as regular SAM header
    lines = String.split(header_text, "\n", trim: true)
    header = SamParser.parse_header(lines)

    # Add reference sequences to header if not present in header text
    sq_names = MapSet.new(Enum.map(header.sq, & &1["SN"]))

    sq_entries = Enum.reduce(ref_seqs, header.sq, fn {name, length}, acc ->
      if MapSet.member?(sq_names, name) do
        acc
      else
        [%{"type" => "SQ", "SN" => name, "LN" => Integer.to_string(length)} | acc]
      end
    end)

    %{header | sq: sq_entries}
  end

  @doc """
  Parses reference sequences from BAM file.
  Returns {reference_sequences, remaining_data}
  """
  def parse_reference_sequences(data, 0, acc) do
    {Enum.reverse(acc), data}
  end

  def parse_reference_sequences(data, n_ref, acc) do
    <<l_name::little-32, rest::binary>> = data
    <<name::binary-size(l_name), rest::binary>> = rest

    # Remove null terminator
    name = binary_part(name, 0, l_name - 1)

    <<l_ref::little-32, rest::binary>> = rest

    parse_reference_sequences(rest, n_ref - 1, [{name, l_ref} | acc])
  end

  @doc """
  Parses BAM alignments from binary data.
  """
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
  Returns {alignment, remaining_data}
  """
  def parse_single_alignment(data, header) do
    case data do
      <<block_size::little-32, rest::binary>> ->
        if byte_size(rest) < block_size do
          :eof
        else
          <<aln_block::binary-size(block_size), rest::binary>> = rest
          alignment = parse_alignment_block(aln_block, header)
          {alignment, rest}
        end
      _ ->
        :eof
    end
  end

  @doc """
  Parses a BAM alignment block into an Alignment struct.
  """
  def parse_alignment_block(block, header) do
    <<
      ref_id::little-signed-32,
      pos::little-signed-32,
      l_read_name::unsigned-8,
      mapq::unsigned-8,
      _bin::little-16,  # Unused variable, prefixed with underscore
      n_cigar_op::little-16,
      flag::little-16,
      l_seq::little-32,
      next_ref_id::little-signed-32,
      next_pos::little-signed-32,
      tlen::little-signed-32,
      rest::binary
    >> = block

    # Get reference name
    rname = if ref_id >= 0 and ref_id < length(header.sq) do
      Enum.at(header.sq, ref_id)["SN"]
    else
      "*"
    end

    # Get next reference name
    rnext = cond do
      next_ref_id == -1 -> "*"
      next_ref_id == ref_id -> "="
      next_ref_id >= 0 and next_ref_id < length(header.sq) ->
        Enum.at(header.sq, next_ref_id)["SN"]
      true -> "*"
    end

    # Parse read name
    <<read_name::binary-size(l_read_name - 1), 0, rest::binary>> = rest

    # Parse CIGAR
    {cigar_string, rest} = parse_cigar_data(rest, n_cigar_op)

    # Parse sequence and quality
    {seq, qual, rest} = parse_seq_and_qual(rest, l_seq)

    # Parse auxiliary fields
    tags = parse_auxiliary_fields(rest, %{})

    # Create alignment struct
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

  @doc """
  Parses CIGAR operations from binary data.
  """
  def parse_cigar_data(data, n_cigar_op) do
    {cigar_ops, rest} = parse_cigar_ops(data, n_cigar_op, [])
    cigar_string = cigar_ops_to_string(cigar_ops)
    {cigar_string, rest}
  end

  defp parse_cigar_ops(data, 0, acc) do
    {Enum.reverse(acc), data}
  end

  defp parse_cigar_ops(data, n, acc) do
    <<cigar::little-32, rest::binary>> = data
    op_len = cigar >>> 4
    op_code = cigar &&& 0xF

    op_char = case op_code do
      0 -> "M"
      1 -> "I"
      2 -> "D"
      3 -> "N"
      4 -> "S"
      5 -> "H"
      6 -> "P"
      7 -> "="
      8 -> "X"
      _ -> "?"
    end

    parse_cigar_ops(rest, n - 1, [{op_len, op_char} | acc])
  end

  defp cigar_ops_to_string([]) do
    "*"
  end

  defp cigar_ops_to_string(ops) do
    Enum.map_join(ops, "", fn {len, op} -> "#{len}#{op}" end)
  end

  @doc """
  Parses sequence and quality scores from binary data.
  """
  def parse_seq_and_qual(data, l_seq) do
    # Calculate bytes needed for sequence (4-bit packed)
    seq_bytes = div(l_seq + 1, 2)

    <<seq_data::binary-size(seq_bytes), qual_data::binary-size(l_seq), rest::binary>> = data

    if l_seq == 0 do
      {"*", "*", rest}
    else
      seq = decode_sequence(seq_data, l_seq)
      qual = decode_quality(qual_data)
      {seq, qual, rest}
    end
  end

  defp decode_sequence(_data, 0) do
    "*"
  end

  defp decode_sequence(data, length) do
    # Unpack 4-bit encoded sequence
    decode_sequence_bytes(data, length, [])
    |> Enum.reverse()
    |> Enum.join()
  end

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

  defp decode_base(value) do
    case value do
      0 -> "="
      1 -> "A"
      2 -> "C"
      3 -> "M"
      4 -> "G"
      5 -> "R"
      6 -> "S"
      7 -> "V"
      8 -> "T"
      9 -> "W"
      10 -> "Y"
      11 -> "H"
      12 -> "K"
      13 -> "D"
      14 -> "B"
      15 -> "N"
      _ -> "?"
    end
  end

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
  """
  def parse_auxiliary_fields(<<>>, acc) do
    acc
  end

  def parse_auxiliary_fields(data, acc) do
    if byte_size(data) < 3 do
      acc
    else
      <<tag1, tag2, value_type, rest::binary>> = data
      tag = <<tag1, tag2>>

      {value, rest} = parse_tag_value(value_type, rest)

      parse_auxiliary_fields(rest, Map.put(acc, tag, {<<value_type>>, value}))
    end
  end

  defp parse_tag_value(?A, <<value, rest::binary>>) do
    {<<value>>, rest}
  end

  defp parse_tag_value(?c, <<value::signed-8, rest::binary>>) do
    {value, rest}
  end

  defp parse_tag_value(?C, <<value::unsigned-8, rest::binary>>) do
    {value, rest}
  end

  defp parse_tag_value(?s, <<value::little-signed-16, rest::binary>>) do
    {value, rest}
  end

  defp parse_tag_value(?S, <<value::little-unsigned-16, rest::binary>>) do
    {value, rest}
  end

  defp parse_tag_value(?i, <<value::little-signed-32, rest::binary>>) do
    {value, rest}
  end

  defp parse_tag_value(?I, <<value::little-unsigned-32, rest::binary>>) do
    {value, rest}
  end

  defp parse_tag_value(?f, <<value::little-float-32, rest::binary>>) do
    {value, rest}
  end

  defp parse_tag_value(?Z, data) do
    # Null-terminated string
    [string, rest] = :binary.split(data, <<0>>)
    {string, rest}
  end

  defp parse_tag_value(?H, data) do
    # Hex string
    [hex_string, rest] = :binary.split(data, <<0>>)
    {hex_string, rest}
  end

  defp parse_tag_value(?B, <<sub_type, count::little-32, rest::binary>>) do
    {values, rest} = parse_array_values(sub_type, count, rest, [])
    {values, rest}
  end

  defp parse_array_values(_sub_type, 0, rest, acc) do
    {Enum.reverse(acc), rest}
  end

  defp parse_array_values(?c, count, <<value::signed-8, rest::binary>>, acc) do
    parse_array_values(?c, count - 1, rest, [value | acc])
  end

  defp parse_array_values(?C, count, <<value::unsigned-8, rest::binary>>, acc) do
    parse_array_values(?C, count - 1, rest, [value | acc])
  end

  defp parse_array_values(?s, count, <<value::little-signed-16, rest::binary>>, acc) do
    parse_array_values(?s, count - 1, rest, [value | acc])
  end

  defp parse_array_values(?S, count, <<value::little-unsigned-16, rest::binary>>, acc) do
    parse_array_values(?S, count - 1, rest, [value | acc])
  end

  defp parse_array_values(?i, count, <<value::little-signed-32, rest::binary>>, acc) do
    parse_array_values(?i, count - 1, rest, [value | acc])
  end

  defp parse_array_values(?I, count, <<value::little-unsigned-32, rest::binary>>, acc) do
    parse_array_values(?I, count - 1, rest, [value | acc])
  end

  defp parse_array_values(?f, count, <<value::little-float-32, rest::binary>>, acc) do
    parse_array_values(?f, count - 1, rest, [value | acc])
  end
end
