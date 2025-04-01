defmodule SamParserTest do
  use ExUnit.Case
  doctest SamParser

  test "parse SAM header" do
    header_lines = [
      "@HD\tVN:1.6\tSO:coordinate",
      "@SQ\tSN:pstS\tLN:1000",
      "@RG\tID:1\tSM:sample1",
      "@PG\tID:minimap2\tPN:minimap2\tVN:2.24-r1122",
      "@CO\tExample SAM file for testing"
    ]

    header = SamParser.parse_header(header_lines)

    assert header.hd["VN"] == "1.6"
    assert header.hd["SO"] == "coordinate"
    assert Enum.at(header.sq, 0)["SN"] == "pstS"
    assert Enum.at(header.sq, 0)["LN"] == "1000"
    assert Enum.at(header.rg, 0)["ID"] == "1"
    assert Enum.at(header.pg, 0)["ID"] == "minimap2"
    assert Enum.at(header.co, 0) == "Example SAM file for testing"
  end

  test "parse SAM alignment" do
    alignment_line = "read1\t0\tpstS\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t!!!!!!!!!!";

    alignment = SamParser.parse_alignment(alignment_line)

    assert alignment.qname == "read1"
    assert alignment.flag == 0
    assert alignment.rname == "pstS"
    assert alignment.pos == 1
    assert alignment.mapq == 60
    assert alignment.cigar == "10M"
    assert alignment.seq == "ACGTACGTAC"
    assert alignment.qual == "!!!!!!!!!!"
  end

  test "parse BAM file" do
    result = SamParser.parse_bam("test/test.bam")

    assert %SamParser.SamFile{} = result
    assert result.header != nil

    # Verify we can decompress and parse the BAM content
    assert length(result.alignments) > 0

    # Check first alignment has expected fields
    first_alignment = Enum.at(result.alignments, 0)
    assert first_alignment.qname != nil
    assert is_integer(first_alignment.flag)
    assert first_alignment.rname != nil
    assert is_integer(first_alignment.pos)
  end

  test "SAM and BAM parsing produce same results" do
    sam_result = SamParser.parse_sam("test/test.sam")

    bam_result = SamParser.parse_bam("test/test.bam")

    # Compare number of alignments
    assert length(sam_result.alignments) == length(bam_result.alignments)

    # Compare first alignment fields
    first_sam_alignment = Enum.at(sam_result.alignments, 0)
    first_bam_alignment = Enum.at(bam_result.alignments, 0)
    assert first_sam_alignment.qname == first_bam_alignment.qname
    assert first_sam_alignment.flag == first_bam_alignment.flag
    assert first_sam_alignment.rname == first_bam_alignment.rname
    assert first_sam_alignment.pos == first_bam_alignment.pos
    assert first_sam_alignment.mapq == first_bam_alignment.mapq
    assert first_sam_alignment.cigar == first_bam_alignment.cigar
    assert first_sam_alignment.seq == first_bam_alignment.seq
    assert first_sam_alignment.qual == first_bam_alignment.qual

    # Compare header fields
    assert sam_result.header.hd == bam_result.header.hd
    assert sam_result.header.sq == bam_result.header.sq
    assert sam_result.header.rg == bam_result.header.rg
    assert sam_result.header.co == bam_result.header.co

  end

  test "interpret_flags parses all flag bits correctly" do
    # Create test alignments with different flag combinations
    alignment_all_flags = %SamParser.Alignment{flag: 0xFFF}  # All bits set
    alignment_no_flags = %SamParser.Alignment{flag: 0x0}     # No bits set

    flags_all = SamParser.Alignment.Parser.interpret_flags(alignment_all_flags)
    flags_none = SamParser.Alignment.Parser.interpret_flags(alignment_no_flags)

    # Test all flag bits
    assert flags_all.paired == true
    assert flags_all.proper_pair == true
    assert flags_all.unmapped == true
    assert flags_all.next_unmapped == true
    assert flags_all.reversed == true
    assert flags_all.next_reversed == true
    assert flags_all.first == true
    assert flags_all.last == true
    assert flags_all.secondary == true
    assert flags_all.filtered == true
    assert flags_all.duplicate == true
    assert flags_all.supplementary == true

    # Test no flags set
    assert flags_none.paired == false
    assert flags_none.proper_pair == false
    assert flags_none.unmapped == false
    assert flags_none.next_unmapped == false
    assert flags_none.reversed == false
    assert flags_none.next_reversed == false
    assert flags_none.first == false
    assert flags_none.last == false
    assert flags_none.secondary == false
    assert flags_none.filtered == false
    assert flags_none.duplicate == false
    assert flags_none.supplementary == false
  end

  test "analyze_cigar parses different CIGAR operations correctly" do
    alignment = %SamParser.Alignment{
      cigar: "10M2I3D4S5H6N=X",
      seq: "ACGTACGTACAA",
      pos: 1
    }

    result = SamParser.Alignment.Parser.analyze_cigar(alignment.cigar)

    assert result.operations == [
      {10, "M"}, {2, "I"}, {3, "D"}, {4, "S"}, {5, "H"}, {6, "N"}, {1, "="}, {1, "X"}
    ]
    assert result.aligned_ref_bases == 21  # 10M + 3D + 6N + 1= + 1X
    assert result.aligned_read_bases == 17 # 10M + 2I + 4S + 1= + 1X
    assert result.clipped_bases == 9      # 4S + 5H
    assert result.insertions == 2         # 2I
    assert result.deletions == 3          # 3D
    assert result.matches == 1            # 1=
    assert result.mismatches == 1         # 1X
    assert result.match_or_mismatch == 10 # 10M
    assert result.skipped == 6            # 6N
  end

  test "get_end_position calculates alignment end position correctly" do
    alignment = %SamParser.Alignment{
      pos: 100,
      cigar: "10M5D15M"
    }

    end_pos = SamParser.Alignment.Parser.get_end_position(alignment)
    assert end_pos == 129  # 100 + 10 + 5 + 15 - 1 (convert to 0-based)
  end

  test "overlaps_region checks alignment region overlap correctly" do
    alignment = %SamParser.Alignment{
      pos: 100,
      cigar: "30M"
    }

    # Test various overlap scenarios
    assert SamParser.Alignment.Parser.overlaps_region?(alignment, 90, 110) == true   # Overlap at start
    assert SamParser.Alignment.Parser.overlaps_region?(alignment, 110, 120) == true  # Overlap in middle
    assert SamParser.Alignment.Parser.overlaps_region?(alignment, 120, 140) == true  # Overlap at end
    assert SamParser.Alignment.Parser.overlaps_region?(alignment, 90, 140) == true   # Complete containment
    assert SamParser.Alignment.Parser.overlaps_region?(alignment, 50, 90) == false   # Before alignment
    assert SamParser.Alignment.Parser.overlaps_region?(alignment, 140, 160) == false # After alignment
  end

  test "reconstruct_read_sequence handles all CIGAR operations" do
    alignment = %SamParser.Alignment{
      seq: "ACGTACGT",
      cigar: "2M1I3M1D2M",
      pos: 1
    }

    read_seq = SamParser.Alignment.Parser.reconstruct_read_sequence(alignment)
    assert read_seq == "ACGTACGT"  # Original sequence preserved

    # Test with clipping operations
    alignment = %SamParser.Alignment{
      seq: "ACGTACGT",
      cigar: "2S4M2H",
      pos: 1
    }

    read_seq = SamParser.Alignment.Parser.reconstruct_read_sequence(alignment)
    assert read_seq == "ACGTACGT"  # Original sequence preserved with soft/hard clips
  end

  test "extract_quality_scores handles quality strings correctly" do
    alignment = %SamParser.Alignment{
      qual: "!~ABCDEF"  # Phred+33 scores: 0,93,65,66,67,68,69,70
    }

    scores = SamParser.Alignment.Parser.extract_quality_scores(alignment)
    assert scores == [0, 93, 32, 33, 34, 35, 36, 37]

    # Test with missing quality scores
    alignment = %SamParser.Alignment{qual: "*"}
    assert SamParser.Alignment.Parser.extract_quality_scores(alignment) == []
  end

  test "create_alignment_view generates correct alignment visualization" do
    alignment = %SamParser.Alignment{
      seq: "ACGTA",
      cigar: "2M1I2M",
      pos: 1
    }
    reference = "ACTGA"

    view = SamParser.Alignment.Parser.create_alignment_view(alignment, reference)
    assert is_binary(view)
    assert String.contains?(view, "Ref:")
    assert String.contains?(view, "Read:")
  end

  test "header formatting preserves all fields" do
    header = %SamParser.Header{
      hd: %{"VN" => "1.6", "SO" => "coordinate"},
      sq: [%{"SN" => "chr1", "LN" => "1000"}],
      rg: [%{"ID" => "1", "SM" => "sample1"}],
      pg: [%{"ID" => "prog1", "PN" => "test"}],
      co: ["Test comment"]
    }

    lines = SamParser.format_header(header)
    assert length(lines) == 5
    assert Enum.any?(lines, &(&1 =~ "@HD" && &1 =~ "VN:1.6" && &1 =~ "SO:coordinate"))
    assert Enum.any?(lines, &(&1 =~ "@SQ" && &1 =~ "SN:chr1" && &1 =~ "LN:1000"))
    assert Enum.any?(lines, &(&1 =~ "@RG" && &1 =~ "ID:1" && &1 =~ "SM:sample1"))
    assert Enum.any?(lines, &(&1 =~ "@PG" && &1 =~ "ID:prog1" && &1 =~ "PN:test"))
    assert Enum.any?(lines, &(&1 == "@CO\tTest comment"))
  end

  test "flag checking functions work correctly" do
    alignment = %SamParser.Alignment{flag: 0x3}  # paired (0x1) + proper_pair (0x2)
    assert SamParser.is_paired?(alignment) == true
    assert SamParser.is_properly_paired?(alignment) == true
    assert SamParser.is_mapped?(alignment) == true
    assert SamParser.is_reverse?(alignment) == false

    alignment = %SamParser.Alignment{flag: 0x14}  # unmapped (0x4) + reversed (0x10)
    assert SamParser.is_mapped?(alignment) == false
    assert SamParser.is_reverse?(alignment) == true
    assert SamParser.is_paired?(alignment) == false
  end

  test "parse_tag_value handles different tag types" do
    # Test integer tags
    assert SamParser.parse_tag_value("i", "42") == 42
    assert SamParser.parse_tag_value("I", "42") == 42

    # Test float tags
    assert SamParser.parse_tag_value("f", "3.14") == 3.14

    # Test string tags
    assert SamParser.parse_tag_value("Z", "hello") == "hello"
    assert SamParser.parse_tag_value("H", "4142") == "4142"

    # Test array tags
    assert SamParser.parse_tag_value("B", "i,1,2,3") == [1, 2, 3]
    assert SamParser.parse_tag_value("B", "f,1.1,2.2") == [1.1, 2.2]
  end

  test "parse_sam handles invalid files" do
    assert_raise RuntimeError, fn -> SamParser.parse_sam("nonexistent.sam") end

    # Test invalid alignment line
    alignment_line = "insufficient\tfields"
    assert_raise RuntimeError, fn -> SamParser.parse_alignment(alignment_line) end
  end

  test "filtering functions work correctly" do
    # Create test alignments
    alignments = [
      %SamParser.Alignment{rname: "chr1", pos: 100},
      %SamParser.Alignment{rname: "chr1", pos: 200},
      %SamParser.Alignment{rname: "chr2", pos: 150}
    ]
    sam_file = %SamParser.SamFile{alignments: alignments}

    # Test filter_by_reference
    chr1_result = SamParser.filter_by_reference(sam_file, "chr1")
    assert length(chr1_result.alignments) == 2
    assert Enum.all?(chr1_result.alignments, &(&1.rname == "chr1"))

    # Test filter_by_position
    pos_result = SamParser.filter_by_position(sam_file, 120, 180)
    assert length(pos_result.alignments) == 1
    assert Enum.all?(pos_result.alignments, &(&1.pos >= 120 and &1.pos <= 180))
  end

  test "writing SAM files preserves data" do
    # Create a test SAM file structure
    header = %SamParser.Header{
      hd: %{"VN" => "1.6"},
      sq: [%{"SN" => "chr1", "LN" => "1000"}]
    }

    alignment = %SamParser.Alignment{
      qname: "read1",
      flag: 0,
      rname: "chr1",
      pos: 1,
      mapq: 60,
      cigar: "10M",
      seq: "ACGTACGTAC",
      qual: "!!!!!!!!!!",
      tags: %{"NM" => {"i", 0}}
    }

    sam_file = %SamParser.SamFile{
      header: header,
      alignments: [alignment]
    }

    # Write to a temporary file
    path = "test/temp_test.sam"
    SamParser.write_sam(sam_file, path)

    # Read it back and compare
    read_sam = SamParser.parse_sam(path)
    assert read_sam.header.hd == header.hd
    assert read_sam.header.sq == header.sq

    read_aln = hd(read_sam.alignments)
    assert read_aln.qname == alignment.qname
    assert read_aln.flag == alignment.flag
    assert read_aln.seq == alignment.seq
    assert read_aln.tags["NM"] == {"i", 0}

    # Clean up
    File.rm!(path)
  end

  test "format_tag_value handles different tag types" do
    assert SamParser.format_tag_value("i", 42) == "42"
    assert SamParser.format_tag_value("f", 3.14) == "3.14"
    assert SamParser.format_tag_value("Z", "hello") == "hello"
    assert SamParser.format_tag_value("B", [1, 2, 3]) == "i,1,2,3"
    assert SamParser.format_tag_value("B", [1.1, 2.2]) == "f,1.1,2.2"
  end

  test "array type inference works correctly" do
    assert SamParser.infer_array_type([1]) == {"c", "int8"}
    assert SamParser.infer_array_type([1.5]) == {"f", "float"}
    assert SamParser.infer_array_type([]) == {"i", "int32"}  # Default
    assert SamParser.infer_array_type([127]) == {"c", "int8"}
    assert SamParser.infer_array_type([255]) == {"C", "uint8"}
    assert SamParser.infer_array_type([32767]) == {"s", "int16"}
    assert SamParser.infer_array_type([65535]) == {"S", "uint16"}
  end

  test "header validation and conversion" do
    # Test empty header
    header = %SamParser.Header{}
    assert header.hd == nil
    assert header.sq == []
    assert header.rg == []
    assert header.pg == []
    assert header.co == []

    # Test header with all fields populated
    header = %SamParser.Header{
      hd: %{"VN" => "1.6", "SO" => "coordinate"},
      sq: [
        %{"SN" => "chr1", "LN" => "1000"},
        %{"SN" => "chr2", "LN" => "500"}
      ],
      rg: [
        %{"ID" => "sample1", "SM" => "test1"},
        %{"ID" => "sample2", "SM" => "test2"}
      ],
      pg: [%{"ID" => "prog1", "PN" => "test"}],
      co: ["Comment 1", "Comment 2"]
    }

    # Test header formatting
    lines = SamParser.format_header(header)
    assert length(lines) == 8  # 1 HD + 2 SQ + 2 RG + 1 PG + 2 CO
    assert Enum.any?(lines, &(&1 =~ "@HD"))
    assert Enum.count(lines, &(&1 =~ "@SQ")) == 2
    assert Enum.count(lines, &(&1 =~ "@RG")) == 2
    assert Enum.count(lines, &(&1 =~ "@CO")) == 2

    # Test parsing header back
    parsed = SamParser.parse_header(lines)
    assert parsed.hd["VN"] == "1.6"
    assert parsed.hd["SO"] == "coordinate"
    assert length(parsed.sq) == 2
    assert length(parsed.rg) == 2
    assert length(parsed.co) == 2
  end

  test "parse_header handles malformed lines" do
    malformed_lines = [
      "@HD",  # Missing fields
      "@SQ\tSN:chr1",  # Missing LN
      "@RG\tID:1\tSM:",  # Empty value
      "@CO"  # Empty comment
    ]

    header = SamParser.parse_header(malformed_lines)
    assert header.hd == %{}  # Empty map for malformed HD line
    assert header.sq == [%{"SN" => "chr1"}]  # Partial SQ info preserved
    assert header.rg == [%{"ID" => "1", "SM" => ""}]  # Empty values preserved
    assert header.co == [""]  # Empty comment preserved
  end

  test "alignment validation and edge cases" do
    # Test default alignment
    alignment = %SamParser.Alignment{}
    assert alignment.qname == nil
    assert alignment.flag == 0
    assert alignment.rname == "*"
    assert alignment.pos == 0
    assert alignment.mapq == 0
    assert alignment.cigar == "*"
    assert alignment.rnext == "*"
    assert alignment.pnext == 0
    assert alignment.tlen == 0
    assert alignment.seq == "*"
    assert alignment.qual == "*"
    assert alignment.tags == %{}

    # Test with all fields populated
    alignment = %SamParser.Alignment{
      qname: "read1",
      flag: 0x1,
      rname: "chr1",
      pos: 100,
      mapq: 60,
      cigar: "10M",
      rnext: "=",
      pnext: 200,
      tlen: 100,
      seq: "ACGTACGTAC",
      qual: "!!!!!!!!!!",
      tags: %{"NM" => {"i", 0}}
    }

    assert alignment.qname == "read1"
    assert alignment.flag == 0x1
    assert alignment.rname == "chr1"
    assert alignment.pos == 100
    assert alignment.mapq == 60
    assert alignment.cigar == "10M"
    assert alignment.rnext == "="
    assert alignment.pnext == 200
    assert alignment.tlen == 100
    assert alignment.seq == "ACGTACGTAC"
    assert alignment.qual == "!!!!!!!!!!"
    assert alignment.tags["NM"] == {"i", 0}

    # Test flag operations
    assert SamParser.is_paired?(alignment) == true
    assert SamParser.is_mapped?(alignment) == true
    assert SamParser.is_properly_paired?(alignment) == false
    assert SamParser.is_reverse?(alignment) == false
    assert SamParser.is_secondary?(alignment) == false
    assert SamParser.is_supplementary?(alignment) == false
  end
end
