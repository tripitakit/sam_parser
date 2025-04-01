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
end
