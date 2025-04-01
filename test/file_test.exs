defmodule SamParser.FileTest do
  use ExUnit.Case

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
end
