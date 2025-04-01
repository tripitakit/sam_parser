defmodule SamParserTest do
  use ExUnit.Case
  doctest SamParser

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
