defmodule SamParser.Alignment.ParserTest do
  use ExUnit.Case
  alias SamParser.Alignment
  alias SamParser.Alignment.Parser

  test "build_flag creates correct flag values from flag map" do
    # Test with all flags false
    all_false = %{
      paired: false,
      proper_pair: false,
      unmapped: false,
      next_unmapped: false,
      reversed: false,
      next_reversed: false,
      first: false,
      last: false,
      secondary: false,
      filtered: false,
      duplicate: false,
      supplementary: false
    }
    assert Parser.build_flag(all_false) == 0x0

    # Test with all flags true
    all_true = %{
      paired: true,
      proper_pair: true,
      unmapped: true,
      next_unmapped: true,
      reversed: true,
      next_reversed: true,
      first: true,
      last: true,
      secondary: true,
      filtered: true,
      duplicate: true,
      supplementary: true
    }
    assert Parser.build_flag(all_true) == 0xFFF

    # Test specific combinations
    paired_reads = %{
      paired: true,
      proper_pair: true,
      unmapped: false,
      next_unmapped: false,
      reversed: false,
      next_reversed: false,
      first: true,
      last: false,
      secondary: false,
      filtered: false,
      duplicate: false,
      supplementary: false
    }
    assert Parser.build_flag(paired_reads) == 0x43 # 0x1 + 0x2 + 0x40

    # Test unmapped read
    unmapped_read = %{
      paired: true,
      proper_pair: false,
      unmapped: true,
      next_unmapped: false,
      reversed: false,
      next_reversed: false,
      first: false,
      last: true,
      secondary: false,
      filtered: false,
      duplicate: false,
      supplementary: false
    }
    assert Parser.build_flag(unmapped_read) == 0x85 # 0x1 + 0x4 + 0x80
  end

  test "extract_reference_sequence creates correct sequence with various CIGAR operations" do
    reference = "ACTGACTGACTGACTGACTG" # 20-base reference sequence

    # Simple match - extract 10 bases from reference at position 1 (0-based)
    alignment = %Alignment{
      pos: 2, # 1-based position
      cigar: "10M"
    }
    assert Parser.extract_reference_sequence(alignment, reference) == "CTGACTGACT"

    # Test with deletions - should be included in reference
    alignment = %Alignment{
      pos: 2,
      cigar: "3M2D5M" # CTGACTGACT with bases 5-6 deleted
    }
    assert Parser.extract_reference_sequence(alignment, reference) == "CTGTGACT"

    # Test with insertions - should not be in reference sequence
    alignment = %Alignment{
      pos: 2,
      cigar: "3M2I5M" # Insertions don't affect reference
    }
    assert Parser.extract_reference_sequence(alignment, reference) == "CTGACTGA"

    # Test with skipped regions (N)
    alignment = %Alignment{
      pos: 2,
      cigar: "3M2N5M" # CTG...ACTGA with 2 bases skipped (replaced with N)
    }
    assert Parser.extract_reference_sequence(alignment, reference) == "CTGNNTGACT"

    # Test with exact matches (=) and mismatches (X)
    alignment = %Alignment{
      pos: 2,
      cigar: "3=2X5=" # Same as 10M for reference extraction
    }
    assert Parser.extract_reference_sequence(alignment, reference) == "CTGACTGACT"
  end

  test "create_alignment_view handles CIGAR operations" do
    reference = "ACTGACTGACTGACTGACTG" # 20-base reference sequence

    # Test matching bases - "M" CIGAR operation
    alignment = %Alignment{
      pos: 2,
      cigar: "5M",
      seq: "CTGAC"
    }
    view = Parser.create_alignment_view(alignment, reference)
    assert String.contains?(view, "Ref:  CTGAC")
    assert String.contains?(view, "     |||||")
    assert String.contains?(view, "Read: CTGAC")

    # Test with insertions - "I" CIGAR operation
    alignment = %Alignment{
      pos: 2,
      cigar: "3M2I2M",
      seq: "CTGTTAC"
    }
    view = Parser.create_alignment_view(alignment, reference)
    assert String.contains?(view, "Ref:  CTG--AC")
    assert String.contains?(view, "     |||  ||")
    assert String.contains?(view, "Read: CTGTTAC")

    # Test with missing sequence
    alignment = %Alignment{
      pos: 2,
      cigar: "3M",
      seq: "*"
    }
    view = Parser.create_alignment_view(alignment, reference)
    assert view == "No sequence available for alignment view"
  end

  test "analyze_cigar handles edge cases" do
    # Test with empty CIGAR string
    assert Parser.analyze_cigar("*") == %{
      operations: [],
      op_counts: %{},
      aligned_ref_bases: 0,
      aligned_read_bases: 0,
      clipped_bases: 0,
      insertions: 0,
      deletions: 0,
      matches: 0,
      mismatches: 0,
      match_or_mismatch: 0,
      skipped: 0
    }

    # Test with a complex CIGAR string that includes all operation types
    cigar = "5M2I3D1S2H4N2=3X"
    result = Parser.analyze_cigar(cigar)

    assert result.operations == [
      {5, "M"}, {2, "I"}, {3, "D"}, {1, "S"}, {2, "H"}, {4, "N"}, {2, "="}, {3, "X"}
    ]
    assert result.aligned_ref_bases == 17  # 5M + 3D + 4N + 2= + 3X
    assert result.aligned_read_bases == 13 # 5M + 2I + 1S + 2= + 3X
    assert result.clipped_bases == 3      # 1S + 2H
    assert result.insertions == 2         # 2I
    assert result.deletions == 3          # 3D
    assert result.matches == 2            # 2=
    assert result.mismatches == 3         # 3X
    assert result.match_or_mismatch == 5  # 5M
    assert result.skipped == 4            # 4N

    # Test just a simple alignment
    result = Parser.analyze_cigar("10M")
    assert result.operations == [{10, "M"}]
    assert result.op_counts == %{"M" => 10}
    assert result.aligned_ref_bases == 10
    assert result.aligned_read_bases == 10
    assert result.match_or_mismatch == 10
  end

  test "interpret_flags and build_flag are inverse operations" do
    # Test with various flag combinations
    flag_values = [0x0, 0x1, 0x3, 0x5, 0x9, 0x11, 0x43, 0x83, 0x103, 0x203, 0x403, 0x803, 0xFFF]

    for flag_value <- flag_values do
      alignment = %Alignment{flag: flag_value}
      flag_map = Parser.interpret_flags(alignment)
      rebuilt_flag = Parser.build_flag(flag_map)
      assert rebuilt_flag == flag_value, "Flag #{flag_value} should round-trip correctly"
    end
  end

  test "overlaps_region? handles different alignment positions" do
    # Test alignment at beginning
    alignment = %Alignment{
      pos: 1,
      cigar: "10M"
    }
    assert Parser.overlaps_region?(alignment, 1, 5) == true
    assert Parser.overlaps_region?(alignment, 5, 15) == true
    assert Parser.overlaps_region?(alignment, 11, 20) == false

    # Test alignment in middle
    alignment = %Alignment{
      pos: 50,
      cigar: "20M"
    }
    assert Parser.overlaps_region?(alignment, 40, 49) == false
    assert Parser.overlaps_region?(alignment, 40, 50) == true
    assert Parser.overlaps_region?(alignment, 60, 70) == true
    assert Parser.overlaps_region?(alignment, 70, 80) == false

    # Test with complex CIGAR - alignment spans 100-119 (10M), then 5 bases deleted (not in read),
    # then spans 125-134 (10M)
    alignment = %Alignment{
      pos: 100,
      cigar: "10M5D10M"
    }
    assert Parser.overlaps_region?(alignment, 90, 99) == false
    assert Parser.overlaps_region?(alignment, 90, 100) == true
    assert Parser.overlaps_region?(alignment, 120, 125) == true  # Overlaps deletion
    # The alignment end position is 100 + 10 + 5 + 10 - 1 = 124
    assert Parser.overlaps_region?(alignment, 124, 140) == true
    assert Parser.overlaps_region?(alignment, 125, 140) == false

    # Test exact boundary cases
    alignment = %Alignment{
      pos: 100,
      cigar: "10M"
    }
    assert Parser.overlaps_region?(alignment, 100, 100) == true  # Left boundary
    assert Parser.overlaps_region?(alignment, 109, 109) == true  # Right boundary
    assert Parser.overlaps_region?(alignment, 110, 110) == false # Just outside right boundary
    assert Parser.overlaps_region?(alignment, 99, 99) == false   # Just outside left boundary
  end

  test "extract_quality_scores handles various quality strings" do
    # Test with range of ASCII quality values
    alignment = %Alignment{qual: "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"}
    scores = Parser.extract_quality_scores(alignment)

    # Phred+33 encoding means '!' is 0, '"' is 1, etc.
    expected_scores = Enum.to_list(0..41)
    assert scores == expected_scores

    # Test with single quality score
    alignment = %Alignment{qual: "I"}  # ASCII 73, Phred+33 = 40
    assert Parser.extract_quality_scores(alignment) == [40]

    # Test with empty quality string
    alignment = %Alignment{qual: ""}
    assert Parser.extract_quality_scores(alignment) == []
  end

  test "reconstruct_read_sequence with different CIGAR operations" do
    # Basic test - should just return the sequence
    alignment = %Alignment{
      seq: "ACTGACTGAC",
      cigar: "10M"
    }
    assert Parser.reconstruct_read_sequence(alignment) == "ACTGACTGAC"

    # Test with complex CIGAR
    alignment = %Alignment{
      seq: "ACTGTTGAC",
      cigar: "4M2I3M"  # Insertion in middle
    }
    assert Parser.reconstruct_read_sequence(alignment) == "ACTGTTGAC"

    # Test with clipping
    alignment = %Alignment{
      seq: "ACTGACTGAC",
      cigar: "2S6M2S"  # Soft-clipped at both ends
    }
    assert Parser.reconstruct_read_sequence(alignment) == "ACTGACTGAC"
  end
end
