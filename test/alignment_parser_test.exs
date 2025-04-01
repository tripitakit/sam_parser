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

  describe "extract_reference_sequence with independent validation" do
    setup do
      # Create a reference sequence with distinct bases at each position to easily verify extraction
      reference = "ACGTACGTACGTACGTACGT" # 20-base reference with repeating pattern

      # Create a map showing the 0-based index of each position to aid debugging
      # Index:  0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9
      # Base:   A C G T A C G T A C G T A C G T A C G T

      %{reference: reference}
    end

    test "simple match extraction - 10M", %{reference: reference} do
      alignment = %Alignment{
        pos: 2, # 1-based position (C at index 1)
        cigar: "10M"
      }

      # Expected: Extract 10 bases starting from position 2 (1-based)
      # Reference indices 1-10 (0-based): CGTACGTACG
      expected = String.slice(reference, 1, 10)

      result = Parser.extract_reference_sequence(alignment, reference)
      assert result == expected
      assert result == "CGTACGTACG"
    end

    test "extraction with deletions - 3M2D5M", %{reference: reference} do
      alignment = %Alignment{
        pos: 2, # 1-based position (C at index 1)
        cigar: "3M2D5M"
      }

      # Manual calculation:
      # 1. 3M: Take 3 bases from reference at position 2 (1-based)
      #    => Indices 1-3 (0-based): "CGT"
      # 2. 2D: Skip 2 bases in reference (positions 5-6, 1-based)
      #    => Skip indices 4-5 (0-based): "AC"
      # 3. 5M: Take 5 bases from reference at position 7 (1-based)
      #    => Indices 6-10 (0-based): "GTACG"
      # Result: "CGT" + "GTACG" = "CGTGTACG"

      first_part = String.slice(reference, 1, 3)  # CGT
      second_part = String.slice(reference, 6, 5) # GTACG
      expected = first_part <> second_part

      result = Parser.extract_reference_sequence(alignment, reference)
      assert result == expected
      assert result == "CGTGTACG"
    end

    test "extraction with skipped regions - 3M2N5M", %{reference: reference} do
      alignment = %Alignment{
        pos: 2, # 1-based position (C at index 1)
        cigar: "3M2N5M"
      }

      # Manual calculation:
      # 1. 3M: Take 3 bases from reference at position 2 (1-based)
      #    => Indices 1-3 (0-based): "CGT"
      # 2. 2N: Add "NN" for skipped region, but keep advancing in reference
      #    => Skip indices 4-5 (0-based): "AC" (replaced with "NN")
      # 3. 5M: Take 5 bases from reference at position 7 (1-based)
      #    => Indices 6-10 (0-based): "GTACG"
      # Result: "CGT" + "NN" + "GTACG" = "CGTNNGTACG"

      first_part = String.slice(reference, 1, 3)  # CGT
      skipped = "NN"
      third_part = String.slice(reference, 6, 5)  # GTACG
      expected = first_part <> skipped <> third_part

      result = Parser.extract_reference_sequence(alignment, reference)
      assert result == expected
      assert result == "CGTNNGTACG"
    end

    test "extraction with insertions - 3M2I5M", %{reference: reference} do
      alignment = %Alignment{
        pos: 2, # 1-based position (C at index 1)
        cigar: "3M2I5M"
      }

      # Manual calculation:
      # 1. 3M: Take 3 bases from reference at position 2 (1-based)
      #    => Indices 1-3 (0-based): "CGT"
      # 2. 2I: Insertions don't consume reference bases,
      #    => No change to reference position
      # 3. 5M: Take 5 bases from reference at position 5 (1-based)
      #    => Indices 4-8 (0-based): "ACGTA"
      # Result: "CGT" + "ACGTA" = "CGTACGTA"

      first_part = String.slice(reference, 1, 3)  # CGT
      second_part = String.slice(reference, 4, 5) # ACGTA
      expected = first_part <> second_part

      result = Parser.extract_reference_sequence(alignment, reference)
      assert result == expected
      assert result == "CGTACGTA"
    end

    test "extraction with exact matches and mismatches - 3=2X5=", %{reference: reference} do
      alignment = %Alignment{
        pos: 2, # 1-based position (C at index 1)
        cigar: "3=2X5="
      }

      # Manual calculation:
      # For reference extraction, =, X, and M all consume the reference in the same way
      # 1. 3=: Take 3 bases from reference at position 2 (1-based)
      #    => Indices 1-3 (0-based): "CGT"
      # 2. 2X: Take 2 bases from reference at position 5 (1-based)
      #    => Indices 4-5 (0-based): "AC"
      # 3. 5=: Take 5 bases from reference at position 7 (1-based)
      #    => Indices 6-10 (0-based): "GTACG"
      # Result: "CGT" + "AC" + "GTACG" = "CGTACGTACG"

      first_part = String.slice(reference, 1, 3)  # CGT
      second_part = String.slice(reference, 4, 2) # AC
      third_part = String.slice(reference, 6, 5)  # GTACG
      expected = first_part <> second_part <> third_part

      result = Parser.extract_reference_sequence(alignment, reference)
      assert result == expected
      assert result == "CGTACGTACG"
    end
  end

  describe "extract_reference_sequence edge cases" do
    setup do
      %{reference: "ACGTACGTACGTACGTACGT"}
    end

    test "handles invalid CIGAR strings", %{reference: reference} do
      # Empty CIGAR string
      alignment = %Alignment{
        pos: 2,
        cigar: "*"
      }
      assert Parser.extract_reference_sequence(alignment, reference) == ""

      # Invalid CIGAR operation
      alignment = %Alignment{
        pos: 2,
        cigar: "3M2Z5M" # 'Z' is not a valid CIGAR operation
      }
      assert_raise ArgumentError, fn -> Parser.extract_reference_sequence(alignment, reference) end

      # Invalid CIGAR format (missing number before operation)
      alignment = %Alignment{
        pos: 2,
        cigar: "3MM5M"
      }
      assert_raise ArgumentError, fn -> Parser.extract_reference_sequence(alignment, reference) end
    end

    test "handles out-of-bounds reference access", %{reference: reference} do
      # Position too large
      alignment = %Alignment{
        pos: 30, # Outside of reference bounds
        cigar: "5M"
      }
      assert_raise ArgumentError, fn -> Parser.extract_reference_sequence(alignment, reference) end

      # Position + CIGAR extends past reference
      alignment = %Alignment{
        pos: 18, # Within bounds, but CIGAR extends past end
        cigar: "10M"
      }
      # Should either raise or truncate, depending on implementation
      result = Parser.extract_reference_sequence(alignment, reference)
      assert String.length(result) <= 3 # At most can return "CGT" (positions 18-20)
    end

    test "handles edge positions", %{reference: reference} do
      # Position at start of reference
      alignment = %Alignment{
        pos: 1,
        cigar: "5M"
      }
      result = Parser.extract_reference_sequence(alignment, reference)
      assert result == "ACGTA"

      # Position at end of reference
      alignment = %Alignment{
        pos: 20, # Last position (1-based)
        cigar: "1M"
      }
      result = Parser.extract_reference_sequence(alignment, reference)
      assert result == "T"
    end
  end

  describe "create_alignment_view with exact formatting validation" do
    setup do
      reference = "ACGTACGTACGTACGTACGT"
      %{reference: reference}
    end

    test "validates simple match alignment (5M)", %{reference: reference} do
      alignment = %Alignment{
        pos: 2,
        cigar: "5M",
        seq: "CGTAC"
      }

      expected = "Ref:  CGTAC\n      |||||\nRead: CGTAC\n"
      result = Parser.create_alignment_view(alignment, reference)

      assert result == expected
    end

    test "validates insertion alignment (3M2I3M)", %{reference: reference} do
      alignment = %Alignment{
        pos: 2,
        cigar: "3M2I3M",
        seq: "CGTTTACG"  # CGT + TT + ACG
      }

      # Manually calculated:
      # 3M: CGT matches the reference exactly
      # 2I: TT is inserted in read but not in reference (represented as --)
      # 3M: ACG matches the reference exactly
      expected = "Ref:  CGT--ACG\n      |||  |||\nRead: CGTTTACG\n"
      result = Parser.create_alignment_view(alignment, reference)

      assert result == expected
    end

    test "validates deletion alignment (3M2D3M)", %{reference: reference} do
      alignment = %Alignment{
        pos: 2,
        cigar: "3M2D3M",
        seq: "CGTGTA" # CGT + GTA (skipping AC in reference)
      }

      # Manually calculated:
      # 3M: CGT matches the reference exactly
      # 2D: AC in reference is deleted (represented as --)
      # 3M: GTA matches the reference exactly
      expected = "Ref:  CGTACGTA\n      |||  |||\nRead: CGT--GTA\n"
      result = Parser.create_alignment_view(alignment, reference)

      assert result == expected
    end

    test "validates soft-clipping (2S3M2S)", %{reference: reference} do
      alignment = %Alignment{
        pos: 2,
        cigar: "2S3M2S",
        seq: "TTCGTTT" # TT + CGT + TT
      }

      # Manually calculated:
      # 2S: TT is soft-clipped (not aligned to reference)
      # 3M: CGT matches the reference exactly
      # 2S: TT is soft-clipped (not aligned to reference)
      expected = "Ref:    CGT  \n      |||\nRead: TTCGTTT\n"
      result = Parser.create_alignment_view(alignment, reference)

      assert result == expected
    end

    test "validates complex alignment with multiple operations", %{reference: reference} do
      alignment = %Alignment{
        pos: 2,
        cigar: "2M1I1D2M1X1=",
        seq: "CGTATGT" # CG + T + - + AC + G + T
      }

      # Manually calculated:
      # 2M: CG matches the reference exactly
      # 1I: T is inserted in read
      # 1D: T in reference is deleted
      # 2M: AC matches the reference exactly
      # 1X: G in read is a mismatch to G in reference
      # 1=: T matches the eference exactly
      expected = "Ref:  CG-TACGT\n      ||  |  |\nRead: CGT-ATGT\n"
      result = Parser.create_alignment_view(alignment, reference)

      assert result == expected
    end
  end

  describe "create_alignment_view edge cases" do
    setup do
      reference = "ACGTACGTACGTACGTACGT"
      %{reference: reference}
    end

    test "handles missing sequence", %{reference: reference} do
      alignment = %Alignment{
        pos: 2,
        cigar: "3M",
        seq: "*"
      }
      result = Parser.create_alignment_view(alignment, reference)
      assert result == "No sequence available for alignment view"
    end

    test "handles empty CIGAR string", %{reference: reference} do
      alignment = %Alignment{
        pos: 2,
        cigar: "*",
        seq: "CGT"
      }
      result = Parser.create_alignment_view(alignment, reference)
      assert result == "No CIGAR string available for alignment view"
    end

    test "handles sequence and CIGAR length mismatch", %{reference: reference} do
      # Sequence shorter than needed for CIGAR
      alignment = %Alignment{
        pos: 2,
        cigar: "5M",
        seq: "CGT" # Only 3 bases for 5M
      }
      assert_raise ArgumentError, fn -> Parser.create_alignment_view(alignment, reference) end

      # Sequence longer than needed for CIGAR
      alignment = %Alignment{
        pos: 2,
        cigar: "3M",
        seq: "CGTAC" # 5 bases for 3M
      }
      assert_raise ArgumentError, fn -> Parser.create_alignment_view(alignment, reference) end
    end

    test "handles out-of-bounds reference access", %{reference: reference} do
      alignment = %Alignment{
        pos: 18,
        cigar: "5M", # Goes beyond reference bounds
        seq: "CGTAC"
      }
      assert_raise ArgumentError, fn -> Parser.create_alignment_view(alignment, reference) end
    end
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
  end
end
