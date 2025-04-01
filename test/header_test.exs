defmodule SamParser.HeaderTest do
  use ExUnit.Case

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
end
