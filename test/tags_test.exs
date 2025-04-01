defmodule SamParser.TagsTest do
  use ExUnit.Case

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
end
