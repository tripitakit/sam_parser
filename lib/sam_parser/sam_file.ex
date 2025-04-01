defmodule SamParser.SamFile do
  @moduledoc "Structure for storing a complete SAM file with header and alignments"
  defstruct header: %SamParser.Header{}, alignments: []
end
