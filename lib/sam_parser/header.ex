defmodule SamParser.Header do
  @moduledoc "Structure for storing SAM header information"
  defstruct hd: nil, sq: [], rg: [], pg: [], co: []

  @type t() :: %__MODULE__{
    hd: map() | nil,      # Header line information
    sq: list(map()),      # Reference sequence dictionary
    rg: list(map()),      # Read groups
    pg: list(map()),      # Programs
    co: list(String.t())  # Comments
  }
end
