defmodule SamParser.SamFile do
  @moduledoc "Structure for storing a complete SAM file with header and alignments"
  alias SamParser.{Header, Alignment}
  defstruct header: %SamParser.Header{}, alignments: []

  @type t() :: %__MODULE__{
    header: Header.t(),
    alignments: list(Alignment.t())
  }
end
