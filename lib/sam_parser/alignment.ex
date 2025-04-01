defmodule SamParser.Alignment do
  @moduledoc "Structure for storing SAM alignment records"
  defstruct [
    # Mandatory fields
    qname: nil,   # Query template NAME
    flag: 0,      # Bitwise FLAG
    rname: "*",   # Reference sequence NAME
    pos: 0,       # 1-based leftmost mapping POSition
    mapq: 0,      # MAPping Quality
    cigar: "*",   # CIGAR string
    rnext: "*",   # Ref. name of the mate/next read
    pnext: 0,     # Position of the mate/next read
    tlen: 0,      # Observed Template LENgth
    seq: "*",     # Segment SEQuence
    qual: "*",    # ASCII of Phred-scaled base QUALity+33
    # Optional field storage
    tags: %{}     # Map to store TAG:TYPE:VALUE optional fields
  ]
end