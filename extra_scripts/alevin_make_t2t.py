from pathlib import Path
import sys

import pandas as pd


DEFAULT_INDEX = Path("/gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu/")


def main() -> None:
    """Create a transcript-to-transcript mapping for a spliceu index."""
    index = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_INDEX

    t2g = pd.read_csv(index / "spliceu_t2g.tsv", sep="\t", names=["t", "g"])
    t2t = t2g.copy()
    t2t.g = t2t.t
    t2t.to_csv(index / "t2t.tsv", sep="\t", header=False, index=False)


if __name__ == "__main__":
    main()
