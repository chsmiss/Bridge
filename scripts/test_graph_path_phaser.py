#!/usr/bin/env python3
import gzip
import tempfile
import unittest
from pathlib import Path

import graph_path_phaser as gp


class GraphPhaserTest(unittest.TestCase):
    def test_full_read_context_selects_supported_branch(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            gfa = td / "g.gfa"
            gfa.write_text(
                "H\tVN:Z:1.0\n"
                "S\tu0\t" + "A" * 31 + "C" * 31 + "\tLN:i:62\tKC:f:10\n"
                "S\tu1\t" + "C" * 31 + "G" * 31 + "\tLN:i:62\tKC:f:10\n"
                "S\tu2\t" + "C" * 31 + "T" * 31 + "\tLN:i:62\tKC:f:3\n"
                "S\tu3\t" + "G" * 31 + "A" * 31 + "\tLN:i:62\tKC:f:10\n"
                "L\tu0\t+\tu1\t+\t31M\tDR:i:5\tGR:i:0\tPE:i:0\n"
                "L\tu0\t+\tu2\t+\t31M\tDR:i:2\tGR:i:0\tPE:i:0\n"
                "L\tu1\t+\tu3\t+\t31M\tDR:i:5\tGR:i:0\tPE:i:0\n"
            )
            reads = td / "r.fq.gz"
            seq = "A" * 31 + "C" * 31 + "G" * 31 + "A" * 31
            with gzip.open(reads, "wt") as out:
                for i in range(5):
                    record = seq[i : i + 101]
                    out.write(
                        f"@r{i}\n{record}\n+\n" + "I" * len(record) + "\n"
                    )
            graph = gp.Graph.from_gfa(gfa)
            index = gp.KmerIndex(graph, 31)
            raw, _ = gp.collect_read_contexts(graph, index, reads, None, None, 6)
            paths, _ = gp.resolve_paths(
                graph,
                raw,
                gp.Counter(),
                gp.Counter(),
                gp.Counter(),
                0.72,
                4,
                50,
            )
            self.assertTrue(any(path == ["u0", "u1", "u3"] for path in paths))
            self.assertFalse(
                any(path[:2] == ["u0", "u2"] for path in paths if len(path) > 1)
            )


if __name__ == "__main__":
    unittest.main()
