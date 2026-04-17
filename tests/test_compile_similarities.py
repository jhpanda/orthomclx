from __future__ import annotations

import struct
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from orthomcl.compile_similarities import RECORD_STRUCT, compile_similarities


def write_file(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


class CompileSimilaritiesTest(unittest.TestCase):
    def test_compile_similarities_writes_integer_coded_binary(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            similar_sequences = tmp_path / "similarSequences.txt"
            compiled_dir = tmp_path / "compiled"
            write_file(
                similar_sequences,
                "\n".join(
                    [
                        "a1\tb1\taaa\tbbb\t1\t-50\t95.0\t90.0",
                        "b1\ta1\tbbb\taaa\t1\t-50\t95.0\t90.0",
                        "a2\tb2\taaa\tbbb\t2\t-10\t80.0\t50.0",
                    ]
                )
                + "\n",
            )

            summary = compile_similarities(similar_sequences, compiled_dir)

            self.assertEqual(summary.record_count, 3)
            self.assertEqual(summary.protein_count, 4)
            self.assertEqual(summary.taxon_count, 2)
            self.assertEqual(
                summary.proteins_path.read_text(),
                "0\ta1\n1\tb1\n2\ta2\n3\tb2\n",
            )
            self.assertEqual(
                summary.taxa_path.read_text(),
                "0\taaa\n1\tbbb\n",
            )

            data = summary.binary_path.read_bytes()
            self.assertEqual(len(data), 3 * RECORD_STRUCT.size)
            first = RECORD_STRUCT.unpack(data[: RECORD_STRUCT.size])
            self.assertEqual(first[:4], (0, 1, 0, 1))
            self.assertAlmostEqual(first[4], 1.0)
            self.assertEqual(first[5], -50)
            self.assertAlmostEqual(first[6], 95.0)
            self.assertAlmostEqual(first[7], 90.0)

    def test_compile_similarities_rewrites_zero_evalues_to_underflow_exponent(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            similar_sequences = tmp_path / "similarSequences.txt"
            compiled_dir = tmp_path / "compiled"
            write_file(
                similar_sequences,
                "\n".join(
                    [
                        "a1\tb1\taaa\tbbb\t0\t0\t95.0\t90.0",
                        "b1\ta1\tbbb\taaa\t0\t0\t95.0\t90.0",
                        "a2\tb2\taaa\tbbb\t5\t-20\t80.0\t50.0",
                    ]
                )
                + "\n",
            )

            summary = compile_similarities(similar_sequences, compiled_dir)

            data = summary.binary_path.read_bytes()
            first = RECORD_STRUCT.unpack(data[: RECORD_STRUCT.size])
            second = RECORD_STRUCT.unpack(data[RECORD_STRUCT.size : 2 * RECORD_STRUCT.size])
            third = RECORD_STRUCT.unpack(data[2 * RECORD_STRUCT.size : 3 * RECORD_STRUCT.size])
            self.assertEqual(first[5], -21)
            self.assertEqual(second[5], -21)
            self.assertEqual(third[5], -20)
