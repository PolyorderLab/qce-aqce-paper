import importlib.util
import math
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "scripts/write_bvk2_boundary_crossing_panels.py"
SPEC = importlib.util.spec_from_file_location("boundary_crossings", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class BoundaryCrossingPanelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.summaries, cls.endpoints = MODULE.build_dataset()
        cls.by_transition = {
            row["transition"]: row for row in cls.summaries
        }

    def test_dataset_contains_four_accepted_crossings(self):
        self.assertEqual(
            set(self.by_transition), {"S/DIS", "S/C", "G/C", "G/L"}
        )
        self.assertEqual(len(self.endpoints), 8)
        for row in self.summaries:
            self.assertEqual(row["publication_status"], "accepted")
            self.assertLessEqual(row["bracket_width"], 0.0010000001)
            self.assertLess(row["deltaF_lower"] * row["deltaF_upper"], 0.0)

    def test_roots_reproduce_canonical_sources(self):
        self.assertAlmostEqual(
            self.by_transition["S/DIS"]["root_fA"], 0.30209467256842853
        )
        self.assertAlmostEqual(
            self.by_transition["S/C"]["root_fA"], 0.207827654974695
        )
        self.assertAlmostEqual(
            self.by_transition["G/C"]["root_fA"], 0.3624211535554627
        )
        self.assertAlmostEqual(
            self.by_transition["G/L"]["root_fA"], 0.36334347773613607
        )

    def test_sdis_keeps_survival_edge_semantics(self):
        row = self.by_transition["S/DIS"]
        self.assertEqual(row["root_method"], "survival_bracket_midpoint")
        self.assertEqual(row["strict_signed_endpoint_root"], "false")
        self.assertGreater(
            abs(row["root_fA"] - row["diagnostic_secant_fA"]), 1.0e-6
        )
        self.assertLess(
            abs(row["root_fA"] - row["diagnostic_secant_fA"]),
            row["root_uncertainty"],
        )

    def test_chi54_shared_gyr_state_is_below_both_competitors(self):
        shared = [
            row
            for row in self.endpoints
            if row["transition"] in {"G/C", "G/L"}
            and math.isclose(row["coordinate"], 0.363)
        ]
        self.assertEqual(len(shared), 2)
        self.assertAlmostEqual(
            shared[0]["phase_a_energy_density"],
            shared[1]["phase_a_energy_density"],
        )
        self.assertTrue(
            all(row["deltaF_phase_a_minus_phase_b"] < 0.0 for row in shared)
        )

    def test_writer_produces_parseable_svg_and_tables(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            summaries, endpoints = MODULE.write_outputs(output)
            self.assertEqual(len(summaries), 4)
            self.assertEqual(len(endpoints), 8)
            ET.parse(output / "bvk2_boundary_crossings.svg")
            svg = (output / "bvk2_boundary_crossings.svg").read_text()
            for transition in ("S/DIS", "S/C", "G/C", "G/L"):
                self.assertIn(transition, svg)
            self.assertTrue((output / "crossing_summary.csv").is_file())
            self.assertTrue((output / "crossing_endpoints.csv").is_file())
            self.assertTrue((output / "README.md").is_file())


if __name__ == "__main__":
    unittest.main()
