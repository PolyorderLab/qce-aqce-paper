import csv
import hashlib
import importlib.util
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from PIL import Image


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "scripts/write_bvk2_si_convergence_morphology.py"
SPEC = importlib.util.spec_from_file_location("si_convergence_morphology", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SIConvergenceMorphologyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.states = MODULE.select_states()
        cls.by_phase = {row["phase"]: row for row in cls.states}

    def test_states_are_canonical_accepted_endpoints(self):
        self.assertEqual(set(self.by_phase), {"LAM", "CYL", "BCC", "GYR"})
        for row in self.states:
            self.assertEqual(row["endpoint_status"], "accepted")
            self.assertEqual(row["phase_status"], "accepted")
            self.assertEqual(row["shell_pass"], "true")
            self.assertEqual(len(row["canonical_endpoint_source_sha256"]), 64)
            self.assertLess(abs(row["mean_error"]), 2.0e-10)
            self.assertGreater(row["contrast"], 0.99)

    def test_selection_uses_shared_g_pocket_coordinate_and_target_bcc(self):
        for phase in ("LAM", "CYL", "GYR"):
            row = self.by_phase[phase]
            self.assertEqual(row["chiN"], 54.0)
            self.assertEqual(row["fA"], 0.363)
        bcc = self.by_phase["BCC"]
        self.assertEqual(bcc["transition"], "S/C")
        self.assertEqual(bcc["chiN"], 47.5)
        self.assertEqual(bcc["fA"], 0.2085)
        self.assertLess(bcc["cell_stress"], 1.0e-3)

    def test_audit_energy_and_trace_energy_are_kept_distinct(self):
        lam = self.by_phase["LAM"]
        self.assertGreater(abs(lam["audit_minus_trace_energy_density"]), 1.0e-3)
        self.assertEqual(lam["energy_audit_factor"], "3")
        self.assertLess(
            self.by_phase["BCC"]["selected_plateau_energy_gap_density"], 1.0e-6
        )

    def test_writer_is_deterministic_and_outputs_valid_packet(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            out1, out2 = root / "one", root / "two"
            states, convergence, morphology = MODULE.write_outputs(out1)
            MODULE.write_outputs(out2)

            expected = {
                "representative_states.csv",
                "convergence_traces.csv",
                "morphology_sections.csv",
                "bvk2_si_convergence_morphology.svg",
                "bvk2_si_convergence_morphology.png",
                "README.md",
            }
            self.assertEqual({path.name for path in out1.iterdir()}, expected)
            for name in expected:
                self.assertEqual(digest(out1 / name), digest(out2 / name))

            self.assertEqual(len(states), 4)
            self.assertGreater(len(convergence), 20)
            self.assertEqual(len(morphology), 1024 + 96 * 168 + 48**2 + 112**2)
            ET.parse(out1 / "bvk2_si_convergence_morphology.svg")
            with Image.open(out1 / "bvk2_si_convergence_morphology.png") as preview:
                self.assertEqual(preview.size, (1500, 1030))
                self.assertEqual(preview.mode, "RGB")
            svg = (out1 / "bvk2_si_convergence_morphology.svg").read_text()
            for phase in ("LAM", "CYL", "BCC", "GYR"):
                self.assertIn(phase, svg)
            self.assertIn("diagnostic", svg)

            with (out1 / "representative_states.csv").open(
                newline="", encoding="utf-8"
            ) as handle:
                public = list(csv.DictReader(handle))
            self.assertTrue(all(row["claim_kind"] == "phase_endpoint" for row in public))
            self.assertTrue(all(row["phase_status"] == "accepted" for row in public))


if __name__ == "__main__":
    unittest.main()
