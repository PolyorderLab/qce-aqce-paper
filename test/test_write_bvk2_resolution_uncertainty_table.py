import copy
import importlib.util
import math
import tempfile
import unittest
from collections import Counter
from pathlib import Path
from unittest import mock

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "scripts/write_bvk2_resolution_uncertainty_table.py"
SPEC = importlib.util.spec_from_file_location("resolution_table", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class ResolutionUncertaintyTableTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.points = MODULE.build_pointwise_table()
        cls.contracts = MODULE.summarize_contracts(cls.points)

    def filtered_extra_reader(self, f_a: float, delta_f: float):
        original = MODULE.read_rows
        endpoint_path = (
            MODULE.RESULTS / "bvk2_ud_theta_gc_filtered_r6/endpoint_pairs.csv"
        )
        manifest_path = MODULE.RESULTS / "bvk2_ud_theta_gc_filtered_r6/manifest.csv"
        endpoints = original(endpoint_path)
        manifests = original(manifest_path)
        template = next(
            row
            for row in endpoints
            if row["chiN"] == "30.0" and row["fA"] == "0.351592857"
        )
        template_manifest = next(
            row for row in manifests if row["job_id"] == template["job_id"]
        )
        synthetic_artifact = (
            MODULE.resolve_repository_path(template_manifest["outdir"])
            / f"synthetic_fA{f_a:.9f}.csv"
        )
        synthetic_job = f"synthetic-chi30-fA{f_a:.9f}"
        endpoint = copy.deepcopy(template)
        endpoint.update(
            {
                "artifact": str(synthetic_artifact),
                "fA": str(f_a),
                "job_id": synthetic_job,
                "delta_f": str(delta_f),
            }
        )
        cyl_energy = float(endpoint["cyl_energy_density"])
        endpoint["gyr_energy_density"] = str(cyl_energy + delta_f)
        manifest = copy.deepcopy(template_manifest)
        manifest.update(
            {
                "artifact": str(synthetic_artifact),
                "fA": str(f_a),
                "job_id": synthetic_job,
            }
        )
        raw = copy.deepcopy(original(MODULE.resolve_repository_path(template["artifact"]))[0])
        raw.update(
            {
                "fA": str(f_a),
                "delta_f": str(delta_f),
                "gyr_energy_density": endpoint["gyr_energy_density"],
            }
        )

        def altered(path):
            if path.resolve() == synthetic_artifact.resolve():
                return [copy.deepcopy(raw)]
            rows = copy.deepcopy(original(path))
            if path.resolve() == endpoint_path.resolve():
                rows.append(copy.deepcopy(endpoint))
            elif path.resolve() == manifest_path.resolve():
                rows.append(copy.deepcopy(manifest))
            return rows

        return altered

    def test_all_canonical_points_are_covered(self):
        self.assertEqual(len(self.points), 106)
        self.assertEqual(
            Counter(row["transition"] for row in self.points),
            {"S/DIS": 33, "S/C": 27, "G/C": 22, "G/L": 24},
        )
        self.assertTrue(
            all(row["canonical_status"] == "accepted" for row in self.points)
        )
        self.assertTrue(
            all(not Path(row["canonical_source"]).is_absolute() for row in self.points)
        )

    def test_source_routing_and_canonical_matching_fail_closed(self):
        direct_gl_sources = {
            MODULE.resolve_repository_path(row["canonical_source"])
            for row in self.points
            if row["audit_status"] == "current_direct_contract"
        }
        for row in self.points:
            family = MODULE.classify_canonical_source(
                row["transition"],
                row["canonical_source"],
                direct_gl_sources=direct_gl_sources,
            )
            self.assertIsInstance(family, str)
        with self.assertRaisesRegex(ValueError, "unrecognized canonical source"):
            MODULE.classify_canonical_source(
                "S/C",
                "results/unrelated/boundary.csv",
                direct_gl_sources=direct_gl_sources,
            )
        with self.assertRaisesRegex(ValueError, "canonical source mismatch"):
            MODULE.require_same_source(
                "results/a.csv", "results/b.csv", context="mutation"
            )
        legacy = PROJECT_ROOT / "results/bvk2_gc_gl_boundary_campaign"
        self.assertEqual(
            MODULE.resolve_repository_path(MODULE.relative(legacy)), legacy.resolve()
        )

    def test_current_stress_validator_rejects_missing_nonfinite_and_above_limit(self):
        self.assertAlmostEqual(
            MODULE.validate_current_stress([1e-4, -2e-3], context="valid"),
            2e-3,
        )
        for mutation in ([], [""], [None], [float("nan")], [float("inf")], [0.00201]):
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                MODULE.validate_current_stress(mutation, context="mutation")

    def test_richardson_match_requires_chin_coordinate_and_active_provenance(self):
        publication = {
            "transition": "S/C",
            "chiN": "12.0",
            "fA": "0.3",
            "source": "results/bvk2_gc_gl_boundary_campaign",
        }
        correction = {
            "transition": "S/C",
            "chiN": "12.0",
            "root_fA_corrected": "0.3",
        }
        provenance = {
            "transition": "S/C",
            "chiN": "12.0",
            "fA_corrected": "0.3",
            "status": "accepted",
            "gate_reason": "accepted",
            "source": "results/bvk2_gc_gl_boundary_campaign",
        }
        path = PROJECT_ROOT / "results/bvk2_richardson_correction/test.csv"
        self.assertEqual(
            MODULE.match_correction(publication, [(correction, path)], [provenance]),
            (correction, path),
        )
        wrong_chi = dict(correction, chiN="13.0")
        with self.assertRaisesRegex(ValueError, "uniquely match canonical correction"):
            MODULE.match_correction(publication, [(wrong_chi, path)], [provenance])
        with self.assertRaisesRegex(ValueError, "active Richardson provenance"):
            MODULE.match_correction(publication, [(correction, path)], [provenance, provenance])
        wrong_source = dict(provenance, source="results/unrelated")
        with self.assertRaisesRegex(ValueError, "active Richardson provenance"):
            MODULE.match_correction(publication, [(correction, path)], [wrong_source])

    def test_full_build_rejects_nonfinite_active_correction_values(self):
        original = MODULE.read_rows

        for nonfinite in ("NaN", "Inf", "-Inf"):
            def altered(path, replacement=nonfinite):
                rows = copy.deepcopy(original(path))
                if path.name == "corrected_roots.csv":
                    for row in rows:
                        if row["transition"] == "S/C" and row["chiN"] == "12.0":
                            row["shift"] = replacement
                return rows

            with self.subTest(nonfinite=nonfinite), mock.patch.object(
                MODULE, "read_rows", altered
            ), self.assertRaisesRegex(ValueError, "shift is nonfinite"):
                MODULE.build_pointwise_table()

    def test_full_build_rejects_duplicate_filtered_gc_identity(self):
        original = MODULE.read_rows
        endpoint_path = (
            MODULE.RESULTS / "bvk2_ud_theta_gc_filtered_r6/endpoint_pairs.csv"
        )

        def altered(path):
            rows = copy.deepcopy(original(path))
            if path.resolve() == endpoint_path.resolve():
                rows.append(copy.deepcopy(rows[0]))
            return rows

        with mock.patch.object(
            MODULE, "read_rows", altered
        ), self.assertRaisesRegex(ValueError, "identities are not unique"):
            MODULE.build_pointwise_table()

    def test_real_filtered_gc_ledger_validates_all_29_endpoints_and_off_bracket_rows(self):
        root = MODULE.RESULTS / "bvk2_ud_theta_gc_filtered_r6"
        endpoints = MODULE.read_rows(root / "endpoint_pairs.csv")
        selected = MODULE.validate_filtered_gc_contract(
            MODULE.read_rows(root / "boundary_points.csv"),
            endpoints,
            MODULE.read_rows(root / "manifest.csv"),
        )
        self.assertEqual(len(endpoints), 29)
        self.assertEqual(len(selected), 9)
        self.assertEqual(
            [float(row["fA"]) for row in selected[30.0]],
            [0.351592857, 0.353592857],
        )

    def test_fully_provenanced_same_sign_off_bracket_endpoint_is_allowed(self):
        altered = self.filtered_extra_reader(0.357592857, -0.0015)
        with mock.patch.object(MODULE, "read_rows", altered):
            self.assertEqual(len(MODULE.build_pointwise_table()), 106)

    def test_filtered_gc_stale_delta_and_secant_are_rejected(self):
        original = MODULE.read_rows
        boundary_path = (
            MODULE.RESULTS / "bvk2_ud_theta_gc_filtered_r6/boundary_points.csv"
        )
        for field, replacement in (
            ("delta_lower", "0.123"),
            ("coordinate", "0.352"),
        ):
            def altered(path, key=field, value=replacement):
                rows = copy.deepcopy(original(path))
                if path.resolve() == boundary_path.resolve():
                    rows[0][key] = value
                return rows

            with self.subTest(field=field), mock.patch.object(
                MODULE, "read_rows", altered
            ), self.assertRaisesRegex(ValueError, "is stale"):
                MODULE.build_pointwise_table()

    def test_valid_narrower_interior_sign_change_makes_summary_stale(self):
        altered = self.filtered_extra_reader(0.352592857, -0.0001)
        with mock.patch.object(
            MODULE, "read_rows", altered
        ), self.assertRaisesRegex(ValueError, "is stale"):
            MODULE.build_pointwise_table()

    def test_filtered_gc_missing_and_duplicate_manifest_joins_fail(self):
        original = MODULE.read_rows
        manifest_path = MODULE.RESULTS / "bvk2_ud_theta_gc_filtered_r6/manifest.csv"
        for mode in ("missing", "duplicate"):
            def altered(path, mutation=mode):
                rows = copy.deepcopy(original(path))
                if path.resolve() == manifest_path.resolve():
                    if mutation == "missing":
                        rows.pop()
                    else:
                        rows.append(copy.deepcopy(rows[0]))
                return rows

            with self.subTest(mode=mode), mock.patch.object(
                MODULE, "read_rows", altered
            ), self.assertRaisesRegex(ValueError, "manifest"):
                MODULE.build_pointwise_table()

    def test_filtered_gc_raw_contract_and_provenance_mutations_fail(self):
        original = MODULE.read_rows
        endpoint_path = (
            MODULE.RESULTS / "bvk2_ud_theta_gc_filtered_r6/endpoint_pairs.csv"
        )
        endpoint = original(endpoint_path)[0]
        raw_path = MODULE.resolve_repository_path(endpoint["artifact"])
        for field, replacement in (
            ("schema", "wrong-schema"),
            ("gyr_grid", "64x64x64"),
            ("sensor_filter_ratio", "5"),
            ("status", "rejected"),
            ("gyr_cell_root", endpoint["cyl_cell_root"]),
        ):
            def altered(path, key=field, value=replacement):
                rows = copy.deepcopy(original(path))
                if path.resolve() == raw_path.resolve():
                    rows[0][key] = value
                return rows

            with self.subTest(field=field), mock.patch.object(
                MODULE, "read_rows", altered
            ), self.assertRaises(ValueError):
                MODULE.build_pointwise_table()

        def artifact_mismatch(path):
            rows = copy.deepcopy(original(path))
            if path.resolve() == endpoint_path.resolve():
                rows[0]["artifact"] = rows[1]["artifact"]
            return rows

        with mock.patch.object(
            MODULE, "read_rows", artifact_mismatch
        ), self.assertRaisesRegex(ValueError, "artifact"):
            MODULE.build_pointwise_table()

    def test_filtered_gc_tied_minimum_sign_changes_fail(self):
        altered = self.filtered_extra_reader(0.349592857, -0.0001)
        with mock.patch.object(
            MODULE, "read_rows", altered
        ), self.assertRaisesRegex(ValueError, "is tied"):
            MODULE.build_pointwise_table()

    def test_endpoint_coordinate_sets_and_gl15_record_are_unique(self):
        rows = [{"fA": "0.2"}, {"fA": "0.3"}]
        selected = MODULE.exact_coordinate_rows(
            rows,
            coordinates=(0.2, 0.3),
            coordinate_key="fA",
            per_coordinate=1,
            context="mutation",
        )
        self.assertEqual(len(selected), 2)
        with self.assertRaisesRegex(ValueError, "unique rows"):
            MODULE.exact_coordinate_rows(
                rows + [{"fA": "0.2"}],
                coordinates=(0.2, 0.3),
                coordinate_key="fA",
                per_coordinate=1,
                context="mutation",
                reject_extra=True,
            )
        with self.assertRaisesRegex(ValueError, "coordinate set is not exact"):
            MODULE.exact_coordinate_rows(
                rows + [{"fA": "0.4"}],
                coordinates=(0.2, 0.3),
                coordinate_key="fA",
                per_coordinate=1,
                context="mutation",
                reject_extra=True,
            )
        with self.assertRaisesRegex(ValueError, "expected one"):
            MODULE.unique_record([{"chiN": "15"}, {"chiN": "15"}], context="G/L15")

    def test_identity_coverage_rejects_duplicates_missing_and_extra_rows(self):
        publication = {
            "transition": "G/L",
            "chiN": "15.0",
            "fA": "0.4",
            "source": "results/source.csv",
        }
        generated = {
            "transition": "G/L",
            "chiN": 15.0,
            "fA": 0.4,
            "canonical_source": "results/source.csv",
        }
        MODULE.validate_identity_coverage([publication], [generated])
        with self.assertRaisesRegex(ValueError, "identities are not unique"):
            MODULE.validate_identity_coverage([publication], [generated, generated])
        mutated = dict(generated, fA=0.41)
        with self.assertRaisesRegex(ValueError, "one-to-one"):
            MODULE.validate_identity_coverage([publication], [mutated])

    def test_contract_summary_rejects_any_invariant_mutation(self):
        contract_id = self.points[0]["contract_id"]
        same_contract = [
            copy.deepcopy(row)
            for row in self.points
            if row["contract_id"] == contract_id
        ]
        self.assertGreater(len(same_contract), 1)
        for field in (
            "phase_a_resolution",
            "phase_b",
            "energy_evaluation",
            "resolution_treatment",
            "uncertainty_scope",
            "stress_policy",
        ):
            mutated = copy.deepcopy(same_contract)
            mutated[0][field] = f"mutated-{field}"
            with self.subTest(field=field), self.assertRaisesRegex(
                ValueError, "not homogeneous"
            ):
                MODULE.summarize_contracts(mutated)

    def test_sdis_direction_and_uncertainty_units_are_explicit(self):
        rows = [row for row in self.points if row["transition"] == "S/DIS"]
        self.assertEqual(
            Counter(row["root_direction"] for row in rows),
            {"fA": 24, "chiN": 9},
        )
        chi_roots = [row for row in rows if row["root_direction"] == "chiN"]
        self.assertTrue(
            all(
                math.isclose(row["coordinate_half_width"], 0.00375, abs_tol=1e-12)
                for row in chi_roots
            )
        )
        self.assertTrue(
            all(row["root_semantics"].startswith("ordered-branch") for row in rows)
        )

    def test_legacy_richardson_uncertainty_is_not_overclaimed(self):
        corrected = [row for row in self.points if row["richardson_used"] == "true"]
        self.assertEqual(len(corrected), 19)
        self.assertTrue(
            all(row["strict_total_coordinate_bound"] == "false" for row in corrected)
        )
        shifted = [
            row
            for row in corrected
            if row["audit_status"]
            == "resolution_shift_exceeds_stored_half_width"
        ]
        self.assertEqual(len(shifted), 14)
        self.assertEqual(
            sum(row["transition"] == "G/C" for row in shifted), 8
        )

    def test_new_source_families_use_current_grids_and_provenance(self):
        fixed_sc = [
            row
            for row in self.points
            if row["contract_id"] == "SC_FIXED_FA_BCC48_CYL48x84"
        ]
        self.assertEqual(len(fixed_sc), 6)
        self.assertTrue(all(row["root_direction"] == "chiN" for row in fixed_sc))
        self.assertTrue(
            all("bvk2_sc_fixed_fa_chi_continuation" in row["canonical_source"] for row in fixed_sc)
        )

        filtered_gc = [
            row
            for row in self.points
            if row["contract_id"] == "GC_FILTERED_R6_GYR96_CYL128x224"
        ]
        self.assertEqual(len(filtered_gc), 9)
        self.assertTrue(all(row["phase_a_resolution"] == "96x96x96" for row in filtered_gc))
        self.assertTrue(all(row["phase_b_resolution"] == "128x224" for row in filtered_gc))
        self.assertTrue(
            all("bvk2_ud_theta_gc_filtered_r6" in row["evidence_sources"] for row in filtered_gc)
        )

        matched_gl = [
            row
            for row in self.points
            if row["audit_status"] == "current_unified_matched_spacing_contract"
        ]
        self.assertEqual(len(matched_gl), 20)
        self.assertTrue(
            all("bvk2_ud_theta_gl_matched_spacing_audit" in row["canonical_source"] for row in matched_gl)
        )

    def test_later_uniform_stress_audit_is_bounded(self):
        failures = [
            row
            for row in self.points
            if row["stress_tier_under_current_contract"] == "above_current_maximum"
        ]
        self.assertEqual(len(failures), 10)
        self.assertEqual(
            {(row["transition"], row["chiN"]) for row in failures},
            {
                ("S/C", 40.0), ("S/C", 45.0), ("S/C", 49.0),
                ("S/C", 52.0), ("S/C", 55.0), ("S/C", 58.0),
                ("S/C", 60.0), ("G/C", 52.0), ("G/C", 55.0),
                ("G/C", 60.0),
            },
        )

    def test_current_gl_contract_is_direct_and_within_stress_maximum(self):
        rows = [row for row in self.points if row["transition"] == "G/L"]
        self.assertEqual(len(rows), 24)
        self.assertTrue(all(row["richardson_used"] == "false" for row in rows))
        self.assertTrue(
            all(
                row["stress_tier_under_current_contract"] in {"target", "fallback"}
                for row in rows
            )
        )
        self.assertEqual(
            {row["phase_b_resolution"] for row in rows}, {"26", "40", "46", "1024"}
        )

    def test_writer_emits_pointwise_and_contract_tables(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            points, contracts = MODULE.write_outputs(output)
            self.assertEqual(len(points), 106)
            self.assertGreaterEqual(len(contracts), 10)
            self.assertTrue(
                (output / "phase_boundary_resolution_uncertainty.csv").is_file()
            )
            self.assertTrue((output / "resolution_contracts.csv").is_file())
            self.assertTrue((output / "resolution_contract_map.svg").is_file())
            self.assertTrue((output / "acceptance_rejection_audit.svg").is_file())
            resolution_svg = (output / "resolution_contract_map.svg").read_text()
            acceptance_svg = (output / "acceptance_rejection_audit.svg").read_text()
            self.assertIn("A-block fraction", resolution_svg)
            self.assertIn("attempted bracket width", acceptance_svg)
            self.assertIn("Matplotlib + mpltex ACS via uv", resolution_svg)
            self.assertIn("Matplotlib + mpltex ACS via uv", acceptance_svg)
            self.assertNotIn(
                "Resolution contracts for accepted BVK2 phase-boundary points",
                resolution_svg,
            )
            self.assertNotIn(
                "Acceptance-status and rejected-candidate audit",
                acceptance_svg,
            )
            readme = (output / "README.md").read_text()
            self.assertIn("automatically a total continuum uncertainty", readme)
            self.assertIn("close the manuscript's numerical-certification gate", readme)
            self.assertIn("106 canonical accepted", readme)
            self.assertIn("19 legacy Richardson-corrected", readme)
            self.assertIn("10** direct", readme)


if __name__ == "__main__":
    unittest.main()
