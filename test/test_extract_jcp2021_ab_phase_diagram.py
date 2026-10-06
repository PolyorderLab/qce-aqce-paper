import csv
import importlib.util
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "results" / "jcp2021_ab_phase_diagram_digitization" / "extract_from_eps.py"


def load_extractor():
    spec = importlib.util.spec_from_file_location("extract_jcp2021_ab", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_vector_extraction_retains_complete_ab_topology():
    module = load_extractor()
    postscript = module.extract_postscript(module.DEFAULT_SOURCE)
    paths = module.parse_top_panel_paths(postscript)
    segments = module.build_segments(paths)

    assert len(segments) == 18
    assert {segment.transition for segment in segments} == {
        "DIS/S_cp",
        "S_cp/S",
        "DIS/S",
        "S/C",
        "C/G",
        "C/O70",
        "G/O70",
        "G/L",
        "O70/L",
    }
    for transition in module.COLORS:
        assert {segment.side for segment in segments if segment.transition == transition} == {"left", "right"}


def test_topology_nodes_match_vector_source_coordinates():
    module = load_extractor()
    nodes = {(row["node"], row["side"]): row for row in module.topology_nodes()}

    critical = nodes[("critical", "center")]
    assert critical["f_A"] == 0.5
    assert abs(critical["chiN"] - 10.537542662116) < 1.0e-12

    sphere_left = nodes[("sphere_triple", "left")]
    assert abs(sphere_left["f_A"] - 0.237241379310) < 1.0e-12
    assert abs(sphere_left["chiN"] - 17.747440273038) < 1.0e-12


def test_checked_output_has_small_independent_overlay_residual():
    audit = SCRIPT.parent / "jcp2021_ab_extraction_audit.csv"
    assert audit.exists()
    with audit.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 18
    assert max(float(row["raster_residual_p90_px"]) for row in rows) <= 3.0
