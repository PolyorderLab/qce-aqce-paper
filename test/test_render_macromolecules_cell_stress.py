from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import pytest


PROJECT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT / "scripts" / "render_macromolecules_figures.py"
sys.path.insert(0, str(SCRIPT.parent))
SPEC = importlib.util.spec_from_file_location("render_macromolecules_figures", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_cell_stress_renderer_rejects_mixed_schema(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        MODULE,
        "rows",
        lambda _path: [
            {
                "schema": "lamellar-cell-stress-mechanism-v1",
                "case_id": "f0.5_chiN20",
                "model": "bvk2",
            }
        ],
    )
    with pytest.raises(ValueError, match="one v2 schema"):
        MODULE.render_cell_stress()


def test_cell_stress_renderer_rejects_in_progress_generation(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    result_dir = tmp_path / "lamellar_cell_stress_mechanism"
    result_dir.mkdir()
    (result_dir / ".promotion-in-progress").write_text(
        "interrupted-generation\n", encoding="utf-8"
    )
    monkeypatch.setattr(MODULE, "RESULTS", tmp_path)
    with pytest.raises(ValueError, match="promotion is incomplete"):
        MODULE.render_cell_stress()


def test_cell_stress_renderer_is_atomic_and_data_only() -> None:
    source = SCRIPT.read_text(encoding="utf-8")
    function = source[source.index("def render_cell_stress") : source.index(
        "def render_phase_diagram"
    )]
    assert "lamellar-cell-stress-mechanism-v2" in function
    assert "bvk2_fixed" in function
    assert '"liu2019_opf"' in function
    assert "plot_models = [" in function
    plot_model_block = function.split("plot_models = [", 1)[1].split("]", 1)[0]
    assert '"liu2019_opf"' not in plot_model_block
    assert ".replace(output)" in function
    assert "source.is_absolute()" in function
    assert ".promotion-in-progress" in function
    assert "generation_manifest.csv" in function
    assert "hashlib.sha256" in function
    assert 'axes[0, 0].set_ylabel(r"$10^2[F/V-(F/V)_{\\min}]/s_m$")' in function
    assert 'axes[1, 0].set_ylabel(r"$\\sigma/s_m$")' in function
    assert 'FormatStrFormatter("%.1f")' in function
    assert "MultipleLocator(0.5)" in function
    assert "handles, labels = axes[0, 0].get_legend_handles_labels()" in function
    assert 'axes[1, 1].legend(handles, labels, ncol=2, loc="upper left")' in function
    assert "axes[0, 0].legend(" not in function
    assert "axes[0, col].set_ylabel" not in function
    assert "axes[1, col].set_ylabel" not in function


def test_figure2_manuscript_matches_matched_control() -> None:
    manuscript = (
        PROJECT / "docs/manuscript/macromolecules/manuscript.md"
    ).read_text(encoding="utf-8")
    figure2 = manuscript[manuscript.index("### 3.2 ") : manuscript.index(
        "### 3.3 "
    )]
    figure2_flat = " ".join(figure2.split())
    assert "OPF" not in figure2
    assert "from 2.655% to 0.194%" in figure2_flat
    assert "from 4.266% to 0.711%" in figure2_flat
    assert "K_\\psi=K_{\\psi,0}" not in figure2
    assert "QCE provides a matched control" in figure2_flat
    assert "isolates the effect of AQCE's adaptive stiffness (eq 21)" in figure2_flat
    assert "symmetric state \\((f_A,\\chi N)=(0.50,20)\\)" in figure2_flat
    assert "asymmetric state \\((0.35,30)\\)" in figure2_flat
    assert "\\sigma(D)=\\left.\\partial(F/V)/\\partial\\ln D\\right|_{\\phi_D}" in figure2_flat
    assert "same period (Figures 2c and 2d)" in figure2_flat
    assert "same cell size" not in figure2_flat
    assert "Root ratios and branch validation are documented in Section S3.3" in figure2_flat
    assert "without implying uniform improvement" not in figure2_flat
    assert "subtract each model's minimum free-energy density" in figure2_flat
    assert "scale both energy and stress by the same positive factor" in figure2_flat
    assert "preserving the locations of the minima and stress zeros" in figure2_flat
    assert "not absolute free energies between functionals" not in figure2_flat

    caption = figure2.split("**Figure 2.", 1)[1].split("\n\n", 1)[0]
    assert len(caption.split()) <= 100
    assert (
        "\\(s_m=\\max\\{|F_m(D_{0,m})/V|,1\\}\\) is the within-model plotting scale"
        in " ".join(caption.split())
    )
    assert figure2.count("s_m=") == 1
    assert "stress-free ratios" not in caption
    assert "Filled symbols mark each model's" in caption
    assert "fixed-stiffness control for AQCE" not in caption
    assert "K_{\\psi,2}" not in caption

    figure3 = manuscript[manuscript.index("### 3.3 Lamellar") : manuscript.index(
        "### 3.4 Diblock"
    )]
    figure3_flat = " ".join(figure3.split())
    assert "condition-specific coefficients are determined from SCFT force and stress data" in figure3_flat
    assert "strongest literature baseline for equilibrium-period accuracy [12,21]" in figure3_flat
    assert "base nonlinear closure as the principal source of improved profile fidelity" in figure3_flat
    assert "bounded adaptive stiffness further corrects period selection" in figure3_flat

    figure2_svg = (
        PROJECT
        / "results/lamellar_cell_stress_mechanism/lamellar_cell_stress_mechanism.svg"
    ).read_text(encoding="utf-8")
    assert "<!-- OPF -->" not in figure2_svg
    assert figure2_svg.count(r"<!-- $10^2[F/V-(F/V)_{\min}]/s_m$ -->") == 1
    assert figure2_svg.count(r"<!-- $\sigma/s_m$ -->") == 1
    assert r"\partial(F/V)/\partial\ln D" not in figure2_svg


def test_figure2_uses_public_burp_name() -> None:
    result_dir = PROJECT / "results" / "lamellar_cell_stress_mechanism"
    legacy_name = "BURP" + "-TI"
    for name in ("scan.csv", "summary.csv", "README.md"):
        text = (result_dir / name).read_text(encoding="utf-8")
        assert legacy_name not in text
    readme = (result_dir / "README.md").read_text(encoding="utf-8")
    assert "Figure 2 plots UD, BURP, QCE, and AQCE" in readme
    assert "source ledger also retains OPF" in readme


def test_figure2_uses_explicit_fixed_stiffness_notation() -> None:
    assert MODULE.model_label("bvk2_fixed") == "QCE"
    assert MODULE.model_label("bvk2") == "AQCE"
    result_dir = PROJECT / "results" / "lamellar_cell_stress_mechanism"
    for name in ("scan.csv", "summary.csv", "README.md"):
        text = (result_dir / name).read_text(encoding="utf-8")
        assert "fixed Kψ" not in text
    readme = (result_dir / "README.md").read_text(encoding="utf-8")
    assert "`K_{psi,2}=K_{psi,0}`" in readme
    si = (
        PROJECT / "docs/manuscript/macromolecules/supporting_information.md"
    ).read_text(encoding="utf-8")
    assert "\\(K_{\\psi,2}=K_{\\psi,0}\\)" in si
