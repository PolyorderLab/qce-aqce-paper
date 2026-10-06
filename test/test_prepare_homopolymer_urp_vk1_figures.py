#!/usr/bin/env python3
"""Regression checks for the homopolymer URP/VK1 publication figures."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
RESULT = PROJECT / "results/homopolymer_urp_vk1_validation"

for name in ("urp_vk1_error_map.svg", "urp_vk1_saddle_profile.svg"):
    path = RESULT / name
    assert path.is_file(), path
    root = ET.parse(path).getroot()
    assert root.tag.endswith("svg")
    view_box = [float(value) for value in root.attrib["viewBox"].split()]
    assert view_box[2] > 200 and view_box[3] > 150, (name, view_box)

readme = (RESULT / "README.md").read_text(encoding="utf-8")
assert "auxiliary-potential profiles" in readme
assert "not\nindependently relaxed density profiles" in readme
for relative in (
    "source/accuracy_saddle_approx/"
    "accuracy_saddle_approx-eps-converted-to.pdf",
    "source/accuracy_saddle_approx/error_K_L_a0.01.dat",
    "source/accuracy_saddle_approx/error_K_L_a0.8.dat",
    "source/accuracy_saddle_approx/error_K_a_L0.5.dat",
    "source/accuracy_saddle_approx/error_K_a_L8.dat",
    "source/accuracy_saddle_approx/plot.py",
    "source/iw_saddle_approx/iw_saddle_approx-eps-converted-to.pdf",
    "source/iw_saddle_approx/phitype5_L3_a0.5_b1_c1.mat",
    "source/iw_saddle_approx/input_phitype5.txt",
    "source/iw_saddle_approx/plot.py",
):
    assert (RESULT / relative).is_file(), relative

print("homopolymer URP/VK1 validation figures passed")
