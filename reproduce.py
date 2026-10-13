#!/usr/bin/env python3
"""Reproduce the computational results reported in the article."""

import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT / "examples"), str(ROOT / "research")]

from automatic_recurrences import (
    OriginalLabyrinthArray,
    OriginalNegativeArray,
    learn_kernel,
    verify_recurrence,
)
from compute_examples import minimize, rational_automaton
from png_read import read_rgb


def check_discovery():
    predecessors = [(1, 0), (1, 2), (2, -2)]
    negative = OriginalNegativeArray(forcing=1)
    rejected = []
    for window, witness in [(1, (2, 2)), (2, (4, 4)), (3, (52, 16))]:
        candidate = learn_kernel(negative, 2, 2, window)
        verdict = verify_recurrence(
            candidate, predecessors, [1, 1, 1], 2, forcing=1
        )
        assert not verdict["valid"]
        assert tuple(verdict["witness"]) == witness
        rejected.append({"window": window, "witness": list(witness)})

    recovered = []
    for prefix, provider, window, states in [
        ("negative_forced", negative, 4, 53),
        ("labyrinth", OriginalLabyrinthArray(), 6, 254),
    ]:
        candidate = json.loads(json.dumps(learn_kernel(provider, 2, 2, window)))
        stored = json.loads(
            (ROOT / "research" / "generated" / f"{prefix}_machine.json").read_text()
        )
        assert len(candidate["outputs"]) == states
        for key in ("base", "dimension", "digits", "initial", "transitions", "outputs"):
            assert candidate[key] == stored[key], (prefix, key)
        recovered.append({"example": prefix, "window": window, "states": states})

    print("Binary samples reproduce the three failing indices and both stored machines.", flush=True)
    return {"rejected_negative_samples": rejected, "recovered_machines": recovered}


def check_small_delannoy_powers():
    results = []
    for exponent, expected in [(1, 3), (2, 58)]:
        machine = rational_automaton(3, exponent)
        classes, rounds = minimize(machine)
        count = len(set(classes))
        assert count == expected, (exponent, count)
        results.append({
            "modulus": 3**exponent,
            "polynomial_states": len(machine["outputs"]),
            "minimal_lsd_states": count,
            "minimization_rounds": rounds,
        })
    print("Delannoy polynomial identities and minimal state counts modulo 3 and 9 verified.", flush=True)
    return results


def main():
    if not __debug__:
        raise SystemExit("Run Python with assertions enabled, without -O or -OO.")

    results = {
        "binary_discovery": check_discovery(),
        "delannoy_small_powers": check_small_delannoy_powers(),
    }
    for script in (
        "research/certify_original_examples.py",
        "research/verify_mod27.py",
        "examples/nonabelian_examples.py",
    ):
        print(f"Running {script}", flush=True)
        subprocess.run([sys.executable, str(ROOT / script)], cwd=ROOT, check=True)

    binary = json.loads((ROOT / "research/generated/original_examples_report.json").read_text())
    for example, expected in (
        ("negative", (53, 107, 44, 44, 39, 686)),
        ("labyrinth", (254, 710, 204, 212, 151, 552493)),
    ):
        record = binary[example]
        actual = tuple(record[key] for key in (
            "minimal_lsd_states", "minimal_msd_states", "numerical_block_rules",
            "colored_substitution_states", "distinct_rgb_colors",
        )) + (record["recurrence_certificate"]["reachable_product_states"],)
        assert actual == expected, (example, actual)

    group_counts = json.loads(
        (ROOT / "examples/generated/nonabelian/verification.json").read_text()
    )
    assert group_counts["genuine_homomorphism"]["distinct_values"] == 39
    assert group_counts["negative_ordered_product"]["distinct_values"] == 27
    assert group_counts["negative_ordered_product"]["cells_changed_by_reversing_multiplication"] == 17737

    for name in ("hom_central.png", "hom_delannoy.png", "negative_central.png", "negative_reversed.png"):
        generated = ROOT / "examples" / "generated" / "nonabelian" / name
        assert read_rgb(generated) == read_rgb(ROOT / "images" / name), name

    results["paper_images_match_regenerated_pixels"] = True
    destination = ROOT / "research" / "generated" / "reproduction_results.json"
    destination.write_text(json.dumps(results, indent=2) + "\n")
    print("All reported computational checks passed.", flush=True)


if __name__ == "__main__":
    main()
