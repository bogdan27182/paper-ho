#!/usr/bin/env python3
"""Exact group recurrences and figures. Python standard library only.

Images are raw data: n increases to the right, m increases downwards.
The finite checks verify the implementations and concrete witnesses;
automaticity is proved in the accompanying article.
"""
from itertools import product
from pathlib import Path
import json
import struct
import zlib

OUT = Path(__file__).resolve().parent / "generated" / "nonabelian"
P = 3
IDENTITY = (0, 0, 0)
X, Y = (1, 0, 0), (0, 1, 0)
PALETTE3 = ((247, 245, 239), (40, 116, 142), (204, 91, 70))
PALETTE9 = (
    (247, 245, 239), (168, 211, 213), (44, 123, 149),
    (232, 199, 132), (145, 151, 107), (59, 85, 109),
    (221, 135, 96), (169, 77, 76), (73, 43, 72),
)


def mul(g, h):
    x, y, z = g[:3]
    u, v, w = h[:3]
    first = ((x + u) % P, (y + v) % P, (z + w + x * v) % P)
    return first if len(g) == 3 else first + ((g[3] + h[3]) % P,)


def inverse(g):
    x, y, z = g[:3]
    first = (-x % P, -y % P, (x * y - z) % P)
    return first if len(g) == 3 else first + (-g[3] % P,)


def alpha(g):
    x, y, z, t = g
    return y, x, (x * y - z) % P, t


def beta(g):
    return 0, 0, g[3], g[3]


def tau(g):
    x, y, z = g
    return x, y, (x * y - z) % P


def hom_rule(g, h, k):
    return mul(mul(alpha(g), beta(h)), beta(k))


def genuine_hom_array(size):
    a = [[None] * size for _ in range(size)]
    for n in range(size):
        a[0][n] = (n % P, 1, 0, 1)
    for m in range(1, size):
        a[m][0] = alpha(a[m - 1][0])
        for n in range(1, size):
            a[m][n] = hom_rule(a[m - 1][n], a[m][n - 1], a[m - 1][n - 1])
    return a


def negative_array(size, reverse=False):
    # Row m needs entries through n=size-1 from row m-1 through n=size.
    # This triangular extension prevents truncation from changing the square.
    a = [[IDENTITY] + [Y] * (2 * size - 2)]
    for m in range(1, size):
        row = [X]
        for n in range(1, 2 * size - 1 - m):
            left, above = a[m - 1][n + 1], row[n - 1]
            row.append(mul(above, left) if reverse else mul(left, above))
        a.append(row)
    return [row[:size] for row in a]


def negative_by_weight(size):
    # Independently ordered evaluation on 2m+n <= 3(size-1).
    # Each predecessor has smaller weight; there is no row-wise extension.
    bound = 3 * (size - 1)
    values = {}
    for weight in range(bound + 1):
        for m in range(weight // 2 + 1):
            n = weight - 2 * m
            if m == 0:
                values[m, n] = Y if n else IDENTITY
            elif n == 0:
                values[m, n] = X
            else:
                values[m, n] = mul(values[m - 1, n + 1], values[m, n - 1])
    return [[values[m, n] for n in range(size)] for m in range(size)]


def png(path, rows, palette, scale=3):
    width, height = len(rows[0]) * scale, len(rows) * scale
    raw = bytearray()
    for row in rows:
        pixels = b"".join(bytes(palette[v]) * scale for v in row)
        raw.extend((b"\0" + pixels) * scale)

    def chunk(tag, data):
        return struct.pack(">I", len(data)) + tag + data + struct.pack(">I", zlib.crc32(tag + data))

    path.write_bytes(
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
        + chunk(b"IDAT", zlib.compress(raw, 9))
        + chunk(b"IEND", b"")
    )


def commutator(g, h):
    return mul(mul(mul(g, h), inverse(g)), inverse(h))


def count_values(a):
    return len({v for row in a for v in row})


def check_groups():
    h3 = list(product(range(P), repeat=3))
    g81 = list(product(range(P), repeat=4))
    for g in h3:
        assert mul(g, inverse(g)) == IDENTITY
        for h in h3:
            assert tau(mul(g, h)) == mul(tau(h), tau(g))
            for k in h3:
                assert mul(mul(g, h), k) == mul(g, mul(h, k))
    for g in g81:
        assert alpha(alpha(g)) == g
        assert mul(g, inverse(g)) == (0, 0, 0, 0)
        for h in g81:
            assert alpha(mul(g, h)) == mul(alpha(g), alpha(h))
            assert beta(mul(g, h)) == mul(beta(g), beta(h))
            assert mul(g, beta(h)) == mul(beta(h), g)
    # The preceding identities imply hom_rule is a homomorphism G^3 -> G.
    return {"H3_associativity_triples": len(h3) ** 3,
            "endomorphism_pairs_per_map": len(g81) ** 2,
            "central_image_pairs": len(g81) ** 2}


def check_hom_array(a):
    size = len(a)
    # Independent scalar decomposition plus the closed Delannoy digit formula.
    digit = ((1, 1, 1), (1, 0, 2), (1, 2, 1))
    for m in range(size):
        for n in range(size):
            x, y, z, t = a[m][n]
            assert (x, y) == ((n % 3, 1) if m % 2 == 0 else (1, n % 3))
            u, v, value = m, n, 1
            while u or v:
                value = value * digit[u % 3][v % 3] % 3
                u //= 3
                v //= 3
            assert t == value
            if m and n:
                assert z == (n - a[m - 1][n][2] + a[m][n - 1][3] + a[m - 1][n - 1][3]) % 3
    g, h = a[1][1], a[1][2]
    assert commutator(g, h) == (0, 0, 1, 0)
    assert count_values(a) > 27
    return {"size": size, "distinct_values": count_values(a),
            "interior_noncommuting_positions": [[1, 1], [1, 2]],
            "values": [g, h], "commutator": commutator(g, h),
            "central_nonzero_cells": sum(v[2] != 0 for row in a for v in row)}


def check_negative(a, rev):
    size = len(a)
    independent = negative_by_weight(35)
    assert [row[:35] for row in a[:35]] == independent
    changed = 0
    for m in range(size):
        for n in range(size):
            assert rev[m][n] == tau(a[m][n])
            changed += a[m][n] != rev[m][n]
    assert a[2][2] == (0, 2, 2)
    assert commutator(a[1][1], a[1][2]) == (0, 0, 1)
    return {"size": size, "distinct_values": count_values(a),
            "independent_weight_order_square": 35,
            "cells_changed_by_reversing_multiplication": changed,
            "nonzero_central_witness": {"position": [2, 2], "value": a[2][2]},
            "interior_noncommuting_positions": [[1, 1], [1, 2]],
            "values": [a[1][1], a[1][2]],
            "commutator": commutator(a[1][1], a[1][2])}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    checks = {"group_checks": check_groups()}
    size = 243
    hom = genuine_hom_array(size)
    signed = negative_array(size)
    reverse = negative_array(size, reverse=True)
    checks["genuine_homomorphism"] = check_hom_array(hom)
    checks["negative_ordered_product"] = check_negative(signed, reverse)
    images = {
        "hom_central": ([[g[2] for g in row] for row in hom], PALETTE3),
        "hom_delannoy": ([[g[3] for g in row] for row in hom], PALETTE3),
        "negative_quotient": ([[3*g[0]+g[1] for g in row] for row in signed], PALETTE9),
        "negative_central": ([[g[2] for g in row] for row in signed], PALETTE3),
        "negative_reversed": ([[g[2] for g in row] for row in reverse], PALETTE3),
        "negative_order_effect": ([[int(g != h) for g, h in zip(row, other)]
                                    for row, other in zip(signed, reverse)], PALETTE3[:2]),
    }
    for name, (rows, palette) in images.items():
        png(OUT / (name + ".png"), rows, palette)
    checks["palette3"] = PALETTE3
    checks["coordinate_convention"] = "m down; n right; 0 through 242 inclusive"
    (OUT / "verification.json").write_text(json.dumps(checks, indent=2) + "\n")
    (OUT / "small_arrays.json").write_text(json.dumps({
        "genuine_homomorphism": [r[:10] for r in hom[:10]],
        "negative_ordered_product": [r[:10] for r in signed[:10]],
    }, indent=2) + "\n")
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
