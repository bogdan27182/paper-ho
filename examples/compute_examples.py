#!/usr/bin/env python3
"""Shared polynomial and automaton routines for the paper's verification programs.

Imported by reproduce.py and the programs in research/. Run python3 reproduce.py
from the repository root to reproduce the reported computational results.
The mathematical arguments are in the accompanying article.
"""

from itertools import product


def normalize(poly, modulus):
    return {exponent: value % modulus for exponent, value in poly.items()
            if value % modulus}


def multiply(a, b, modulus):
    out = {}
    for ex, value in a.items():
        for ey, other in b.items():
            pos = tuple(x + y for x, y in zip(ex, ey))
            out[pos] = (out.get(pos, 0) + value * other) % modulus
    return normalize(out, modulus)


def power(poly, exponent, modulus, dimension=2):
    result = {(0,) * dimension: 1}
    while exponent:
        if exponent & 1:
            result = multiply(result, poly, modulus)
        exponent >>= 1
        if exponent:
            poly = multiply(poly, poly, modulus)
    return result


def freeze(poly):
    return tuple(sorted(poly.items()))


def rational_automaton(prime, exponent):
    modulus = prime**exponent
    E = prime**(exponent - 1)
    Q = normalize({(0, 0): 1, (1, 0): -1, (0, 1): -1, (1, 1): -1}, modulus)
    digits = list(product(range(prime), repeat=2))
    factor = power(Q, (prime - 1) * E, modulus)
    # Independently check the denominator identity on which the construction rests.
    lhs = power(Q, prime * E, modulus)
    substituted = {(prime * i, prime * j): v for (i, j), v in Q.items()}
    assert lhs == power(substituted, E, modulus)
    initial = freeze(power(Q, E - 1, modulus))
    states = [initial]
    identifiers = {initial: 0}
    transitions = []
    for state in states:
        expanded = multiply(dict(state), factor, modulus)
        buckets = {digit: {} for digit in digits}
        for (i, j), value in expanded.items():
            buckets[i % prime, j % prime][i // prime, j // prime] = value
        row = []
        for digit in digits:
            nxt = freeze(buckets[digit])
            assert all(i <= E and j <= E for (i, j), _ in nxt)
            if nxt not in identifiers:
                identifiers[nxt] = len(states)
                states.append(nxt)
            row.append(identifiers[nxt])
        transitions.append(row)
    outputs = [dict(state).get((0, 0), 0) for state in states]
    # Verify each full polynomial identity from the stored state/transition table.
    for state, row in zip(states, transitions):
        reconstructed = {}
        for (r, s), nxt in zip(digits, row):
            for (i, j), value in states[nxt]:
                reconstructed[prime * i + r, prime * j + s] = value
        assert reconstructed == multiply(dict(state), factor, modulus)
    return {"base": prime, "modulus": modulus, "denominator_power": E,
            "digits": digits, "states": states, "transitions": transitions,
            "outputs": outputs, "initial": 0}


def minimize(machine):
    """Moore refinement; every source state is reachable by construction."""
    outputs = machine["outputs"]
    transitions = machine["transitions"]
    colors = outputs[:]
    rounds = 0
    while True:
        classes = {}
        new = []
        for i, row in enumerate(transitions):
            signature = (outputs[i], tuple(colors[j] for j in row))
            if signature not in classes:
                classes[signature] = len(classes)
            new.append(classes[signature])
        rounds += 1
        if len(set(new)) == len(set(colors)):
            return new, rounds
        colors = new


def reverse_to_msd(machine):
    """A state is the output function on the original LSD state set."""
    trans = machine["transitions"]
    vectors = [tuple(machine["outputs"])]
    ids = {vectors[0]: 0}
    transitions = []
    for vector in vectors:
        row = []
        for d in range(len(machine["digits"])):
            nxt = tuple(vector[trans[q][d]] for q in range(len(trans)))
            if nxt not in ids:
                ids[nxt] = len(vectors)
                vectors.append(nxt)
            row.append(ids[nxt])
        transitions.append(row)
    result = {"base": machine["base"], "modulus": machine["modulus"],
              "digits": machine["digits"], "states": vectors,
              "transitions": transitions,
              "outputs": [v[machine["initial"]] for v in vectors], "initial": 0}
    assert transitions[0][0] == 0, "The MSD start state must loop on zero padding."
    return result


def evaluate(machine, i, j, msd=False):
    p = machine["base"]
    digits = []
    while i or j:
        digits.append((i % p) * p + j % p)
        i //= p
        j //= p
    state = machine["initial"]
    for digit in reversed(digits) if msd else digits:
        state = machine["transitions"][state][digit]
    return machine["outputs"][state]


def recurrence_grid(side, modulus):
    grid = [[1] * side for _ in range(side)]
    for i in range(1, side):
        for j in range(1, side):
            grid[i][j] = (grid[i - 1][j] + grid[i][j - 1] + grid[i - 1][j - 1]) % modulus
    return grid


def block_from_vector(lsd, vector, side):
    """Values at side*(m,n)+(i,j) from the MSD state for (m,n)."""
    p = lsd["base"]
    k = 0
    scale = 1
    while scale < side:
        scale *= p
        k += 1
    assert scale == side
    block = []
    for i, j in product(range(side), repeat=2):
        q = lsd["initial"]
        x, y = i, j
        for _ in range(k):
            q = lsd["transitions"][q][(x % p) * p + y % p]
            x //= p
            y //= p
        block.append(vector[q])
    return tuple(block)


def certify_substitution(lsd, msd, side, project=False):
    p = lsd["base"]
    pairs = {}
    for vector in msd["states"]:
        small = block_from_vector(lsd, vector, side)
        large = block_from_vector(lsd, vector, p * side)
        if project:
            small = tuple(int(v != 0) for v in small)
            large = tuple(int(v != 0) for v in large)
        if small in pairs:
            assert pairs[small] == large, "Block output does not determine its substitution."
        pairs[small] = large
    blocks = list(pairs)
    ids = {block: i for i, block in enumerate(blocks)}
    rules = []
    for block in blocks:
        large = pairs[block]
        children = []
        for r, s in product(range(p), repeat=2):
            child = tuple(large[(r * side + i) * (p * side) + s * side + j]
                          for i, j in product(range(side), repeat=2))
            assert child in ids, "Substitution is not closed."
            children.append(ids[child])
        rules.append(children)
    start = block_from_vector(lsd, msd["states"][0], side)
    if project:
        start = tuple(int(v != 0) for v in start)
    start_id = ids[start]
    assert rules[start_id][0] == start_id
    return {"block_side": side, "factor": p, "projection": project,
            "blocks": blocks, "rules": rules, "initial": start_id}
