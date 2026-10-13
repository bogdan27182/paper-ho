#!/usr/bin/env python3
"""Verify a finite polynomial certificate and the all-length modulo-27 count."""
from pathlib import Path
from collections import deque
import json
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'examples'))
from compute_examples import normalize,multiply,power,minimize,recurrence_grid,evaluate,rational_automaton
OUT=Path(__file__).with_name('generated')


def main():
    OUT.mkdir(exist_ok=True)
    path=OUT/'mod27_rational_machine.json'
    if not path.exists():
        path.write_text(json.dumps(rational_automaton(3,3),separators=(',',':'))+'\n')
    m=json.loads(path.read_text())
    assert (m['base'],m['modulus'],m['denominator_power'],m['initial'])==(3,27,9,0)
    assert m['digits']==[[i,j] for i in range(3) for j in range(3)]
    polys=[{tuple(ex):value for ex,value in state} for state in m['states']]
    assert len(polys)==len(m['transitions'])==len(m['outputs'])==1362
    assert len({tuple(sorted(p.items())) for p in polys})==len(polys)
    Q=normalize({(0,0):1,(1,0):-1,(0,1):-1,(1,1):-1},27)
    factor=power(Q,18,27)
    assert polys[0]==power(Q,8,27)
    assert power(Q,27,27)==power({(3*i,3*j):v for (i,j),v in Q.items()},9,27)
    for q,poly in enumerate(polys):
        assert all(0<=i<=9 and 0<=j<=9 and 0<v<27 for (i,j),v in poly.items())
        assert m['outputs'][q]==poly.get((0,0),0)
        assert len(m['transitions'][q])==9
        reconstructed={}
        for (r,s),nxt in zip(m['digits'],m['transitions'][q]):
            assert 0<=nxt<len(polys)
            for (i,j),v in polys[nxt].items():
                at=(3*i+r,3*j+s)
                assert at not in reconstructed
                reconstructed[at]=v
        assert reconstructed==multiply(poly,factor,27),q
    reachable={0};pending=deque([0])
    while pending:
        for nxt in m['transitions'][pending.popleft()]:
            if nxt not in reachable:reachable.add(nxt);pending.append(nxt)
    assert len(reachable)==len(polys)
    classes,rounds=minimize(m)
    assert len(set(classes))==1362
    def step(v):return [sum(v[nxt] for nxt in row) for row in m['transitions']]
    vectors=[[int(v!=0) for v in m['outputs']]]
    for _ in range(7):vectors.append(step(vectors[-1]))
    residual=[a-24*b+192*c-512*d for a,b,c,d in zip(vectors[7],vectors[6],vectors[5],vectors[4])]
    assert not any(residual)
    initial=[v[0] for v in vectors]
    assert initial[4:7]==[6165,54384,477888]
    for r in range(4,8):assert initial[r]==8**(r-4)*(18*r*r+471*r+3993)
    # Direct recurrence arithmetic uses neither the automaton nor transfer matrices.
    grid=recurrence_grid(243,27)
    for i,row in enumerate(grid):
        for j,v in enumerate(row):assert evaluate(m,i,j)==v
    for r in range(6):
        side=3**r
        assert sum(grid[i][j]!=0 for i in range(side) for j in range(side))==initial[r]
    certificate=dict(states=1362,transitions=1362*9,minimal_lsd_states=1362,minimization_rounds=rounds,
        vectors_M_power_v=vectors,identity='M^4 (M-8I)^3 v = 0',
        initial_nonzero_counts=initial,formula='8^(r-4) * (18*r^2 + 471*r + 3993), r >= 4',
        independent_grid_side=243,independent_grid_values=243**2,
        scope='The finite polynomial identities certify the whole infinite array. The integer vector identity and three initial counts prove the formula for all r >= 4.')
    (OUT/'mod27_count_certificate.json').write_text(json.dumps(certificate,separators=(',',':'))+'\n')
    print('Modulo 27: 1362 polynomial states, 12258 transition identities, M^4(M-8I)^3 v=0, and 59049 independent values verified.')


if __name__=='__main__':main()
