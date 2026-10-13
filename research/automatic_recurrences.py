"""Finite-kernel discovery with exhaustive recurrence verification.

Discovery proposes machines from finite samples. Verification proves the
recurrence and boundary identities for every index by finite reachability.
"""
from itertools import product


def machine(base, dimension, transitions, outputs, initial=0):
    return dict(base=base,dimension=dimension,digits=list(product(range(base),repeat=dimension)),
                transitions=transitions,outputs=outputs,initial=initial)


def evaluate(m, point, start=None):
    q=m['initial'] if start is None else start
    x=list(point)
    assert min(x)>=0
    digit_ids={tuple(r):i for i,r in enumerate(m['digits'])}
    while any(x):
        digit=tuple(v%m['base'] for v in x)
        q=m['transitions'][q][digit_ids[digit]]
        x=[v//m['base'] for v in x]
    return m['outputs'][q]


def shifted(m, offset, outside=0):
    """Compute x -> m(x+offset), extended by outside at negative indices."""
    base=m['base']
    digits=m['digits']
    digit_ids={tuple(r):i for i,r in enumerate(digits)}
    states=[(m['initial'],tuple(offset))]
    ids={states[0]:0}
    transitions=[]
    outputs=[]
    for q,carry in states:
        outputs.append(evaluate(m,carry,q) if min(carry)>=0 else outside)
        row=[]
        for digit in digits:
            total=tuple(a+b for a,b in zip(digit,carry))
            residue=tuple(v%base for v in total)
            nxt=(m['transitions'][q][digit_ids[residue]],tuple(v//base for v in total))
            assert all(min(c,0)<=v<=max(c,1) for v,c in zip(nxt[1],offset))
            if nxt not in ids:
                ids[nxt]=len(states)
                states.append(nxt)
            row.append(ids[nxt])
        transitions.append(row)
    return machine(base,m['dimension'],transitions,outputs)


def interior_machine(base, widths):
    digits=list(product(range(base),repeat=len(widths)))
    cap=max(1,max(widths))
    states=[(tuple(0 for _ in widths),1)]
    ids={states[0]:0}
    transitions=[]
    outputs=[]
    for values,scale in states:
        outputs.append(int(all(v>=d for v,d in zip(values,widths))))
        row=[]
        for digit in digits:
            nxt=(tuple(min(d,v+scale*r) for v,r,d in zip(values,digit,widths)),min(cap,base*scale))
            if nxt not in ids:
                ids[nxt]=len(states)
                states.append(nxt)
            row.append(ids[nxt])
        transitions.append(row)
    return machine(base,len(widths),transitions,outputs)


def constant_machine(base,dimension,value):
    return machine(base,dimension,[[0]*(base**dimension)],[value])


def verify_recurrence(candidate,predecessors,coefficients,modulus,boundary=1,forcing=0,max_states=None):
    """Prove an additive recurrence everywhere, or return an exact witness.

    boundary and forcing may be machines or constant values. Boundary data
    are interpreted only on the prescribed faces.
    """
    base=candidate['base']
    dim=candidate['dimension']
    digits=[tuple(x) for x in candidate['digits']]
    for q,row in enumerate(candidate['transitions']):
        assert len(row)==len(digits)
        assert all(0<=v<len(candidate['outputs']) for v in row)
        assert candidate['outputs'][q]==candidate['outputs'][row[0]], 'Zero padding changes output'
    widths=tuple(max(0,*(v[i] for v in predecessors)) for i in range(dim))
    machines=[candidate]+[shifted(candidate,tuple(-x for x in v)) for v in predecessors]
    inside_index=len(machines)
    machines.append(interior_machine(base,widths))
    for data in (boundary,forcing):
        machines.append(data if isinstance(data,dict) else constant_machine(base,dim,data))
    assert all([tuple(x) for x in m['digits']]==digits for m in machines)
    initial=tuple(m['initial'] for m in machines)
    states=[initial]
    ids={initial:0}
    parents=[None]
    for state_index,state in enumerate(states):
        outputs=[m['outputs'][q] for m,q in zip(machines,state)]
        expected=(sum(c*v for c,v in zip(coefficients,outputs[1:inside_index]))+outputs[-1])%modulus if outputs[inside_index] else outputs[-2]
        if outputs[0]!=expected:
            word=[]
            at=state_index
            while parents[at] is not None:
                at,d=parents[at]
                word.append(d)
            word.reverse()
            point=tuple(sum((base**j)*digits[d][i] for j,d in enumerate(word)) for i in range(dim))
            return {'valid':False,'witness':point,'actual':outputs[0],'expected':expected,
                    'interior':bool(outputs[inside_index]),'reachable_before_failure':len(states)}
        for d in range(len(digits)):
            nxt=tuple(m['transitions'][q][d] for m,q in zip(machines,state))
            if nxt not in ids:
                ids[nxt]=len(states)
                states.append(nxt)
                parents.append((state_index,d))
                if max_states is not None and len(states)>max_states:
                    raise RuntimeError(f'Verification exceeded its explicit operational limit of {max_states} states')
    return {'valid':True,'reachable_product_states':len(states),
            'component_state_counts':[len(m['outputs']) for m in machines]}


def learn_kernel(value,base,dimension,window):
    digits=list(product(range(base),repeat=dimension))
    sample=list(product(range(window),repeat=dimension))
    def fingerprint(scale,residue):
        return tuple(value(tuple(scale*v+r for v,r in zip(x,residue))) for x in sample)
    states=[(1,(0,)*dimension)]
    ids={fingerprint(*states[0]):0}
    transitions=[]
    outputs=[]
    for scale,residue in states:
        outputs.append(value(residue))
        row=[]
        for digit in digits:
            nxt=(base*scale,tuple(r+scale*d for r,d in zip(residue,digit)))
            key=fingerprint(*nxt)
            if key not in ids:
                ids[key]=len(states)
                states.append(nxt)
            row.append(ids[key])
        transitions.append(row)
    result=machine(base,dimension,transitions,outputs)
    result['discovery']={'window':window,'section_representatives':states}
    return result


class OriginalNegativeArray:
    """Exact bits for P={(1,0),(1,2),(2,-2)}, with two unit boundary rows/columns."""
    def __init__(self,maximum_side=32768,forcing=0):
        self.rows=[]
        self.side=0
        self.maximum_side=maximum_side
        self.forcing=forcing
    def ensure(self,side):
        if side<=self.side: return
        size=1
        while size<side: size*=2
        if size>self.maximum_side:
            raise RuntimeError(f'Sampling needs side {size}, exceeding the explicit limit {self.maximum_side}')
        top=2*size-2
        rows=[(1<<(top+1))-1]
        if size>1: rows.append((1<<top)-1)
        for i in range(2,size):
            row=rows[-1]^(rows[-1]<<2)^(rows[-2]>>2)
            mask=(1<<(top-i+1))-1
            if self.forcing: row^=mask
            row=(row|3)&mask
            rows.append(row)
        self.rows=rows
        self.side=size
    def __call__(self,point):
        m,n=point
        self.ensure(max(m,n)+1)
        return (self.rows[m]>>n)&1


class OriginalLabyrinthArray:
    """Exact bits for the five-input nonnegative recurrence in the source."""
    def __init__(self,maximum_side=32768):
        self.rows=[]
        self.side=0
        self.maximum_side=maximum_side
    def ensure(self,side):
        if side<=self.side: return
        size=1
        while size<side: size*=2
        if size>self.maximum_side:
            raise RuntimeError(f'Sampling needs side {size}, exceeding the explicit limit {self.maximum_side}')
        mask=(1<<size)-1
        rows=[mask]*min(4,size)
        for i in range(4,size):
            row=(rows[i-1]<<3)^(rows[i-3]<<1)^rows[i-4]^(rows[i-4]<<4)^mask
            row=(row|15)&mask
            step=4
            while step<size:
                row^=(row<<step)&mask
                step*=2
            rows.append(row)
        self.rows=rows
        self.side=size
    def __call__(self,point):
        m,n=point
        self.ensure(max(m,n)+1)
        return (self.rows[m]>>n)&1
