#!/usr/bin/env python3
"""Certify automata, block rules and pixel agreement for the original figures."""
from pathlib import Path
import json
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'examples'))
from compute_examples import reverse_to_msd,minimize,certify_substitution,block_from_vector
from automatic_recurrences import OriginalNegativeArray,OriginalLabyrinthArray,learn_kernel,verify_recurrence
from png_read import read_rgb
OUT=Path(__file__).with_name('generated')


def save(name,value):
    (OUT/(name+'.json')).write_text(json.dumps(value,separators=(',',':'))+'\n')


def certify(prefix,figure,provider,predecessors,expected_numeric_rules):
    path=OUT/(prefix+'_machine.json')
    if path.exists():
        m=json.loads(path.read_text())
    else:
        history=[]
        for window in (1,2,3,4,6,8,12,16):
            m=learn_kernel(provider,2,2,window)
            verdict=verify_recurrence(m,predecessors,[1]*len(predecessors),2,forcing=1)
            history.append(dict(window=window,states=len(m['outputs']),verification=verdict))
            if verdict['valid']:break
        else:raise RuntimeError('No verified machine within the explicit search windows')
        save(prefix+'_discovery_history',history)
    m['modulus']=2
    verdict=verify_recurrence(m,predecessors,[1]*len(predecessors),2,forcing=1)
    assert verdict['valid']
    save(prefix+'_machine',m);save(prefix+'_certificate',verdict)
    msd=reverse_to_msd(m)
    substitution=certify_substitution(m,msd,16)
    assert len(substitution['blocks'])==expected_numeric_rules
    save(prefix+'_msd',msd);save(prefix+'_blocks16',substitution)
    block_ids={tuple(b):i for i,b in enumerate(substitution['blocks'])}
    state_blocks=[block_ids[block_from_vector(m,v,16)] for v in msd['states']]
    pixels=read_rgb(ROOT/'original'/(figure+'.png'))
    colors=read_rgb(ROOT/'original'/(figure+'Col.png'))
    assert len(pixels)==len(colors)==512
    assert all(len(row)==512 for row in pixels+colors)
    assert {rgb for row in pixels for rgb in row}<={(0,0,0),(255,255,255)}
    provider.ensure(512)
    assert all(provider((x,y))==int(pixels[y][x]==(0,0,0)) for y in range(512) for x in range(512))
    state_colors={}
    for y in range(512):
        for x in range(512):
            q=msd['initial']
            for k in reversed(range(9)):
                q=msd['transitions'][q][((x>>k)&1)*2+((y>>k)&1)]
            assert msd['outputs'][q]==int(pixels[y][x]==(0,0,0))
            if q in state_colors:assert state_colors[q]==colors[y][x]
            state_colors[q]=colors[y][x]
    assert len(state_colors)==len(msd['states'])
    colored=dict(msd,outputs=[state_colors[q] for q in range(len(msd['states']))])
    classes,rounds=minimize(colored)
    number=max(classes)+1
    groups=[[] for _ in range(number)]
    for q,c in enumerate(classes):groups[c].append(q)
    block_by_class=[];rules=[];rgb=[]
    for group in groups:
        assert len({state_blocks[q] for q in group})==1
        assert len({tuple(classes[r] for r in msd['transitions'][q]) for q in group})==1
        q=group[0]
        block_by_class.append(state_blocks[q])
        rgb.append(state_colors[q])
        rules.append([classes[r] for r in msd['transitions'][q]])
    initial=classes[msd['initial']]
    assert rules[initial][0]==initial
    # Every color-state rule decodes to the certified numerical block rule.
    for c,row in enumerate(rules):
        assert [block_by_class[j] for j in row]==substitution['rules'][block_by_class[c]]
    save(prefix+'_colored_substitution',dict(initial=initial,rules=rules,rgb=rgb,
         numerical_block_ids=block_by_class,msd_state_classes=classes))
    result=dict(predecessors=predecessors,forcing=1,boundary_value=1,
        minimal_lsd_states=len(set(minimize(m)[0])),minimal_msd_states=len(set(minimize(msd)[0])),
        numerical_block_rules=len(substitution['blocks']),colored_substitution_states=number,
        distinct_rgb_colors=len(set(state_colors.values())),colored_minimization_rounds=rounds,
        verified_numeric_pixels=512**2,verified_colored_pixels=512**2,
        all_msd_states_represented_in_color_image=True,recurrence_certificate=verdict)
    print(prefix,result,flush=True)
    return result


def main():
    OUT.mkdir(exist_ok=True)
    negative=certify('negative_forced','figExc',OriginalNegativeArray(forcing=1),
                     [(1,0),(1,2),(2,-2)],44)
    labyrinth=certify('labyrinth','figLab',OriginalLabyrinthArray(),
                     [(0,4),(1,3),(3,1),(4,0),(4,4)],204)
    literal=learn_kernel(OriginalNegativeArray(),2,2,1)
    literal_proof=verify_recurrence(literal,[(1,0),(1,2),(2,-2)],[1,1,1],2,forcing=0)
    assert literal_proof['valid'] and literal['outputs']==[1]
    report=dict(negative=negative,labyrinth=labyrinth,
                literal_negative_caption=dict(solution='constant 1',certificate=literal_proof),
                correction='The original source caption omits +1. The paper includes the corrected rule, which matches every supplied numerical pixel and admits the certified substitution.',
                image_coordinates='first recurrence coordinate horizontal, second vertical; black is 1')
    save('original_examples_report',report)


if __name__=='__main__':main()
