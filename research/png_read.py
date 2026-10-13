"""Read the exact RGB pixels of the supplied noninterlaced 8-bit PNG files."""
from pathlib import Path
import struct
import zlib


def read_rgb(path):
    blob=Path(path).read_bytes()
    assert blob[:8]==b'\x89PNG\r\n\x1a\n'
    position=8
    compressed=bytearray()
    while position<len(blob):
        n=struct.unpack('>I',blob[position:position+4])[0]
        kind=blob[position+4:position+8]
        data=blob[position+8:position+8+n]
        assert zlib.crc32(kind+data)&0xffffffff==struct.unpack('>I',blob[position+8+n:position+12+n])[0]
        if kind==b'IHDR':
            width,height,depth,color,compression,filtering,interlace=struct.unpack('>IIBBBBB',data)
            assert (depth,color,compression,filtering,interlace)==(8,2,0,0,0)
        elif kind==b'IDAT': compressed.extend(data)
        position+=12+n
    raw=zlib.decompress(compressed)
    stride=3*width
    assert len(raw)==height*(stride+1)
    rows=[]
    previous=bytearray(stride)
    for y in range(height):
        kind=raw[y*(stride+1)]
        assert 0<=kind<=4
        row=bytearray(raw[y*(stride+1)+1:(y+1)*(stride+1)])
        for x in range(stride):
            a=row[x-3] if x>=3 else 0
            b=previous[x]
            c=previous[x-3] if x>=3 else 0
            if kind==0: correction=0
            elif kind==1: correction=a
            elif kind==2: correction=b
            elif kind==3: correction=(a+b)//2
            else:
                p=a+b-c
                distances=(abs(p-a),abs(p-b),abs(p-c))
                correction=(a,b,c)[distances.index(min(distances))]
            row[x]=(row[x]+correction)%256
        rows.append([tuple(row[x:x+3]) for x in range(0,stride,3)])
        previous=row
    return rows
