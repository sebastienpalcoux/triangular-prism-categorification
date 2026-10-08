#!/usr/bin/env python3
"""Verify separately regenerated certificates without altering the frozen files."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import verify_f210

def main():
    if not __debug__:raise RuntimeError('Run without -O')
    data=json.loads((ROOT/'data/f210-equations.json').read_text())
    equations=verify_f210.verify_elimination(data)
    basis=verify_f210.verify_membership(data,equations,ROOT/'generated-certificates')
    result=verify_f210.verify_unit(basis,data['prime'])
    print('Regenerated certificate: PASS',result)
if __name__=='__main__':main()
