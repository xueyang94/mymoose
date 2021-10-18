#!/usr/bin/env python3

import sys
import os
import subprocess

if len(sys.argv) < 2:
    print("Usage:\nhit_wrapper.py gold_output [hit arguments]\n")

MOOSE_DIR = os.getenv('MOOSE_DIR', os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', '..', '..')))
hit = os.path.join(MOOSE_DIR, 'framework', 'contrib', 'hit', 'hit')
if not os.path.exists(hit):
    print('Failed to locate hit executable.')
    sys.exit(1)

gold_file = sys.argv[1]
hit_args = sys.argv[2:]

with open(gold_file, 'rb') as f:
    gold = f.read()

result = subprocess.run([hit] + hit_args, capture_output=True)
out = result.stdout

if out == gold:
    sys.exit(0)

print("Output:\n", out, "\ndoes not match gold:\n", gold)
sys.exit(1)
