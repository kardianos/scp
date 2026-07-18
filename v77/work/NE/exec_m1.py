#!/usr/bin/env python3
import runpy
import sys
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.argv = ["sandbox_m1_2d.py", "--quick"]
runpy.run_path("sandbox_m1_2d.py", run_name="__main__")
