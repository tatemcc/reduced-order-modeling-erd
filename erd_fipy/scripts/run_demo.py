"""
CLI driver. Run with:
    python -m erd_fipy.run_demo
"""
import sys
try:
    from ..erd_fipy.stepping import run_sim
except Exception as e:
    print("Could not import or run FiPy-based simulator. Ensure FiPy is installed.")
    print("Error:", e)
    sys.exit(1)

if __name__ == "__main__":
    outpath = run_sim(seed=0)
    print("Wrote:", outpath)
