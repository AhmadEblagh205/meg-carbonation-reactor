"""Backward-compatible entry point for the old `runaway.py` command.

The analysis was renamed because a steady-state PFR parameter sweep is a
thermal-sensitivity screen, not a dynamic runaway or stability calculation.
"""
from thermal_sensitivity import main

if __name__ == "__main__":
    print("Note: `runaway.py` now runs the thermal-sensitivity screen; it does not claim a dynamic runaway boundary.")
    main()
