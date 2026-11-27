# Preserve the original local momentum interval under hard refinement

The original local daughter substitution uses two cells covering
the same physical momentum interval that has eight cells on
the fourfold hard mesh. Use the exact recorded continuous
projection in those cells and the original native projection
everywhere else. Match the cold seed and full hard edges within
each paired archive; verify the incoming raw native source.
This evaluates one shared local physical interval, not one
fixed number of cells that changes size with resolution.

The original initial local-minus-native cold rate is -0.01551056,
compared with -0.01167084 after hard refinement. All eight
initial cold-bin shifts remain negative; the total suppression
weakens by about 24.75%. This is a static discretization
sensitivity. It does not yet test the finite local rebound on
the more resolved hard mesh or imply another choice of physical
momentum interval behaves the same.

Use `early_pole_daughter_original_interval_archive_control.py`
with both compatible fine archive bundles.
