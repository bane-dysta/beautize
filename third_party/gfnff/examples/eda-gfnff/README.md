# EDA-GFNFF examples

After building with `-Dbuild_exe=ON`:

```bash
../../build/eda-gfnff water_dimer.xyz --frag-charges 0,0
../../build/eda-gfnff water_dimer.gjf
../../build/eda-gfnff iodomethane_water.xyz --frag-charges 0,0
```

Each run writes `INPUT.eda.extxyz` by default. The file contains the fragment ID
and the electrostatic, repulsion, dispersion, H-bond, X-bond, three-term, and total
NCI atomic contributions in kcal/mol by default. Add `--unit kJ/mol` for kJ/mol.
Use `-o FILE` to change the output path.

Both inputs describe the same geometry, fragment partition, and neutral fragment
charges, and therefore produce the same EDA-GFNFF values.

The iodomethane-water example exercises the separately reported GFN-FF halogen-bond
correction.
