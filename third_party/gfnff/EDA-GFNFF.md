# EDA-GFNFF command-line program

`eda-gfnff` performs a molecular, interfragment decomposition using
the same geometry, topology, EEQ charges, repulsion parameters, and dispersion
parameters as the GFN-FF single-point calculation.

For fragments A and B, the reported quantity is

```text
E_NCI(A,B) = E_EEQ,pair(A,B) + E_rep,nonbonded(A,B) + E_disp(A,B)
           + E_HB(A,B) + E_XB(A,B)
```

The electrostatic, repulsion, and dispersion terms are accumulated from the corresponding
GFN-FF atom-pair expressions. The native GFN-FF hydrogen-bond and halogen-bond corrections
are accumulated separately from their three-center interaction lists. This is a
frozen-complex decomposition: it does not subtract separately recomputed monomer energies.
The CLI is intended for non-periodic, Cartesian molecular inputs whose fragments are not
covalently connected.

## Build

CMake:

```bash
cmake -S . -B build -Dbuild_exe=ON
cmake --build build
```

This builds both `gfnff` and `eda-gfnff` in the build directory.

Meson:

```bash
meson setup build -Dbuild_exe=true
ninja -C build
```

## XYZ input

The XYZ comment line can begin with Gaussian-style charge and multiplicity:

```text
0 1 optional comment
```

It can also contain keyed values:

```text
charge=-1 spin=2
chrg=0 multiplicity=1
```

Fragment IDs can be supplied as the fifth column of every atom line:

```text
6
0 1 water dimer
O   0.000000   0.000000   0.000000  1
H   0.758602   0.000000   0.504284  1
H  -0.758602   0.000000   0.504284  1
O   0.000000   0.000000   2.900000  2
H   0.758602   0.000000   3.404284  2
H  -0.758602   0.000000   3.404284  2
```

Run it with one integer net charge per fragment:

```bash
eda-gfnff water_dimer.xyz --frag-charges "0,0"
```

Alternatively, use a normal four-column XYZ and define each fragment with the same
repeatable selection syntax as `baneda`. Atom indices are one-based, ranges are
inclusive, and an optional name may precede `=`:

```bash
eda-gfnff water_dimer_plain.xyz \
  --frag water1=1-3 0 \
  --frag water2=4-6 0
```

The integer following each selection is that fragment's net charge. Instead of
attaching charges to the individual definitions, all fragment charges may be supplied
in option order with `--frag-charges`:

```bash
eda-gfnff water_dimer_plain.xyz \
  --frag water1=1-3 \
  --frag water2=4-6 \
  --frag-charges "0,0"
```

Every atom must belong to exactly one explicit fragment. Overlapping selections,
unassigned atoms, or a mixture in which only some `--frag` options carry charges are
rejected. The fragment charges must sum to the total charge read from the XYZ comment
line.

For compatibility with the former per-atom-array syntax, use `--frag-ids` with a
comma/space-separated list or a text file:

```bash
eda-gfnff water_dimer_plain.xyz \
  --frag-ids "1,1,1,2,2,2" \
  --frag-charges "0,0"
```

Fragment IDs from `--frag-ids` or XYZ column 5 are normalized in ascending numeric
order. Explicit repeated `--frag` definitions are numbered in command-line order.

## Gaussian GJF/COM input

The minimal Gaussian reader accepts Cartesian molecule specifications. Without CLI
fragment definitions, every atom must carry a `Fragment=N` field, and the
charge/multiplicity line must contain the total pair followed by one pair per fragment:

```text
%chk=water.chk
#p hf/3-21g

Water dimer

0 1  0 1  0 1
O(Fragment=1)   0.000000   0.000000   0.000000
H(Fragment=1)   0.758602   0.000000   0.504284
H(Fragment=1)  -0.758602   0.000000   0.504284
O(Fragment=2)   0.000000   0.000000   2.900000
H(Fragment=2)   0.758602   0.000000   3.404284
H(Fragment=2)  -0.758602   0.000000   3.404284
```

Repeated `--frag` definitions override Gaussian `Fragment=N` metadata and also allow a
plain Cartesian GJF/COM containing only the total charge/multiplicity pair. In that
mode, provide charges after every `--frag` or through `--frag-charges`.

Run directly:

```bash
eda-gfnff water_dimer.gjf
```

For an ion pair with fragment charges +1 and -1, the charge line is, for example:

```text
0 1  1 1  -1 1
```

The parser supports element symbols, atomic numbers, Gaussian atom metadata such as
`Element(Fragment=N)`, and an optional Gaussian freeze-code column before Cartesian
coordinates. Coordinates are interpreted in Angstrom.

The Gaussian parser is a minimal Fortran adaptation of the banelib INPUT/GJF parsing
approach, included with permission from the banelib author.

## Output

For every fragment pair the CLI prints electrostatic, nonbonded repulsion, dispersion,
GFN-FF hydrogen-bond correction, GFN-FF halogen-bond correction, and their total in the
selected output unit (`kcal/mol` by default or `kJ/mol` with `--unit`). It also prints
the three-term subtotal and the five-term all-pairs total in Hartree, kcal/mol, and
kJ/mol. The full GFN-FF single-point energy and gradient norm are reported separately.

The CLI also prints Multiwfn-style atomic contributions. For every cross-fragment
two-center electrostatic, repulsion, or dispersion interaction, one half of the pair
energy is assigned to each participating atom. Therefore, summing any atomic column
over the whole system exactly reproduces the corresponding all-fragment-pair energy.
The native GFN-FF HB/XB corrections are three-center interactions; their energies are
divided equally among the unique atoms participating in each correction.

An extended XYZ file is written automatically. By default, the input suffix is replaced
with `.eda.extxyz`; `-o` or `--output` selects another path. Coordinates are written in
Angstrom and all atomic energy properties are written in kcal/mol. The properties are:

```text
fragment
eda_electrostatic
eda_repulsion
eda_dispersion
eda_hbond
eda_xbond
eda_three_term
eda_total_nci
```

For example:

```bash
eda-gfnff water_dimer.xyz --frag-charges 0,0 -o water_atomic_eda.extxyz
```

The electrostatic entry is the cross-fragment pair part of the EEQ energy expression.
The one-center EEQ electronegativity/hardness terms are not assigned to fragment pairs.
The repulsion entry contains only the nonbonded GFN-FF repulsion loop. The dispersion
entry is the direct cross-fragment GFN-FF pair-dispersion contribution.

The H-bond and X-bond entries are the native GFN-FF three-center correction energies.
When an interaction contains atoms from exactly two fragments, its complete energy is
assigned to that fragment pair. If a three-center interaction spans three different
fragments, one third is assigned to each of the three fragment pairs so that the matrix
sum remains equal to the GFN-FF HB/XB contribution. Bonded, angle, torsion, bonded-ATM,
and other terms remain outside the reported interfragment NCI decomposition.

## CLI options

```text
-i, --input FILE
--frag [name=]SELECTION [charge]   (repeatable)
--frag-ids LIST_OR_FILE
--frag-charges LIST_OR_FILE
-o, --output FILE
--unit kcal/mol|kJ/mol
-T, --threads N
-v, --verbose
-q, --quiet
-h, --help
```

## Library opt-in

The new matrices are also available from the Fortran library. Set the exact partition
and fragment charges before initialization, and enable decomposition before the
single-point call:

```fortran
allocate(calc%userinput)
calc%userinput%fraglist = fragment_ids
calc%userinput%fragcharges = real(fragment_charges, real64)
calc%do_eda = .true.

call calc%init(...)
call calc%singlepoint(...)

! Upper-triangle matrices, i < j:
calc%res%eda%electrostatic(i,j)
calc%res%eda%repulsion(i,j)
calc%res%eda%dispersion(i,j)
calc%res%eda%hydrogen_bond(i,j)
calc%res%eda%halogen_bond(i,j)

! Atomic arrays; each sum equals the corresponding upper-triangle matrix sum:
calc%res%eda%atom_electrostatic(i)
calc%res%eda%atom_repulsion(i)
calc%res%eda%atom_dispersion(i)
calc%res%eda%atom_hydrogen_bond(i)
calc%res%eda%atom_halogen_bond(i)
```

EDA evaluation is opt-in, so ordinary GFN-FF library calls retain their original cost
and memory behavior.
