# Third-party sources and notices

## Bundled GFN-FF

`third_party/gfnff/` is derived from the user-supplied `ead-gfnff.tgz`, not a
freshly downloaded or unspecified upstream revision. Its original source
notices, license, README and EDA documentation are preserved. Its upstream
README identifies the standalone project as `pprcht/gfnff`, based on the
GFN-FF method of S. Spicher and S. Grimme and code adapted from `xtb`.
Please see that README and the individual source headers for full authorship,
method references and component-specific attribution.

Original uploaded archive SHA-256:

```text
b6dbfe67913179bcc7704b208b92e09c684d52e3857437a8119e1d353c1048f2  ead-gfnff.tgz
```

The GFN-FF sources and modifications remain **LGPL-3.0-or-later**, according
to the existing source notices. The LGPL version 3 text is at
`third_party/gfnff/LICENSE`; the incorporated GPL version 3 text is also
provided at `licenses/GPL-3.0.txt`.

The following changes were made for beautize and are supplied as reviewable
patches, applied in this order:

1. `patches/gfnff-explicit-topology.patch`: adds an optional authoritative
   molecular bond adjacency input; validates it before topology initialization;
   bypasses distance guessing for covalent edges; retains the explicit final
   neighbor list; disables restart for the external graph path.
2. `patches/gfnff-serial-hb-buffers.patch`: makes hydrogen-bond list buffer
   declarations/allocation/copy operations available in a serial, non-OpenMP
   build, while retaining the OpenMP parallel directives.

The vendored tree already contains both modifications. Applying these two
patches to the extracted original archive has been checked to reproduce the
vendored tree exactly. Existing bundled build auxiliaries and component files
retain their own notices; the root beautize MIT license does not relicense them.

## Independent beautize implementation

The original C++17 application, small Fortran adapter, native tests and
project documentation are provided under the root `LICENSE` (MIT). They
contain no copied banelib implementation and are not a transcription of ASE.
Calling the patched library does not make banelib, ASE or Python dependencies.

The prebuilt executable statically links the supplied GFN-FF library, and
uses dynamic system libraries for Fortran/C++/OpenMP and BLAS/LAPACK. The
prebuilt distribution includes the complete matching source archive, build
instructions, licenses and patches, so recipients can modify the library
and rebuild/relink the application. No proprietary-only object or build step
is required. System shared libraries are not bundled in the binary archive.

## Design reference archives (not bundled dependencies)

The ase-3.29.0 source was consulted for the constraint API behavior. 
**Not included or linked** into beautize and not needed to build or run it.



