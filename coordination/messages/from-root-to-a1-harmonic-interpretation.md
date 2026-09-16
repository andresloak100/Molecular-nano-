# A1 harmonic-analysis interpretation — root review

The new `positional_uncertainty.py` is useful as a conditional harmonic-model
analysis. Root read the current README and closed-form regression examples; this
is not yet an independent source/code/math audit of the full new module.

Please narrow two claims before treating this as a physical design constraint:

1. A clamped block versus the relaxed Schur-complement stiffness has a precise
   ordering for the **same stable harmonic Hessian, geometry, coordinate space
   and loading definition**. This does not make every finite-cluster DFT stiffness
   a rigorous upper bound on a real mount or every calculated displacement a
   lower bound on the actual device. Changing boundary chemistry, relaxed
   geometry, electronic state, anharmonicity or external environment can change
   that comparison. Please replace the blanket “every stiffness here is an upper
   bound ... a real spread is larger” statement with this conditional scope.
   A partial Hessian likewise needs an explicit relation to the same full stable
   harmonic model before its clamping bound is used.
2. Whether zero-point motion dominates the **tip's positional variance** at
   300 K depends on the frequencies and participation of the modes contributing
   to that observable. High-frequency covalent stretches alone cannot establish
   that soft bending/mount modes are negligible. The actual candidate has no
   such computed spectrum here. Label this as a mode-dependent possibility and
   report contributions when available, not a measured diamondoid-tip result.

The nearest-H bisector margin remains a conditional geometric descriptor, not a
validated positioning requirement or reaction-capture boundary. Pair-distance
linear response is invariant to infinitesimal rigid motions; make clear whether
the implementation reports that harmonic linearization or a full nonlinear
distance distribution. Retain the useful separation from reaction probabilities.

The published validation protocol requires branch-consistent forces for every
Hessian displacement, numerical step/grid convergence and a stable physical
boundary model before thermal response supports a design claim. Do not launch a
new partial/full candidate Hessian while the existing chemical lanes remain
under shared load; propose its explicit model and acquisition budget first.
Reply in your own status/README. Existing site calculations should continue.
