"""Producer occupation/density conventions with a fake, already-solved object."""
from copy import deepcopy

import numpy as np
import pytest

from nanodesign import electronic_state as state
from e4_fixtures import build_snapshot
from e4_reference_cases import restricted_case, unrestricted_case


class FakeMolecule:
    def __init__(self, snapshot, case):
        self._snapshot = snapshot
        self._basis = deepcopy(snapshot["basis"]["expanded_basis"])
        self.ecp, self.pseudo = {}, {}
        self.cart, self.natm, self.nbas = False, 1, 3
        self.charge = snapshot["geometry"]["charge"]
        self.spin = snapshot["geometry"]["spin"]
        self.nelec = case["nelec"]

    def nao_nr(self):
        return 3

    def atom_symbol(self, index):
        return self._snapshot["geometry"]["symbols"][index]

    def atom_charges(self):
        return np.array(self._snapshot["geometry"]["atomic_numbers"])

    def atom_coords(self, unit):
        assert unit == "Bohr"
        return np.array(self._snapshot["geometry"]["positions_bohr"])

    def ao_labels(self, fmt):
        assert fmt is False
        return deepcopy(self._snapshot["basis"]["ao_labels"])

    def bas_atom(self, index):
        return 0

    def bas_angular(self, index):
        return 0

    def bas_kappa(self, index):
        return 0

    def bas_exp(self, index):
        return np.array(self._snapshot["basis"]["shells"][index]["exponents"])

    def bas_ctr_coeff(self, index):
        return np.array([[1.0]])

    def _libcint_ctr_coeff(self, index):
        return np.array([[1.0]])


class FakeMeanField:
    def __init__(self, reference, case):
        self.reference = reference
        snapshot = build_snapshot(case, reference=reference)
        self.mol = FakeMolecule(snapshot, case)
        self.settings = snapshot["settings"]
        self.converged, self.e_tot = True, -1.0
        self.mo_coeff, self.mo_occ = case["mo_coeff"].copy(), case["mo_occ"].copy()
        self.overlap = case["overlap"].copy()

    def istype(self, reference):
        return self.reference == reference

    def get_ovlp(self):
        return self.overlap.copy()

    def kernel(self, *args, **kwargs):
        raise AssertionError("E4 must never call a solver")


@pytest.mark.parametrize("reference,case", [("RKS", restricted_case()), ("UKS", unrestricted_case())])
def test_capture_preserves_correct_spin_electron_counts_without_solver(reference, case):
    mf = FakeMeanField(reference, case)
    captured = state.capture_snapshot(mf, settings=mf.settings, call_id="e4-capture", coordinate_frame_id="e4")
    assert captured["phase"] == "converged_scf_only"
    diagnostics = state.validate_snapshot(captured)
    for name, expected in zip(("alpha", "beta"), case["nelec"]):
        assert diagnostics[name]["density_electron_trace"] == pytest.approx(expected)
        assert captured["channels"][name]["n_electrons"] == expected
    if reference == "RKS":
        assert captured["source_occupations"]["spatial"] == [2.0, 0.0, 0.0]
        assert captured["channels"]["alpha"] == captured["channels"]["beta"]
    # Capture must own its evidence instead of aliasing mutable solver/settings data.
    mf.mo_coeff[:] = 0
    mf.settings["xc"] = "changed-after-capture"
    state.validate_snapshot(captured)
    assert captured["settings"]["xc"] == "pbe0"
    assert captured["electronic_state_identity_verified"] is False


def test_single_electron_empty_beta_channel_is_well_defined():
    case = unrestricted_case()
    case["nelec"] = (1, 0)
    case["mo_occ"] = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    mf = FakeMeanField("UKS", case)
    captured = state.capture_snapshot(mf, settings=mf.settings, call_id="e4-one-electron")
    report = state.compare_snapshots(captured, deepcopy(captured))
    assert report["comparison_status"] == "compared"
    assert report["metrics"]["beta"]["n_occupied"] == 0
    assert report["metrics"]["beta"]["singular_values"] == []
    assert report["metrics"]["beta"]["density_distance_squared"] == 0.0


@pytest.mark.parametrize("unsupported", ["complex", "fractional"])
def test_capture_rejects_unsupported_wavefunction_conventions(unsupported):
    mf = FakeMeanField("RKS", restricted_case())
    if unsupported == "complex":
        mf.mo_coeff = mf.mo_coeff.astype(complex)
        mf.mo_coeff[0, 0] += 1e-12j
    else:
        mf.mo_occ = np.array([1.5, 0.5, 0.0])
    with pytest.raises(state.ElectronicStateError):
        state.capture_snapshot(mf, settings=mf.settings, call_id="e4-unsupported")
