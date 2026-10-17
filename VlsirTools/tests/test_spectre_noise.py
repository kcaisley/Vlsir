"""Spectre noise analyses: netlisting, and reading their PSF ASCII results."""

import io
from pathlib import Path

import numpy as np
import pytest
import vlsir.spice_pb2 as vsp

from vlsirtools.netlist.spectre import SpectreNetlister
from vlsirtools.spice.spectre import Units, parse_psfascii, parse_psfascii_file


def _netlist_noise(an: vsp.NoiseInput) -> str:
    dest = io.StringIO()
    SpectreNetlister(dest=dest).write_noise(an)
    return dest.getvalue().strip()


def test_spectre_noise_netlist_single_ended_and_differential():
    single = vsp.NoiseInput(
        analysis_name="noise",
        output_p="out",
        input_source="vin",
        fstart=1e3,
        fstop=1e9,
        npts=10,
    )
    assert _netlist_noise(single) == (
        "noise (xtop.out 0) noise start=1000.0 stop=1000000000.0 dec=10 iprobe=xtop.vin"
    )
    differential = vsp.NoiseInput(
        analysis_name="nz",
        output_p="outp",
        output_n="outn",
        input_source="vinp",
        fstart=1.0,
        fstop=1e6,
        npts=5,
    )
    assert _netlist_noise(differential) == (
        "nz (xtop.outp xtop.outn) noise start=1.0 stop=1000000.0 dec=5 iprobe=xtop.vinp"
    )
    with pytest.raises(ValueError):
        _netlist_noise(
            vsp.NoiseInput(
                analysis_name="nz", output_p="out", fstart=1.0, fstop=1e6, npts=5
            )
        )


# Spectre's PSF ASCII for a resistor divider's noise (trimmed to two frequencies)
NOISE_PSF = """HEADER
"PSFversion" "1.00"
"simulator" "spectre"
"analysis type" "noise"
"analysis name" "noise"
"positive output signal" "out"
TYPE
"sweep" FLOAT DOUBLE PROP(
"key" "sweep"
)
"V/sqrt(Hz)" FLOAT DOUBLE PROP(
"units" "V/sqrt(Hz)"
"key" "noise"
)
"V/V" FLOAT DOUBLE PROP(
"units" "V/V"
"key" "noise"
)
"resistor" STRUCT(
"rn" FLOAT DOUBLE PROP(
"units" "V^2/Hz"
"noise_type" "thermal"
)
"total" FLOAT DOUBLE PROP(
"units" "V^2/Hz"
)
) PROP(
"key" "inst"
"master" "resistor"
)
SWEEP
"freq" "sweep" PROP(
"sweep_direction" 0
"units" "Hz"
)
TRACE
"r1" "resistor"
"out" "V/sqrt(Hz)"
"in" "V/sqrt(Hz)"
"gain" "V/V"
VALUE
"freq" 1.000000000000000e+03
"r1" (
2.652120789679248e-18
2.652120789679248e-18
)
"out" 2.574937275779378e-09
"in" 6.437343189468773e-09
"gain" 4.000000000000000e-01
"freq" 1.000000000000000e+06
"r1" (
2.000000000000000e-18
2.100000000000000e-18
)
"out" 2.000000000000000e-09
"in" 5.000000000000000e-09
"gain" 4.000000000000000e-01
END
"""

# An AC analysis (complex traces, one with a property block) and an operating point
AC_PSF = """HEADER
"analysis type" "ac"
"analysis name" "ac"
TYPE
"sweep" FLOAT DOUBLE PROP(
"key" "sweep"
)
"V" COMPLEX DOUBLE PROP(
"units" "V"
)
"I" COMPLEX DOUBLE PROP(
"units" "A"
)
SWEEP
"freq" "sweep" PROP(
"units" "Hz"
)
TRACE
"xtop.vvdd:p" "I" PROP(
"units" "A"
)
"xtop.out" "V"
VALUE
"freq" 1.000000000000000e+03
"xtop.vvdd:p" (-3.0e-05 -3.7e-09)
"xtop.out" (1.28e+01 -1.6e-04)
"freq" 1.000000000000000e+04
"xtop.vvdd:p" (-3.1e-05 -3.8e-08)
"xtop.out" (1.27e+01 -1.6e-03)
END
"""

OP_PSF = """HEADER
"analysis type" "dc"
"analysis name" "op"
TYPE
"V" FLOAT DOUBLE PROP(
"units" "V"
)
"I" FLOAT DOUBLE PROP(
"units" "A"
)
VALUE
"xtop.vvdd:p" "I" -4.746542429969111e-05 PROP(
"units" "A"
)
"xtop.out" "V" 6.609922997407638e-01
END
"""


def test_psfascii_noise_structs_and_totals():
    an = parse_psfascii_file(NOISE_PSF)
    assert (an.analysis_name, an.analysis_type) == ("noise", "noise")
    assert set(an.data) == {"freq", "r1.rn", "r1.total", "out", "in", "gain"}
    np.testing.assert_allclose(an.data["freq"], [1e3, 1e6])
    np.testing.assert_allclose(an.data["r1.total"], [2.652120789679248e-18, 2.1e-18])
    np.testing.assert_allclose(
        an.data["in"] * an.data["gain"], an.data["out"], rtol=1e-9
    )


def test_psfascii_ac_and_operating_point():
    ac = parse_psfascii_file(AC_PSF)
    np.testing.assert_allclose(ac.data["freq"], [1e3, 1e4])
    assert ac.data["xtop.out"].dtype == np.complex128
    np.testing.assert_allclose(ac.data["xtop.out"], [12.8 - 1.6e-4j, 12.7 - 1.6e-3j])
    np.testing.assert_allclose(
        ac.data["xtop.vvdd:p"], [-3.0e-5 - 3.7e-9j, -3.1e-5 - 3.8e-8j]
    )
    op = parse_psfascii_file(OP_PSF)
    np.testing.assert_allclose(op.data["xtop.vvdd:p"], [-4.746542429969111e-05])
    np.testing.assert_allclose(op.data["xtop.out"], [0.6609922997407638])


def test_psfascii_directory_skips_logs_and_info(tmp_path: Path):
    (tmp_path / "noise.noise").write_text(NOISE_PSF)
    (tmp_path / "ac.ac").write_text(AC_PSF)
    (tmp_path / "op.info").write_text(
        OP_PSF.replace('"dc"', '"info"').replace('"op"', '"opinfo"')
    )
    (tmp_path / "logFile").write_text("spectre log\n")
    assert sorted(parse_psfascii(tmp_path)) == ["ac", "noise"]


def test_nutbin_units_name_noise_quantities():
    for text in ("V^2/Hz", "V/sqrt(Hz)", "V/V", "A^2/Hz", "A/sqrt(Hz)"):
        assert Units.from_str(text).value == text
