"""Spectre Monte Carlo analyses: the statement, its options and its child analyses."""

import io

import vlsir
import vlsir.spice_pb2 as vsp

from vlsirtools.netlist.spectre import SpectreNetlister


def _netlist_monte(an: vsp.MonteInput) -> str:
    dest = io.StringIO()
    SpectreNetlister(dest=dest).write_monte(an)
    return dest.getvalue().strip()


def test_spectre_monte_carlo_options_go_on_its_statement():
    an = vsp.MonteInput(
        analysis_name="mc",
        npts=400,
        seed=12345,
        an=[vsp.Analysis(op=vsp.OpInput(analysis_name="op"))],
        ctrls=[vsp.Control(param=vlsir.Param(name="savefamilyplots", value=vlsir.ParamValue(string_value="yes")))],
    )
    lines = _netlist_monte(an).splitlines()
    assert lines[0] == "mc montecarlo numruns=400 seed=12345 variations=all savefamilyplots=\"yes\" {"
    assert lines[1].strip().startswith("op dc") and lines[-1] == "}"
