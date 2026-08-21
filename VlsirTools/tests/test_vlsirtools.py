"""
# Vlsir Python Tools
Unit Tests
"""

import os
import numpy as np
import pytest
from io import StringIO
from typing import Dict, List, Optional

import vlsir.utils_pb2 as vutils
from vlsir.utils_pb2 import Reference, QualifiedName, ParamValue, Param
import vlsir.circuit_pb2 as vckt
from vlsir.circuit_pb2 import (
    Module,
    Signal,
    Connection,
    ConnectionTarget,
    Port,
    Instance,
    Package,
)
import vlsir.spice_pb2 as vsp

import vlsirtools
from vlsirtools.spice import (
    SupportedSimulators,
    SimOptions,
    ResultFormat,
    sim,
)
import vlsirtools.spice.sim_data as sd
from vlsirtools.spice.sim_data import AnalysisType

"""
## Utility & Helper Functions
"""


def _connections(**kwargs):
    """Create a list of `Connection`s from keyword args of the form `portname=conn_target`, where `conn_target` is a `ConnectionTarget`."""
    return [Connection(portname=key, target=value) for key, value in kwargs.items()]


def _params(**kwargs):
    """Create a list of `Param`s from keyword args of the form
    `r=ParamValue(double_value=1e3)`"""
    return [Param(name=key, value=value) for key, value in kwargs.items()]


def _prim(name: str) -> Reference:
    """Create a `Reference` to primitive `name`"""
    return Reference(external=QualifiedName(domain="vlsir.primitives", name=name))


def test_version():
    assert vlsirtools.__version__ == "7.0.0"  # VLSIR_VERSION


def test_verilog_netlist1():
    """Test netlisting to a handful of formats, including Verilog-compatible contents."""

    # "Verilog Compatibility" requires:
    # * All ports must be directed. No "NONE" directions.
    # * No primitive instances.
    #
    # The test here covers module, signal, and port creation,
    # and basic instantiation of other local Modules.

    def _ports() -> List[vckt.Port]:
        return [
            Port(direction="INPUT", signal="inp"),
            Port(direction="OUTPUT", signal="out"),
            Port(direction="INOUT", signal="io"),
        ]

    def _signals() -> List[vckt.Signal]:
        return [
            Signal(name="inp", width=1),
            Signal(name="out", width=1),
            Signal(name="io", width=1),
        ]

    pkg = Package(
        domain="vlsirtools.tests.test_verilog_netlist1",
        modules=[
            # Inner, content-less Module, with a port of each direction
            Module(
                name="inner",
                parameters=_params(
                    a=ParamValue(int64_value=3),
                    b=ParamValue(double_value=1e-9),
                    d=ParamValue(literal="1+1"),
                    e=ParamValue(
                        prefixed=vutils.Prefixed(prefix="MICRO", string_value="11.11")
                    ),
                ),
                ports=_ports(),
                signals=_signals(),
            ),
            # Outer top Module which instantiates `inner`
            Module(
                name="top",
                ports=[],
                signals=_signals(),
                instances=[
                    Instance(
                        name="inner",
                        module=Reference(local="inner"),
                        parameters=_params(
                            a=ParamValue(int64_value=4),
                            b=ParamValue(double_value=2e-9),
                            d=ParamValue(literal="2+2"),
                            e=ParamValue(
                                prefixed=vutils.Prefixed(
                                    prefix="MICRO", string_value="22.22"
                                )
                            ),
                        ),
                        connections=[
                            Connection(portname=name, target=ConnectionTarget(sig=name))
                            for name in ("inp", "out", "io")
                        ],
                    )
                ],
            ),
        ],
    )
    dest = StringIO()
    vlsirtools.netlist(pkg=pkg, dest=dest, fmt="verilog")
    # While verilog is the point here, the other formats should work too:
    vlsirtools.netlist(pkg=pkg, dest=dest, fmt="spice")
    vlsirtools.netlist(pkg=pkg, dest=dest, fmt="spectre")
    vlsirtools.netlist(pkg=pkg, dest=dest, fmt="xyce")


def test_spice_netlist1():
    """Test spice netlisting, including instantiating primitives"""

    def default_conns() -> Dict:
        # Shorthand for connections between node "vvv" and "VSS", used for many instances below.
        return _connections(
            p=ConnectionTarget(sig="vvv"),
            n=ConnectionTarget(sig="VSS"),
        )

    # Create ExternalModules for a few model-based device-types
    mos1 = vlsirtools.primitives.mos(name="mos1", domain="mymodels")
    d1 = vlsirtools.primitives.diode(name="d1", domain="mymodels")
    bjt1 = vlsirtools.primitives.bipolar(name="bjt1", domain="mymodels")
    t1 = vlsirtools.primitives.tline(name="t1", domain="mymodels")

    def _eref(emod: vckt.ExternalModule) -> Reference:
        return Reference(external=emod.name)

    pkg = Package(
        domain="vlsir.tests.test_netlist1",
        ext_modules=[mos1, d1, bjt1, t1],
        modules=[
            Module(
                name="mid",
                parameters=_params(
                    r=ParamValue(double_value=1e3),
                    l=ParamValue(double_value=1e-9),
                    c=ParamValue(double_value=1e-15),
                ),
                ports=[
                    Port(direction="NONE", signal="vvv"),
                    Port(direction="NONE", signal="VSS"),
                ],
                signals=[Signal(name="vvv", width=1), Signal(name="VSS", width=1)],
                instances=[
                    Instance(
                        name="r",
                        module=_prim("resistor"),
                        connections=default_conns(),
                        parameters=_params(r=ParamValue(double_value=1e3)),
                    ),
                    Instance(
                        name="l",
                        module=_prim("inductor"),
                        connections=default_conns(),
                        parameters=_params(l=ParamValue(double_value=1e-9)),
                    ),
                    Instance(
                        name="c",
                        module=_prim("capacitor"),
                        connections=default_conns(),
                        parameters=_params(c=ParamValue(double_value=1e-15)),
                    ),
                    Instance(
                        name="v",
                        module=_prim("vdc"),
                        connections=default_conns(),
                        parameters=_params(
                            dc=ParamValue(double_value=1.1),
                            ac=ParamValue(double_value=0),
                        ),
                    ),
                    Instance(
                        name="i",
                        module=_prim("isource"),
                        connections=default_conns(),
                        parameters=_params(dc=ParamValue(double_value=1e-6)),
                    ),
                    # Model-based instances
                    Instance(
                        name="m",
                        module=_eref(mos1),
                        connections=_connections(
                            d=ConnectionTarget(sig="VSS"),
                            g=ConnectionTarget(sig="VSS"),
                            s=ConnectionTarget(sig="VSS"),
                            b=ConnectionTarget(sig="VSS"),
                        ),
                        parameters=[],
                    ),
                    Instance(
                        name="q",
                        module=_eref(bjt1),
                        connections=_connections(
                            c=ConnectionTarget(sig="VSS"),
                            b=ConnectionTarget(sig="VSS"),
                            e=ConnectionTarget(sig="VSS"),
                        ),
                        parameters=[],
                    ),
                    Instance(
                        name="d",
                        module=_eref(d1),
                        connections=default_conns(),
                        parameters=[],
                    ),
                    Instance(
                        name="t",
                        module=_eref(t1),
                        connections=_connections(
                            p1p=ConnectionTarget(sig="VSS"),
                            p1n=ConnectionTarget(sig="VSS"),
                            p2p=ConnectionTarget(sig="VSS"),
                            p2n=ConnectionTarget(sig="VSS"),
                        ),
                        parameters=[],
                    ),
                ],
            ),
            Module(
                name="top",
                ports=[
                    Port(direction="NONE", signal="VSS"),
                ],
                signals=[Signal(name="vvv", width=1), Signal(name="VSS", width=1)],
                instances=[
                    Instance(
                        name="imid",
                        module=Reference(local="mid"),
                        connections=_connections(
                            vvv=ConnectionTarget(sig="vvv"),
                            VSS=ConnectionTarget(sig="VSS"),
                        ),
                    )
                ],
            ),
        ],
    )
    dest = StringIO()
    vlsirtools.netlist(pkg=pkg, dest=dest, fmt="spice")


def test_spectre_bit_voltage_source():
    """Netlist the VLSIR bit source using Spectre's native syntax."""

    pkg = Package(
        domain="vlsirtools.tests.test_spectre_bit_voltage_source",
        modules=[
            Module(
                name="bit_source_tb",
                signals=[Signal(name="out", width=1), Signal(name="VSS", width=1)],
                instances=[
                    Instance(
                        name="bits",
                        module=_prim("vbit"),
                        connections=_connections(
                            p=ConnectionTarget(sig="out"),
                            n=ConnectionTarget(sig="VSS"),
                        ),
                        parameters=_params(
                            data=ParamValue(literal="0011010"),
                            period=ParamValue(literal="1u"),
                            val0=ParamValue(int64_value=0),
                            val1=ParamValue(double_value=1.2),
                            delay=ParamValue(literal="1n"),
                            rise=ParamValue(literal="3n"),
                            fall=ParamValue(literal="3n"),
                            rptstart=ParamValue(int64_value=1),
                            rpttimes=ParamValue(int64_value=0),
                        ),
                    )
                ],
            )
        ],
    )

    netlist = StringIO()
    vlsirtools.netlist(pkg=pkg, dest=netlist, fmt="spectre")
    text = netlist.getvalue()

    assert "vsource" in text
    assert "type=bit" in text
    assert 'data="0011010"' in text
    assert "period=1u" in text
    assert "val0=0" in text
    assert "val1=1.2" in text
    assert "delay=1n" in text
    assert "rise=3n" in text
    assert "fall=3n" in text
    assert "rptstart=1" in text
    assert "rpttimes=0" in text

    with pytest.raises(ValueError, match="vbit"):
        vlsirtools.netlist(pkg=pkg, dest=StringIO(), fmt="spice")


def test_spectre_sine_voltage_source():
    """Translate SPICE SIN fields to Spectre's native sine parameters."""

    pkg = Package(
        domain="vlsirtools.tests.test_spectre_sine_voltage_source",
        modules=[
            Module(
                name="sine_source_tb",
                signals=[Signal(name="out", width=1), Signal(name="VSS", width=1)],
                instances=[
                    Instance(
                        name="sine",
                        module=_prim("vsin"),
                        connections=_connections(
                            p=ConnectionTarget(sig="out"),
                            n=ConnectionTarget(sig="VSS"),
                        ),
                        parameters=_params(
                            voff=ParamValue(double_value=0.6),
                            vamp=ParamValue(double_value=0.5),
                            freq=ParamValue(literal="10k"),
                            td=ParamValue(literal="2n"),
                            phase=ParamValue(double_value=90.0),
                        ),
                    )
                ],
            )
        ],
    )

    netlist = StringIO()
    vlsirtools.netlist(pkg=pkg, dest=netlist, fmt="spectre")
    text = netlist.getvalue()

    assert "vsource" in text
    assert "type=sine" in text
    assert "sinedc=0.6" in text
    assert "ampl=0.5" in text
    assert "freq=10k" in text
    assert "delay=2n" in text
    assert "sinephase=90" in text
    assert "voff=" not in text
    assert "vamp=" not in text
    assert "td=" not in text
    assert "phase=" not in text.replace("sinephase=", "")


def test_netlist_hdl21_ideal1():
    """Test that netlisting a (deprecated) `hdl21.ideal` element fails with a `RuntimeError`."""

    pkg = Package(
        domain="vlsir.tests.test_netlist_hdl21_ideal1",
        modules=[
            Module(
                name="mid",
                ports=[
                    Port(direction="NONE", signal="vvv"),
                    Port(direction="NONE", signal="VSS"),
                ],
                signals=[Signal(name="vvv", width=1), Signal(name="VSS", width=1)],
                instances=[
                    Instance(
                        name="r1",
                        module=Reference(
                            external=QualifiedName(
                                domain="hdl21.ideal",
                                name="IdealResistor",  # FIXME: being deprecated
                            )
                        ),
                        connections=[
                            Connection(
                                portname="p", target=ConnectionTarget(sig="vvv")
                            ),
                            Connection(
                                portname="n", target=ConnectionTarget(sig="VSS")
                            ),
                        ],
                    ),
                ],
            )
        ],
    )
    dest = StringIO()
    with pytest.raises(RuntimeError):
        vlsirtools.netlist(pkg=pkg, dest=dest, fmt="spice")


def dummy_testbench_package():
    """Create and return a `circuit.Package` with a single, (near) empty testbench module.
    Some simulators *really* don't like empty DUT content, and others don't like singly-connected nodes.
    So the simplest test-bench is two resistors, in parallel, between ground and a single "other node".
    """

    def _r(name: str) -> Instance:
        # Create a canned instance of `vlsir.primitives.resistor` with instance-name `name`
        return Instance(
            name=name,  # <= Instance name argument here
            module=Reference(
                external=QualifiedName(domain="vlsir.primitives", name="resistor")
            ),
            connections=[
                Connection(portname="p", target=ConnectionTarget(sig="the_other_node")),
                Connection(portname="n", target=ConnectionTarget(sig="VSS")),
            ],
            parameters=[Param(name="r", value=ParamValue(double_value=1e3))],
        )

    return Package(
        domain="vlsirtools.tests.dummy_testbench_package",
        modules=[
            Module(
                name="dummy_testbench",
                ports=[
                    Port(direction="NONE", signal="VSS"),
                ],
                signals=[
                    Signal(name="VSS", width=1),
                    Signal(name="the_other_node", width=1),
                ],
                instances=[
                    _r("r1"),
                    _r("r2"),
                ],
            )
        ],
    )


def test_netlist_dummy_testbench():
    """Test netlisting the empty testbench package, used later for simulation tests"""

    dest = StringIO()
    vlsirtools.netlist(pkg=dummy_testbench_package(), dest=dest, fmt="spice")
    vlsirtools.netlist(pkg=dummy_testbench_package(), dest=dest, fmt="spectre")
    vlsirtools.netlist(pkg=dummy_testbench_package(), dest=dest, fmt="xyce")

    # The testbench package is not verilog-compatible; check that it fails.
    with pytest.raises(RuntimeError):
        vlsirtools.netlist(pkg=dummy_testbench_package(), dest=dest, fmt="verilog")


def dummy_sim(skip: Optional[List[AnalysisType]] = None):
    """Create a dummy `SimInput` for the `dummy_testbench` package."""

    if skip is None:
        skip = []

    # Set up an `analysis` for each analysis-type to run on the empty testbench.
    # FIXME: add noise, sweep, etc.
    analyses_by_type: Dict[AnalysisType, vsp.Analysis] = {
        AnalysisType.OP: vsp.Analysis(
            op=vsp.OpInput(
                analysis_name="op1",
                ctrls=[],
            )
        ),
        AnalysisType.DC: vsp.Analysis(
            dc=vsp.DcInput(
                analysis_name="dc1",
                indep_name="DUMMY",
                sweep=vsp.Sweep(
                    linear=vsp.LinearSweep(
                        start=0.0,
                        stop=1.0,
                        step=0.1,
                    ),
                ),
                ctrls=[],
            )
        ),
        AnalysisType.TRAN: vsp.Analysis(
            tran=vsp.TranInput(
                analysis_name="tr1",
                tstop=1e-9,
                tstep=1e-12,
                ic={},
                ctrls=[],
            )
        ),
        AnalysisType.AC: vsp.Analysis(
            ac=vsp.AcInput(
                analysis_name="ac1",
                fstart=1e3,
                fstop=1e6,
                npts=10,
                ctrls=[],
            )
        ),
    }
    # Rejigger those into a list, skipping any that are in the `skip` list.
    an = [a for t, a in analyses_by_type.items() if t not in skip]

    return vsp.SimInput(
        top="dummy_testbench",
        pkg=dummy_testbench_package(),
        an=an,
        ctrls=[
            vsp.Control(param=Param(name="DUMMY", value=ParamValue(int64_value=0))),
            # This parameter "DUMMY" is the sweep-variable for the DC analysis above.
        ],
    )


def dummy_sim_tests(
    simulator: SupportedSimulators, skip: Optional[List[AnalysisType]] = None
) -> sd.SimResult:
    """Run the `dummy_sim` input for simulator `simulator`.
    Returns the `sim_data.SimResult` object for any further inspection."""
    if skip is None:
        skip = []

    # Create the "dummy" `SimInput`
    inp = dummy_sim(skip=skip)

    # Run it, requesting in-memory `SimData` results.
    sd_results = sim(
        inp=inp,
        opts=SimOptions(
            simulator=simulator,
            fmt=ResultFormat.SIM_DATA,
            rundir="./scratch",
        ),
    )

    # Check a handful of things about what comes back
    assert isinstance(sd_results, sd.SimResult)
    if AnalysisType.OP not in skip:
        assert isinstance(sd_results[AnalysisType.OP], sd.OpResult)
    if AnalysisType.DC not in skip:
        assert isinstance(sd_results[AnalysisType.DC], sd.DcResult)
    if AnalysisType.TRAN not in skip:
        assert isinstance(sd_results[AnalysisType.TRAN], sd.TranResult)
    if AnalysisType.AC not in skip:
        assert isinstance(sd_results[AnalysisType.AC], sd.AcResult)
    converted_results = sd_results.to_proto()
    assert isinstance(converted_results, vsp.SimResult)

    # Run it again, but requesting VLSIR ProtoBuf output.
    proto_results = sim(
        inp=inp,
        opts=SimOptions(
            simulator=simulator,
            fmt=ResultFormat.VLSIR_PROTO,
        ),
    )
    assert isinstance(proto_results, vsp.SimResult)

    # And return the `sim_data` version for any further inspection.
    return sd_results


def test_spectre_transient_parameters_are_inline():
    """Put analysis parameters inline without reinterpreting literal controls."""

    from vlsirtools.netlist import SpectreNetlister

    inp = dummy_sim(skip=[AnalysisType.OP, AnalysisType.DC, AnalysisType.AC])
    inp.an[0].tran.ctrls.extend(
        (
            vsp.Control(
                param=Param(
                    name="strobeperiod",
                    value=ParamValue(literal="1p"),
                )
            ),
            vsp.Control(literal="save foo"),
        )
    )
    output = StringIO()
    SpectreNetlister(dest=output).write_sim_input(inp)

    tran_line = next(line for line in output.getvalue().splitlines() if line.startswith("tr1 tran"))
    assert tran_line.endswith("strobeperiod=1p")
    assert "\nsave foo\n" in output.getvalue()


def test_spectre_nutbin_parses_inline_first_variable(tmp_path):
    """Parse the variable-header layout emitted by Spectre 24.1."""

    from vlsirtools.spice.spectre import NumType, parse_nutbin

    values = np.asarray(
        [
            0.0,
            0.0,
            0.0,
            1.0e-9,
            -1.0e-3,
            1.0,
        ],
        dtype=">f8",
    )
    raw_path = tmp_path / "spectre.raw"
    raw_path.write_bytes(
        b"Title: \n"
        b"Date: 12:44:32 PM, Thur Aug 20, 2026\n"
        b"Plotname: Transient Analysis `tran': time = (0 s -> 1 ns)\n"
        b"Flags: real\n"
        b"No. Variables:        3\n"
        b"No. Points:       2\n"
        b"Variables:\t0\ttime\ts\n"
        b"\t\t1\txtop.v:p\tA plot=0 grid=0\n"
        b"\t\t2\txtop.x\tV plot=0 grid=0\n"
        b"Binary:\n" + values.tobytes()
    )

    with raw_path.open("rb") as raw_file:
        result = parse_nutbin(raw_file)["tran"]

    assert result.numtype == NumType.REAL
    assert tuple(result.data) == ("time", "xtop.v:p", "xtop.x")
    np.testing.assert_array_equal(result.data["time"], [0.0, 1.0e-9])
    np.testing.assert_array_equal(result.data["xtop.v:p"], [0.0, -1.0e-3])
    np.testing.assert_array_equal(result.data["xtop.x"], [0.0, 1.0])


def test_spectre_nutbin_rejects_truncated_payload(tmp_path):
    """Reject a payload ending after complete records despite a larger header count."""

    from vlsirtools.spice.spectre import parse_nutbin

    # Two complete records, while the header declares three.
    values = np.asarray(
        [
            0.0,
            0.0,
            1.0e-9,
            1.0,
        ],
        dtype=">f8",
    )
    raw_path = tmp_path / "truncated.raw"
    raw_path.write_bytes(
        b"Plotname: Transient Analysis `tran': time = (0 s -> 2 ns)\n"
        b"Flags: real\n"
        b"No. Variables: 2\n"
        b"No. Points: 3\n"
        b"Variables:\t0\ttime\ts\n"
        b"\t\t1\tout\tV\n"
        b"Binary:\n" + values.tobytes()
    )

    with raw_path.open("rb") as raw_file, pytest.raises(ValueError, match=r"declares 3 points.*only 4 values"):
        parse_nutbin(raw_file)


def test_sim_without_analyses_does_not_invoke_backend():
    """Return an empty result when there are no declared analyses."""

    from vlsirtools.spice.base import Sim as BaseSim

    class BackendWhichMustNotRun(BaseSim):
        @classmethod
        def enum(cls):
            return SupportedSimulators.NGSPICE

        def run(self):
            raise AssertionError("empty simulations must not invoke a backend")

    inp = vsp.SimInput()
    sim_data_result = BackendWhichMustNotRun.sim(inp, SimOptions(fmt=ResultFormat.SIM_DATA))
    proto_result = BackendWhichMustNotRun.sim(inp, SimOptions(fmt=ResultFormat.VLSIR_PROTO))

    assert isinstance(sim_data_result, sd.SimResult)
    assert sim_data_result.an == []
    assert isinstance(proto_result, vsp.SimResult)
    assert proto_result.an == []


def test_simulator_arguments_are_per_invocation(monkeypatch):
    """Pass command arguments through one simulation without global mutation."""

    from vlsirtools.spice.ngspice import NGSpiceSim
    from vlsirtools.spice.spectre import SpectreSim
    from vlsirtools.spice.xyce import XyceSim

    options = SimOptions(
        simulator=SupportedSimulators.SPECTRE,
        simulator_args=("+preset=mx", "+mt=4", "+log", "spectre.log"),
    )
    spectre = SpectreSim(inp=vsp.SimInput(), opts=options)
    commands = []
    monkeypatch.setattr(spectre, "run_subprocess", commands.append)

    spectre.run_spectre_process()

    assert commands == [
        [
            "spectre",
            "+preset=mx",
            "+mt=4",
            "+log",
            "spectre.log",
            "-E",
            "-format",
            "nutbin",
            "netlist.scs",
        ]
    ]

    commands.clear()
    ngspice = NGSpiceSim(inp=vsp.SimInput(), opts=options)
    monkeypatch.setattr(ngspice, "run_subprocess", commands.append)
    ngspice.run_sim_process()
    assert commands == [
        [
            "ngspice",
            "--no-spiceinit",
            "+preset=mx",
            "+mt=4",
            "+log",
            "spectre.log",
            "-b",
            "netlist.sp",
            "-r",
            "netlist.raw",
        ]
    ]

    commands.clear()
    xyce = XyceSim(inp=vsp.SimInput(), opts=options)
    monkeypatch.setattr(xyce, "run_subprocess", lambda cmd: commands.append(cmd))
    xyce.run_xyce_process("tran")
    assert commands == [["Xyce", "+preset=mx", "+mt=4", "+log", "spectre.log", "tran.sp"]]

    with pytest.raises(TypeError, match="sequence"):
        SimOptions(simulator_args="+mt=4")


def test_spectre_can_run_without_parsing_results(monkeypatch):
    """Leave Spectre's native raw format intact when results are not requested."""

    from vlsirtools.spice.spectre import SpectreSim

    options = SimOptions(simulator=SupportedSimulators.SPECTRE, fmt=ResultFormat.NONE)
    spectre = SpectreSim(inp=vsp.SimInput(), opts=options)
    commands = []
    monkeypatch.setattr(spectre, "write_netlist", lambda: None)
    monkeypatch.setattr(spectre, "run_subprocess", commands.append)
    monkeypatch.setattr(spectre, "parse_results", lambda: pytest.fail("results must not be parsed"))

    result = spectre.run()

    assert isinstance(result, sd.SimResult)
    assert result.an == []
    assert commands == [["spectre", "-E", "netlist.scs"]]

    with pytest.raises(ValueError, match="only supported by the Spectre backend"):
        SimOptions(simulator=SupportedSimulators.NGSPICE, fmt=ResultFormat.NONE)


@pytest.mark.spectre
def test_spectre1():
    """Test an empty-input call to the `vlsir.spice.Sim` interface to `spectre`."""
    dummy_sim_tests(SupportedSimulators.SPECTRE)


@pytest.mark.xyce
def test_xyce1():
    """Test an empty-input call to the `vlsir.spice.Sim` interface to `xyce`."""
    dummy_sim_tests(SupportedSimulators.XYCE)


@pytest.mark.ngspice
def test_ngspice1():
    """Test an empty-input call to the `vlsir.spice.Sim` interface to `ngspice`."""

    res = dummy_sim_tests(
        SupportedSimulators.NGSPICE,
        skip=[
            AnalysisType.DC
        ],  ## Skip the DC sweep; ngspice doesn't support this kinda sweep
    )

    # Check that the AC frequency vector came back as expected
    # A follow-on from https://github.com/Vlsir/Vlsir/issues/66
    freqs = np.array(
        [
            1.00000000e03,
            1.25892541e03,
            1.58489319e03,
            1.99526231e03,
            2.51188643e03,
            3.16227766e03,
            3.98107171e03,
            5.01187234e03,
            6.30957344e03,
            7.94328235e03,
            1.00000000e04,
            1.25892541e04,
            1.58489319e04,
            1.99526231e04,
            2.51188643e04,
            3.16227766e04,
            3.98107171e04,
            5.01187234e04,
            6.30957344e04,
            7.94328235e04,
            1.00000000e05,
            1.25892541e05,
            1.58489319e05,
            1.99526231e05,
            2.51188643e05,
            3.16227766e05,
            3.98107171e05,
            5.01187234e05,
            6.30957344e05,
            7.94328235e05,
            1.00000000e06,
        ]
    )
    ac_result = res[AnalysisType.AC]
    assert np.allclose(ac_result.freq, freqs)


@pytest.mark.xyce
def test_xyce_import():
    # Just test importing from the `vlsirtools.xyce` path
    # FIXME: probably deprecate this
    from vlsirtools.xyce import sim


@pytest.mark.spectre
def test_spectre_import():
    # Just test importing from the `vlsirtools.spectre` path
    # FIXME: probably deprecate this
    from vlsirtools.spectre import sim


@pytest.mark.ngspice
def test_noise1():
    """Test the Noise analysis"""

    # A very complicated testbench: a voltage source and resistor in parallel.
    pkg = vckt.Package(
        domain="vlsirtools.tests.test_noise1",
        modules=[
            vckt.Module(
                name="noisetb",
                ports=[
                    vckt.Port(direction="NONE", signal="VSS"),
                ],
                signals=[
                    vckt.Signal(name="VSS", width=1),
                    vckt.Signal(name="the_other_node", width=1),
                ],
                instances=[
                    Instance(
                        name="r1",
                        module=Reference(
                            external=QualifiedName(
                                domain="vlsir.primitives", name="resistor"
                            )
                        ),
                        connections=[
                            Connection(
                                portname="p",
                                target=ConnectionTarget(sig="the_other_node"),
                            ),
                            Connection(
                                portname="n", target=ConnectionTarget(sig="VSS")
                            ),
                        ],
                        parameters=[
                            Param(name="r", value=ParamValue(double_value=1e3))
                        ],
                    ),
                    Instance(
                        name="v1",
                        module=Reference(
                            external=QualifiedName(
                                domain="vlsir.primitives", name="vdc"
                            )
                        ),
                        connections=[
                            Connection(
                                portname="p",
                                target=ConnectionTarget(sig="the_other_node"),
                            ),
                            Connection(
                                portname="n", target=ConnectionTarget(sig="VSS")
                            ),
                        ],
                        parameters=[
                            Param(name="dc", value=ParamValue(double_value=0)),
                            Param(name="ac", value=ParamValue(double_value=0)),
                        ],
                    ),
                ],
            )
        ],
    )
    sim_input = vsp.SimInput(
        pkg=pkg,
        top="noisetb",
        an=[
            vsp.Analysis(
                noise=vsp.NoiseInput(
                    analysis_name="noise1",
                    output_p="the_other_node",
                    output_n="VSS",
                    input_source="v1",
                    fstart=1,
                    fstop=1e10,
                    npts=10,
                    ctrls=[],
                )
            )
        ],
    )
    vlsirtools.spice.sim(
        sim_input,
        opts=SimOptions(
            simulator=SupportedSimulators.NGSPICE, fmt=ResultFormat.SIM_DATA
        ),
    )


def test_ngspice_forward_flags():
    """# Regression test for #110
    Handle extra flags like 'forward' in nutbin output."""
    from vlsirtools.spice.ngspice import parse_nutbin_analysis, NumType

    # A mock "nutbin" file content with the problematic "Flags: real forward"
    dummy_raw_content = (
        b"Title: demo\n"
        b"Date: Wed Jan  7 2026\n"
        b"Plotname: DC transfer characteristic\n"
        b"Flags: real forward\n"
        b"No. Variables: 1\n"
        b"No. Points: 1\n"
        b"Variables:\n"
        b"\t0\tv(1)\tvoltage\n"
        b"Binary:\n" + np.array([1.0], dtype="<f8").tobytes()
    )

    fname = "test_ngspice_forward_flags.raw"
    with open(fname, "wb") as tmp:
        tmp.write(dummy_raw_content)

    try:
        with open(fname, "rb") as f:
            # Skip title and date like the main parser
            f.readline()
            f.readline()
            plotname = f.readline().decode("ascii")

            res = parse_nutbin_analysis(f, plotname)
            assert res.numtype == NumType.REAL
            assert res.analysis_name == "Plotname: DC transfer characteristic\n"
            assert "v(1)" in res.data
            assert res.data["v(1)"][0] == 1.0
    finally:
        if os.path.exists(fname):
            os.remove(fname)


@pytest.mark.ngspice
def test_integration_ngspice_nutbin():
    """
    Integration test: run ngspice via vlsirtools, generate nutbin output, and parse it.
    """
    import tempfile
    from vlsirtools.spice.ngspice import parse_nutbin

    # Minimal netlist: voltage source and resistor
    netlist = """
V1 1 0 1
R1 1 0 1k
.control
op
set filetype=binary
write ./test.raw
.endc
.end
"""
    with tempfile.TemporaryDirectory() as tmpdir:
        netlist_path = os.path.join(tmpdir, "test.sp")
        raw_path = os.path.join(tmpdir, "test.raw")
        with open(netlist_path, "w") as f:
            f.write(netlist)
        # Run ngspice in the temp dir so output is placed there
        os.system(f"cd {tmpdir} && ngspice -b test.sp")
        # Parse nutbin
        assert os.path.exists(raw_path), f"ngspice did not produce {raw_path}"
        with open(raw_path, "rb") as f:
            results = parse_nutbin(f)
        # Check results
        assert results, "No results parsed from nutbin file"
        for analysis, res in results.items():
            assert hasattr(res, "data"), "Parsed result missing data field"
            assert any(
                len(arr) > 0 for arr in res.data.values()
            ), "No data in parsed result"
