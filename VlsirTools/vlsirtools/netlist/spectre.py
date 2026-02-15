"""
# Spectre-Format Netlister
"""

# Std-Lib Imports
from typing import Union

# Local Imports
import vlsir
import vlsir.circuit_pb2 as vckt
import vlsir.spice_pb2 as vsp

# Import the base-class
from .spectre_spice_shared import SpectreSpiceShared
from .base import ResolvedModule, ResolvedParams, SpiceModelRef, SpiceBuiltin, SpiceType


def map_primitive(rmodule: SpiceBuiltin, paramvals: ResolvedParams) -> str:
    """Map a primitive into Spectre's supported names and parameters.
    Returns the "apparent module name" for instances of the primitive.
    Argument `paramvals` is often modified along the way.

    Note spectre syntax is such that the "apparent module name" can be either of:
    (a) a fixed name per spice-prefix, for basic types (r, c, l, etc), or
    (b) the *model* name for model-based devices, (mos, bjt, diode, etc)"""

    # For voltage sources, add spectre's "type" parameter, and potentially rename several parameters
    if rmodule.spice_type == SpiceType.VSOURCE:
        vname: str = rmodule.module.name.name
        vtypes = dict(
            vdc="dc",
            vpulse="pulse",
            vpwl="pwl",
            vsin="sine",
        )
        if vname not in vtypes:
            msg = f"Invalid or unsupported voltage-source type {vname}"
            raise ValueError(msg)
        paramvals.set("type", vtypes[vname])

        if vname == "vpulse":
            # For pulse sources, need to rename most parameters
            paramvals.rename("v1", "val0")
            paramvals.rename("v2", "val1")
            paramvals.rename("td", "delay")
            paramvals.rename("tr", "rise")
            paramvals.rename("tf", "fall")
            paramvals.rename("tpw", "width")
            paramvals.rename("tper", "period")
        elif vname == "vpwl":
            # Spectre expects `wave=[t0 v0 t1 v1 ...]`.
            wave = paramvals.pop("wave").strip()
            if not wave.startswith("["):
                wave = f"[{wave}]"
            paramvals.set("wave", wave)
        elif vname == "vdc":
            if "ac" in paramvals:
                paramvals.rename("ac", "mag")

    # Mapping from spice-prefix to spectre-name for fixed-name types
    basics = {
        SpiceType.RESISTOR: "resistor",
        SpiceType.CAPACITOR: "capacitor",
        SpiceType.INDUCTOR: "inductor",
        SpiceType.VSOURCE: "vsource",
        SpiceType.ISOURCE: "isource",
        SpiceType.VCVS: "vcvs",
        SpiceType.VCCS: "vccs",
        SpiceType.CCCS: "cccs",
        SpiceType.CCVS: "ccvs",
    }
    if rmodule.spice_type in basics:
        return basics[rmodule.spice_type]

    # Other primitive-types get an "apparent module name" equal to their *model* name.
    model_based = {
        SpiceType.MOS,
        SpiceType.BIPOLAR,
        SpiceType.DIODE,
        SpiceType.TLINE,
    }
    if rmodule.spice_type in model_based:  # Invalid
        msg = f"Internal Error: Invalid model-based primitive {rmodule} should have been resolved to a model"
        raise RuntimeError(msg)

    # Otherwise, unclear what this is or how we got here.
    raise RuntimeError(f"Unsupported or Invalid Primitive {rmodule}")


class SpectreNetlister(SpectreSpiceShared):
    """Spectre-Format Netlister"""

    @property
    def enum(self):
        """Get our entry in the `NetlistFormat` enumeration"""
        from . import NetlistFormat

        return NetlistFormat.SPECTRE

    def write_module_definition(self, module: vckt.Module) -> None:
        """Create a Spectre-format definition for proto-Module `module`"""

        # Create the module name
        module_name = self.get_module_name(module)
        # Check for double-definition
        if module_name in self.module_names:
            raise RuntimeError(f"Module {module_name} doubly defined")
        # Add to our visited lists
        self.module_names.add(module_name)
        self.pmodules[module.name] = module

        # Collect and index vckt.Signals in this Module by name.
        self.collect_signals_by_name(module)

        # Create the sub-circuit definition header
        self.writeln(f"subckt {module_name} ")
        self.indent += 1

        if module.ports:  # Create its ports
            self.writeln(
                "+ "
                + " ".join([self.format_port_decl(pport) for pport in module.ports])
                + " "
            )
        else:
            self.writeln("+ // No ports ")

        # Create its parameters, if defined
        if module.parameters:
            formatted = [self.format_param_decl(pparam) for pparam in module.parameters]
            formatted = " ".join(formatted)
            self.writeln("parameters " + formatted + " ")
        else:
            self.writeln("+ // No parameters ")

        self.writeln("")  # End "header" facets
        # Note nothing need be done for internal signals;
        # spice and spectre create these "out of thin air"

        # Create its instances
        for pinst in module.instances:
            self.write_instance(pinst)
        self.writeln("")

        # Write any netlist-literal content
        self.write_literals(module.literals)

        # Close up the sub-circuit
        self.indent -= 1
        self.writeln("ends \n")

    def write_instance(self, pinst: vckt.Instance) -> None:
        """Create and return a netlist-string for Instance `pinst`"""

        # Initial resolution phase.
        # Start by getting the Module or ExternalModule definition.
        ref = self.resolve_reference(pinst.module)
        module = ref.module

        # Resolve its parameter values, including applying module-level defaults
        resolved_instance_parameters = self.get_instance_params(pinst, module)

        # And sort out its "apparent module name" for netlisting
        if isinstance(ref, ResolvedModule):
            module_name = ref.module_name
        elif isinstance(ref, SpiceModelRef):
            # Spectre format makes "model_name" and "module_name" look identical
            module_name = ref.model_name
        elif isinstance(ref, SpiceBuiltin):
            # Map the primitive into Spectre's built-in element names and parameters
            module_name = map_primitive(ref, resolved_instance_parameters)
        else:
            raise RuntimeError(f"Invalid reference {ref}")

        # Check if compact mode is enabled
        compact = self.opts.compact if self.opts else False

        if compact:
            # Build single-line instance
            parts = [pinst.name]

            # Add ports
            if module.ports:
                port_order = [pport.signal for pport in module.ports]
                connection_targets = {c.portname: c.target for c in pinst.connections}
                pconns = []
                for pname in port_order:
                    pconn = connection_targets.get(pname)
                    if pconn is None:
                        raise RuntimeError(f"Unconnected Port {pname} on {pinst.name}")
                    pconns.append(self.format_connection_target(pconn))
                parts.append("( " + " ".join(pconns) + " )")

            # Add module name
            parts.append(module_name)

            # Add parameters
            if resolved_instance_parameters:
                params = " ".join(
                    f"{k}={v}" for k, v in resolved_instance_parameters.items()
                )
                parts.append(params)

            self.writeln(" ".join(parts))
        else:
            # Original verbose format
            # Create the instance name
            # Note spectre format does not include any `SpiceType` based prefixing
            self.writeln(pinst.name + "")

            if module.ports:
                self.writeln("+ // Ports: ")
                # Get `module`'s port-order
                port_order = [pport.signal for pport in module.ports]
                # And write the Instance ports, in that order
                pconns = []
                connection_targets = {
                    connection.portname: connection.target
                    for connection in pinst.connections
                }
                for pname in port_order:
                    pconn = connection_targets.get(pname, None)
                    if pconn is None:
                        raise RuntimeError(f"Unconnected Port {pname} on {pinst.name}")
                    pconns.append(pconn)
                self.writeln(
                    "+ ( "
                    + " ".join(
                        [self.format_connection_target(pconn) for pconn in pconns]
                    )
                    + " )"
                )
            else:
                self.writeln("+ // No ports ")

            # Write the module-name
            self.writeln("+  " + module_name + " ")

            # Write the instance parameters
            self.write_instance_params(resolved_instance_parameters)

            # And add a post-instance blank line
            self.writeln("")

    def write_instance_params(self, pvals: ResolvedParams) -> None:
        """Write Instance parameters `pvals`"""
        if not pvals:
            return self.writeln("+ // No parameters ")
        # Write its parameter-values
        formatted = " ".join([f"{pname}={pval}" for pname, pval in pvals.items()])
        self.writeln("+  " + formatted + " ")

    def format_concat(self, pconc: vckt.Concat) -> str:
        """Format the Concatenation of several other Connections"""
        out = ""
        for part in pconc.parts:
            out += self.format_connection_target(part) + " "
        return out

    def format_port_decl(self, pport: vckt.Port) -> str:
        """Get a netlist `Port` definition"""
        # In Spectre, as well as most spice, this syntax is the same as referring to the Port.
        return self.format_port_ref(pport)

    def format_port_ref(self, pport: vckt.Port) -> str:
        """Get a netlist `Port` reference"""
        return self.format_signal_ref(self.get_signal(pport.signal))

    @classmethod
    def format_signal_ref(cls, psig: vckt.Signal) -> str:
        """Get a netlist definition for Signal `psig`"""
        if psig.width < 1:
            raise RuntimeError
        if psig.width == 1:  # width==1, i.e. a scalar signal
            return psig.name
        # Vector/ multi "bit" Signal. Creates several spice signals.
        return " ".join(
            [f"{psig.name}{cls.format_bus_bit(k)}" for k in reversed(range(psig.width))]
        )

    @classmethod
    def format_signal_slice(cls, pslice: vckt.Slice) -> str:
        """Get a netlist definition for Signal-Slice `pslice`"""
        base = pslice.signal
        indices = list(reversed(range(pslice.bot, pslice.top + 1)))
        if not len(indices):
            raise RuntimeError(f"Attempting to netlist empty slice {pslice}")
        return " ".join([f"{base}{cls.format_bus_bit(k)}" for k in indices])

    @classmethod
    def format_bus_bit(cls, index: Union[int, str]) -> str:
        """Format-specific string-representation of a bus bit-index"""
        # Spectre netlisting uses an underscore prefix, e.g. `bus_0`
        return "_" + str(index)

    def write_comment(self, comment: str) -> None:
        """While Spectre *can* do a bunch of other comment-styles,
        the canonical one is generally the C-style line comment beginning with `//`."""
        self.writeln(f"// {comment}")

    @classmethod
    def format_sim_dut(cls, module_name: str) -> str:
        """# Format the top-level DUT instance for module name `module_name`."""
        return f"xtop 0 {module_name} // Top-Level DUT \n"

    def write_sim_header(self, inp: vsp.SimInput) -> None:
        """# Write header commentary for a `SimInput`"""

        # Note on this line here:
        self.writeln("simulator lang=spectre \n")
        # There's no need to write `simulator lang=spectre` in a file who's extension is `.scs`.
        # But I guess this doesn't hurt? And continues to work if the file gets copied to another extension?

        # Write the base-class version, it includes some nice commentary
        super().write_sim_header(inp)

        # FIXME: do we *really* want this global-zero?
        self.writeln("global 0")
        self.writeln("")

    def write_include(self, inc: vsp.Include) -> None:
        """# Write an `Include` statement"""
        self.writeln(f'include "{inc.path}"')

    def write_lib_include(self, lib: vsp.LibInclude) -> None:
        """# Write a `LibInclude` statement"""
        txt = f'include "{lib.path}" section={lib.section}'
        self.writeln(txt)

    def write_save(self, save: vsp.Save) -> None:
        """# Write a `Save` statement"""
        # Spectre uses `save` statement for specifying signals to save
        # Mode: NONE means save nothing, ALL means save everything
        if save.mode == vsp.Save.SaveMode.ALL:
            self.writeln("save *")
        elif save.mode == vsp.Save.SaveMode.NONE:
            pass  # Don't save anything
        elif save.signal:
            # Save specific signal
            signals = " ".join(s.strip() for s in save.signal.split(",") if s.strip())
            self.writeln(f"save {signals}")

    def write_meas(self, meas: vsp.Meas) -> None:
        """# Write a SPICE-style `.meas` wrapped in Spectre language switches."""
        self.writeln("simulator lang=spice")
        txt = f".meas {meas.analysis_type} {meas.name} {meas.expr}"
        self.writeln(txt)
        self.writeln("simulator lang=spectre")

    def write_sim_param(self, param: vlsir.Param) -> None:
        """# Write a simulation-level parameter"""
        txt = f"parameters  {param.name}={self.get_param_value(param.value)}"
        self.writeln(txt)

    def write_sim_option(self, opt: vsp.SimOptions) -> None:
        """# Write a simulation option"""
        # Spectre syntax: <name> options <option>=<value>
        self.writeln(f"simOpts options {opt.name}={self.get_param_value(opt.value)}")

    def write_ac(self, an: vsp.AcInput) -> None:
        """# Write an AC analysis."""

        if not an.analysis_name:
            raise RuntimeError(f"Analysis name required for {an}")
        if len(an.ctrls):
            raise NotImplementedError  # FIXME!

        # Unpack the analysis / sweep content
        fstart = an.fstart
        if fstart <= 0:
            raise ValueError(f"Invalid `fstart` {fstart}")
        fstop = an.fstop
        if fstop <= 0:
            raise ValueError(f"Invalid `fstop` {fstop}")
        npts = an.npts
        if npts <= 0:
            raise ValueError(f"Invalid `npts` {npts}")

        # Write the analysis command
        line = f"{an.analysis_name} ac start={fstart} stop={fstop} dec={npts}"
        self.writeln(line)

    def write_dc(self, an: vsp.DcInput) -> None:
        """# Write a DC analysis."""

        if not an.analysis_name:
            raise RuntimeError(f"Analysis name required for {an}")
        if len(an.ctrls):
            raise NotImplementedError  # FIXME!

        # Write the analysis command
        param = an.indep_name
        ## Interpret the sweep
        sweep_type = an.sweep.WhichOneof("tp")
        if sweep_type == "linear":
            sweep = an.sweep.linear
            line = f"{an.analysis_name} dc param={param} start={sweep.start} stop={sweep.stop} step={sweep.step}"
            self.writeln(line)
        elif sweep_type == "points":
            sweep = an.sweep.points
            line = f"{an.analysis_name} dc values=[{sweep.points}]"
            self.writeln(line)
        elif sweep_type == "log":
            raise NotImplementedError
        else:
            raise ValueError("Invalid sweep type")

    def write_op(self, an: vsp.OpInput) -> None:
        """# Write an operating point analysis."""

        if not an.analysis_name:
            raise RuntimeError(f"Analysis name required for {an}")
        if len(an.ctrls):
            raise NotImplementedError  # FIXME!

        self.writeln(f"{an.analysis_name} dc oppoint=rawfile")

    def write_tran(self, an: vsp.TranInput) -> None:
        """# Write a transient analysis."""

        if not an.analysis_name:
            raise RuntimeError(f"Analysis name required for {an}")

        # Build the tran command parts
        parts = [f"{an.analysis_name} tran stop={an.tstop}"]

        # Add tstep if specified
        if an.tstep > 0:
            parts.append(f"step={an.tstep}")

        # Handle initial conditions inline
        if an.ic:
            ic_str = " ".join(f"{sig}={val}" for sig, val in an.ic.items())
            parts.append(f"ic {ic_str}")

        # Add noise option (Spectre-specific)
        if an.noise:
            parts.append("isnoisy=yes")

        self.writeln(" ".join(parts))

        # Write control elements for this analysis
        for ctrl in an.ctrls:
            self.write_control_element(ctrl)

    def write_noise(self, an: vsp.NoiseInput) -> None:
        """# Write a noise analysis."""
        raise NotImplementedError

    def write_monte(self, an: vsp.MonteInput) -> None:
        """# Write a Spectre-native Monte Carlo analysis block.

        Emits:
            <name> montecarlo numruns=<npts> seed=<seed> variations=all {
                <controls>
                <child analyses>
            }
        """
        if not an.analysis_name:
            raise RuntimeError(f"Analysis name required for {an}")

        # Build the montecarlo header
        parts = [f"{an.analysis_name} montecarlo"]
        parts.append(f"numruns={an.npts}")
        if an.seed:
            parts.append(f"seed={an.seed}")
        parts.append("variations=all")
        self.writeln(" ".join(parts) + " {")

        # Indent the block contents
        self.indent += 1

        # Write control elements for the Monte Carlo block
        for ctrl in an.ctrls:
            self.write_control_element(ctrl)

        # Write each child analysis
        for child_an in an.an:
            self.write_analysis(child_an)

        # Close the block
        self.indent -= 1
        self.writeln("}")
        self.writeln("")

    def write_sweep(self, an: vsp.SweepInput) -> None:
        """# Write a Spectre-native sweep analysis block.

        Emits:
            <name> sweep param=<variable> values=[...] or start/stop/step {
                <controls>
                <child analyses>
            }
        """
        if not an.analysis_name:
            raise RuntimeError(f"Analysis name required for {an}")

        # Build the sweep header
        parts = [f"{an.analysis_name} sweep"]
        parts.append(f"param={an.variable}")

        # Interpret the sweep values
        sweep_type = an.sweep.WhichOneof("tp")
        if sweep_type == "linear":
            sweep = an.sweep.linear
            parts.append(f"start={sweep.start}")
            parts.append(f"stop={sweep.stop}")
            parts.append(f"step={sweep.step}")
        elif sweep_type == "log":
            sweep = an.sweep.log
            parts.append(f"start={sweep.start}")
            parts.append(f"stop={sweep.stop}")
            parts.append(f"dec={int(sweep.npts)}")
        elif sweep_type == "points":
            sweep = an.sweep.points
            pts_str = " ".join(str(p) for p in sweep.points)
            parts.append(f"values=[{pts_str}]")
        else:
            raise ValueError("Invalid sweep type")

        self.writeln(" ".join(parts) + " {")

        # Indent the block contents
        self.indent += 1

        # Write control elements for the sweep block
        for ctrl in an.ctrls:
            self.write_control_element(ctrl)

        # Write each child analysis
        for child_an in an.an:
            self.write_analysis(child_an)

        # Close the block
        self.indent -= 1
        self.writeln("}")
        self.writeln("")
