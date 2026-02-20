"""
# GDSFactory YAML Netlister

Exports `vlsir.circuit.Package` hierarchy to a GDSFactory-style recursive
YAML netlist: `{module_name: {instances, placements, ports, nets}}`.
"""

from __future__ import annotations

import json
from typing import Any

import vlsir
import vlsir.circuit_pb2 as vckt

from .base import (
    ModuleLike,
    Netlister,
    ResolvedModule,
    ResolvedRef,
    SpiceBuiltin,
    SpiceModelRef,
)


class GdsfactoryYamlNetlister(Netlister):
    """Netlist backend emitting a GDSFactory-style recursive YAML netlist."""

    @property
    def enum(self):
        """Get our entry in the `NetlistFormat` enumeration."""
        from . import NetlistFormat

        return NetlistFormat.GDSFACTORY_YAML

    def write_package(self, pkg: vckt.Package) -> None:
        """Write a recursive GDSFactory-style netlist for `pkg`."""

        # First collect externally-defined modules for port-order resolution.
        for emod in pkg.ext_modules:
            self.get_external_module(emod)

        # Cache local modules for local-reference resolution.
        for mod in pkg.modules:
            self.pmodules[mod.name] = mod

        recnet: dict[str, dict[str, Any]] = {}
        for mod in pkg.modules:
            mod_name = self.get_module_name(mod)
            if mod_name in recnet:
                raise RuntimeError(f"Module {mod_name} doubly defined")
            recnet[mod_name] = self.write_module_definition(mod)

        dumped = self.dump_yaml(recnet)
        self.write(dumped)
        if not dumped.endswith("\n"):
            self.write("\n")
        self.flush()

    def write_module_definition(self, module: vckt.Module) -> dict[str, Any]:
        """Create one module-entry in recursive GDSFactory netlist format."""

        self.collect_signals_by_name(module)

        endpoints_by_net: dict[str, list[str]] = {}
        instances: dict[str, dict[str, Any]] = {}

        def add_endpoint(net_name: str, endpoint: str) -> None:
            bucket = endpoints_by_net.setdefault(net_name, [])
            if endpoint not in bucket:
                bucket.append(endpoint)

        # Top-level ports are endpoints without commas.
        for pport in module.ports:
            for net_name in self.expand_signal_bits(pport.signal):
                add_endpoint(net_name, net_name)

        for pinst in module.instances:
            resolved = self.resolve_reference(pinst.module)
            child_module, component_name = self.module_and_name(resolved)

            inst_data: dict[str, Any] = {"component": component_name}
            resolved_params = self.get_instance_params(pinst, child_module)
            if resolved_params:
                inst_data["settings"] = dict(resolved_params.inner)
            instances[pinst.name] = inst_data

            port_widths = self.port_widths(child_module)
            for conn in pinst.connections:
                pwidth = port_widths.get(conn.portname)
                if pwidth is None:
                    raise RuntimeError(
                        f"Unknown port {conn.portname} on instance {pinst.name}"
                    )

                target_bits = self.expand_target_bits(conn.target)
                if len(target_bits) != pwidth:
                    raise RuntimeError(
                        f"Width mismatch on {pinst.name}.{conn.portname}: "
                        f"port width {pwidth}, connection width {len(target_bits)}"
                    )

                for bit_index, net_name in enumerate(target_bits):
                    endpoint = f"{pinst.name},{conn.portname}"
                    if pwidth > 1:
                        endpoint = f"{endpoint}[{bit_index}]"
                    add_endpoint(net_name, endpoint)

        ports: dict[str, str] = {}
        edges: set[tuple[str, str]] = set()
        for endpoints in endpoints_by_net.values():
            uniq_endpoints = sorted(set(endpoints))
            top_level = [ep for ep in uniq_endpoints if "," not in ep]
            instance_ports = [ep for ep in uniq_endpoints if "," in ep]

            for top_port in top_level:
                if instance_ports and top_port not in ports:
                    ports[top_port] = instance_ports[0]

            if len(instance_ports) >= 2:
                anchor = instance_ports[0]
                for other in instance_ports[1:]:
                    p1, p2 = sorted((anchor, other))
                    edges.add((p1, p2))

        nets = [{"p1": p1, "p2": p2} for p1, p2 in sorted(edges)]
        placements = {iname: {} for iname in instances}

        return {
            "instances": instances,
            "placements": placements,
            "ports": ports,
            "nets": nets,
        }

    def module_and_name(self, resolved: ResolvedRef) -> tuple[ModuleLike, str]:
        """Resolve a reference into (module-like, printable component name)."""
        if isinstance(resolved, ResolvedModule):
            return resolved.module, resolved.module_name
        if isinstance(resolved, SpiceModelRef):
            return resolved.module, resolved.model_name
        if isinstance(resolved, SpiceBuiltin):
            return resolved.module, resolved.module.name.name
        raise RuntimeError(f"Unsupported module reference {resolved}")

    def port_widths(self, module: ModuleLike) -> dict[str, int]:
        """Map each port-name to its bit-width."""
        signal_widths = {sig.name: sig.width for sig in module.signals}
        return {
            port.signal: signal_widths.get(port.signal, 1)
            for port in module.ports
        }

    def expand_target_bits(self, ptarget: vckt.ConnectionTarget) -> list[str]:
        """Expand a connection target into one bit-level net-name per bit."""
        stype = ptarget.WhichOneof("stype")
        if stype == "sig":
            return self.expand_signal_bits(ptarget.sig)
        if stype == "slice":
            step = -1 if ptarget.slice.top >= ptarget.slice.bot else 1
            stop = ptarget.slice.bot + step
            return [
                f"{ptarget.slice.signal}[{idx}]"
                for idx in range(ptarget.slice.top, stop, step)
            ]
        if stype == "concat":
            out: list[str] = []
            for part in ptarget.concat.parts:
                out.extend(self.expand_target_bits(part))
            return out
        raise RuntimeError(f"Invalid connection target type {stype}")

    def expand_signal_bits(self, signal_name: str) -> list[str]:
        """Expand signal `signal_name` into one net-name per bit."""
        signal = self.get_signal(signal_name)
        if signal.width <= 1:
            return [signal.name]
        return [f"{signal.name}[{idx}]" for idx in range(signal.width)]

    @classmethod
    def format_prefix(cls, pre: vlsir.SIPrefix) -> str:
        """Format a `SIPrefix` to a string."""
        # YAML has no built-in SI syntax; emit decimal exponents.
        map = {
            vlsir.SIPrefix.YOCTO: "e-24",
            vlsir.SIPrefix.ZEPTO: "e-21",
            vlsir.SIPrefix.ATTO: "e-18",
            vlsir.SIPrefix.FEMTO: "e-15",
            vlsir.SIPrefix.PICO: "e-12",
            vlsir.SIPrefix.NANO: "e-9",
            vlsir.SIPrefix.MICRO: "e-6",
            vlsir.SIPrefix.MILLI: "e-3",
            vlsir.SIPrefix.CENTI: "e-2",
            vlsir.SIPrefix.DECI: "e-1",
            vlsir.SIPrefix.UNIT: "",
            vlsir.SIPrefix.DECA: "e1",
            vlsir.SIPrefix.HECTO: "e2",
            vlsir.SIPrefix.KILO: "e3",
            vlsir.SIPrefix.MEGA: "e6",
            vlsir.SIPrefix.GIGA: "e9",
            vlsir.SIPrefix.TERA: "e12",
            vlsir.SIPrefix.PETA: "e15",
            vlsir.SIPrefix.EXA: "e17",
            vlsir.SIPrefix.ZETTA: "e18",
            vlsir.SIPrefix.YOTTA: "e19",
        }
        if pre not in map:
            raise ValueError(f"Invalid or Unsupported SIPrefix {pre}")

        return map[pre]

    @staticmethod
    def dump_yaml(data: dict[str, Any]) -> str:
        """Serialize `data` to YAML (or JSON fallback if PyYAML is unavailable)."""
        try:
            import yaml
        except ImportError:
            return json.dumps(data, indent=2)
        return yaml.safe_dump(data, sort_keys=False)
