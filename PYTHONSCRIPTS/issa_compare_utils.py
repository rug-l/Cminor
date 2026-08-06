"""
Helpers for ISSA full vs reduced NetCDF verification.
"""

from __future__ import annotations

import csv
import json
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable

import netCDF4 as nc
import numpy as np


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------



def parse_namelist_value(text: str, key: str) -> str | None:
    """Return first assignment value for `key = ...` in namelist text."""
    pattern = re.compile(
        rf"^\s*{re.escape(key)}\s*=\s*([\"'])(.*?)\1",
        re.IGNORECASE | re.MULTILINE,
    )
    m = pattern.search(text)
    if m:
        return m.group(2).strip()
    pattern2 = re.compile(
        rf"^\s*{re.escape(key)}\s*=\s*(\S+)",
        re.IGNORECASE | re.MULTILINE,
    )
    m2 = pattern2.search(text)
    if m2:
        return m2.group(1).strip().rstrip(",")
    return None


def parse_run_file(run_path: Path) -> dict[str, str]:
    text = run_path.read_text(errors="replace")
    keys = (
        "RedCtrlFile",
        "NetCdfFile",
        "InitFile",
        "SysFile",
        "tBegin",
        "tEnd",
        "Temperature0",
        "AtolGas",
        "AtolAqua",
    )
    out = {k: parse_namelist_value(text, k) for k in keys}
    out["_path"] = str(run_path)
    return {k: v for k, v in out.items() if v is not None}


def _comment_group_name(line: str) -> str | None:
    """Extract group label from a `# ...` comment line."""
    stripped = line.strip()
    if not stripped.startswith("#"):
        return None
    body = stripped[1:].strip()
    if not body:
        return None
    if body.startswith("####") and body.rstrip().endswith("####"):
        inner = body.strip("#").strip()
        return inner or None
    return body


def parse_target_species_groups(ctrl_path: Path) -> list[tuple[str, list[str]]]:
    """Return ordered (group_name, species_list) from TARGET_SPC block."""
    lines = ctrl_path.read_text(errors="replace").splitlines()
    groups: list[tuple[str, list[str]]] = []
    in_block = False
    current_group: str | None = None
    current_species: list[str] = []

    def flush() -> None:
        nonlocal current_species
        if current_species:
            name = current_group if current_group else current_species[0]
            groups.append((name, list(current_species)))
            current_species = []

    for raw in lines:
        line = raw.strip()
        upper = line.upper()
        if upper == "TARGET_SPC":
            in_block = True
            continue
        if upper == "END_TARGET_SPC":
            flush()
            break
        if not in_block:
            continue
        if line.startswith("#"):
            flush()
            grp = _comment_group_name(line)
            if grp is not None:
                current_group = grp
            continue
        if not line:
            continue
        current_species.append(line.split()[0])

    if not groups and not in_block:
        raise ValueError(f"No TARGET_SPC block in {ctrl_path}")
    return groups


def parse_target_species(ctrl_path: Path) -> list[str]:
    return [sp for _, species in parse_target_species_groups(ctrl_path) for sp in species]


def parse_species_families(ctrl_path: Path) -> list[list[str]]:
    """Return species families from SPC_FAMILIES block (one list per line)."""
    lines = ctrl_path.read_text(errors="replace").splitlines()
    families: list[list[str]] = []
    in_block = False
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        upper = line.upper()
        if upper == "SPC_FAMILIES":
            in_block = True
            continue
        if upper == "END_SPC_FAMILIES":
            break
        if in_block:
            members = line.split()
            if members:
                families.append(members)
    return families


def check_run_parity(full_run: dict, reduced_run: dict) -> list[str]:
    warnings: list[str] = []
    for key in ("tBegin", "tEnd", "Temperature0", "AtolGas", "AtolAqua", "InitFile"):
        fv, rv = full_run.get(key), reduced_run.get(key)
        if fv is not None and rv is not None and fv != rv:
            warnings.append(f"  {key}: full={fv!r}  reduced={rv!r}")
    return warnings


# ---------------------------------------------------------------------------
# NetCDF species → variable name
# ---------------------------------------------------------------------------

def sanitize_species_name(name: str) -> str:
    """Match NetCDF_Mod check_name_slash + check_name_bracket."""
    out = name.strip()
    while "/" in out:
        out = out.replace("/", "_", 1)
    if out.startswith("["):
        out = "_" + out
    return out


def nc_var_candidates(species: str) -> list[str]:
    name = sanitize_species_name(species)
    cands: list[str] = []
    seen: set[str] = set()

    def add(v: str) -> None:
        if v not in seen:
            seen.add(v)
            cands.append(v)

    add(name)
    if name.startswith("a"):
        base = name[1:]
        add(f"{name}_1_m3")
        add(f"{name}_1")
        add(f"{name}_Sum")
        add(f"a{base}_1_m3")
        add(f"a{base}_1")
        add(f"a{base}_Sum")
    else:
        add(f"a{name}_1_m3")
        add(f"a{name}_1")
        add(f"a{name}_Sum")
        add(f"{name}_1_m3")
        add(f"{name}_1")
        add(f"{name}_Sum")
    return cands


def resolve_nc_variable(species: str, keys_full: set[str], keys_red: set[str]) -> str | None:
    for cand in nc_var_candidates(species):
        if cand in keys_full and cand in keys_red:
            return cand
    return None


def read_nc_var_units(ds: nc.Dataset, var: str) -> str:
    return str(getattr(ds.variables[var], "units", "") or "")


def classify_nc_phase(units: str, nc_var: str) -> str:
    """Return 'gas' or 'aqua' from NetCDF metadata."""
    u = units.strip().lower()
    var = nc_var.lower()
    if u == "mol/l" or var.endswith("_l"):
        return "aqua"
    if u in ("molec/cm3", "mol/mol"):
        return "gas"
    if u == "mol/m3" and (var.endswith("_m3") or var.endswith("_sum")):
        return "aqua"
    if u == "mol/m3":
        return "gas"
    if re.search(r"_\d+_(l|m3)$", var) or var.endswith("_sum"):
        return "aqua"
    return "gas"


def format_unit_for_axis(units: str) -> str:
    labels = {
        "molec/cm3": r"molec cm$^{-3}$",
        "mol/l": r"mol L$^{-1}$",
        "mol/m3": r"mol m$^{-3}$",
        "mol/mol": r"mol mol$^{-1}$",
    }
    return labels.get(units.strip().lower(), units)


# ---------------------------------------------------------------------------
# Time series loading and alignment
# ---------------------------------------------------------------------------

def load_time_seconds(ds: nc.Dataset) -> np.ndarray:
    if "time" not in ds.variables:
        raise ValueError("NetCDF missing 'time' variable")
    return np.array(ds["time"][:], dtype=float)


def load_series(ds: nc.Dataset, var: str, n_t: int) -> np.ndarray:
    if var not in ds.variables:
        raise KeyError(var)
    v = ds.variables[var]
    if len(v.shape) != 1 or v.shape[0] != n_t:
        raise ValueError(f"Variable {var!r} is not a 1D time series of length {n_t}")
    return np.array(v[:], dtype=float)


def align_time_series(
    t_full: np.ndarray,
    y_full: np.ndarray,
    t_red: np.ndarray,
    y_red: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Intersect time ranges; interpolate reduced onto full grid inside overlap."""
    t_lo = max(float(t_full[0]), float(t_red[0]))
    t_hi = min(float(t_full[-1]), float(t_red[-1]))
    if t_hi <= t_lo:
        raise ValueError("No overlapping time range between full and reduced NetCDF")

    mask = (t_full >= t_lo) & (t_full <= t_hi)
    t_out = t_full[mask]
    yf_out = y_full[mask]

    if len(t_red) == len(t_full) and np.allclose(t_red, t_full):
        yr_out = y_red[mask]
    else:
        yr_interp = np.interp(t_out, t_red, y_red)
        yr_out = yr_interp

    return t_out, yf_out, yr_out


# NetCDF output floor (NetCDF_Mod.f90 y_Min)
NCDF_Y_MIN = 1e-40


def parse_fortran_float(text: str) -> float:
    return float(text.strip().lower().replace("d", "e"))


def run_comparison_atol(full_run: dict, reduced_run: dict) -> tuple[float, float, float]:
    """Return (atol_gas, atol_aqua, atol) using max across both runs."""
    atol_gas_vals: list[float] = []
    atol_aqua_vals: list[float] = []
    for run in (full_run, reduced_run):
        if run.get("AtolGas"):
            atol_gas_vals.append(parse_fortran_float(run["AtolGas"]))
        if run.get("AtolAqua"):
            atol_aqua_vals.append(parse_fortran_float(run["AtolAqua"]))
    atol_gas = max(atol_gas_vals) if atol_gas_vals else 0.0
    atol_aqua = max(atol_aqua_vals) if atol_aqua_vals else 0.0
    atol = max(atol_gas, atol_aqua)
    return atol_gas, atol_aqua, atol


def peak_concentration(full: np.ndarray, reduced: np.ndarray) -> float:
    return float(max(np.nanmax(np.abs(full)), np.nanmax(np.abs(reduced))))


def comparison_threshold(atol: float, min_peak: float) -> float:
    return max(atol, min_peak)


def classify_tier(peak: float, threshold: float, y_min: float = NCDF_Y_MIN) -> str:
    if peak <= y_min:
        return "negligible"
    if peak < threshold:
        return "trace"
    return "active"


@dataclass
class RobustAxisLimits:
    """Axis limits from robust percentiles (xarray plotting ``robust=True``)."""

    lo: float
    hi: float
    has_lower_outliers: bool
    has_upper_outliers: bool


def robust_axis_limits(
    values: np.ndarray,
    *,
    log_scale: bool = False,
    low_pct: float = 2.0,
    high_pct: float = 98.0,
    margin_frac: float = 0.05,
) -> RobustAxisLimits | None:
    """
    Compute y-limits from percentiles, ignoring extreme outliers.

    Matches xarray's ``robust=True`` (2nd and 98th percentile by default).
    """
    vals = values[np.isfinite(values)]
    if log_scale:
        vals = vals[vals > 0]
    if vals.size == 0:
        return None

    if log_scale:
        logv = np.log10(vals)
        lo_p, hi_p = np.percentile(logv, [low_pct, high_pct])
        if not np.isfinite(lo_p) or not np.isfinite(hi_p) or lo_p >= hi_p:
            lo_p, hi_p = float(np.min(logv)), float(np.max(logv))
        pad = max(margin_frac * (hi_p - lo_p), 0.05)
        lo, hi = 10.0 ** (lo_p - pad), 10.0 ** (hi_p + pad)
        has_lo = bool(np.any(logv < lo_p))
        has_hi = bool(np.any(logv > hi_p))
    else:
        lo_p, hi_p = np.percentile(vals, [low_pct, high_pct])
        if not np.isfinite(lo_p) or not np.isfinite(hi_p) or lo_p >= hi_p:
            lo_p, hi_p = float(np.min(vals)), float(np.max(vals))
        pad = margin_frac * (hi_p - lo_p) if hi_p > lo_p else max(abs(hi_p), 1.0) * margin_frac
        lo, hi = lo_p - pad, hi_p + pad
        has_lo = bool(np.any(vals < lo_p))
        has_hi = bool(np.any(vals > hi_p))

    return RobustAxisLimits(lo=lo, hi=hi, has_lower_outliers=has_lo, has_upper_outliers=has_hi)


def safe_rel_err(full: np.ndarray, reduced: np.ndarray, eps: float) -> np.ndarray:
    denom = np.maximum(np.abs(full), eps)
    return np.abs(full - reduced) / denom


def devmax2(full: np.ndarray, reduced: np.ndarray, time_s: np.ndarray, transient_s: float = 1800.0) -> float:
    """
    Maximum relative deviation of daily maxima (Mauersberger 2005, Eq. 8).
    Skip first `transient_s` when locating daily maxima.
    """
    if full.size == 0:
        return float("nan")

    denom = np.maximum(np.abs(full), np.abs(reduced))
    denom = np.maximum(denom, 1e-300)
    dev = np.abs(full - reduced) / denom

    t = time_s.astype(float)
    t0 = float(t[0]) + transient_s
    day_s = 86400.0

    t_start = t0
    t_end = float(t[-1])
    if t_end <= t_start:
        return float(np.max(dev))

    day_max_dev: list[float] = []
    while t_start < t_end:
        t_day_end = min(t_start + day_s, t_end + 1.0)
        in_day = (t >= t_start) & (t < t_day_end)
        if np.any(in_day):
            day_max_dev.append(float(np.max(dev[in_day])))
        t_start += day_s

    return float(max(day_max_dev)) if day_max_dev else float(np.max(dev))


@dataclass
class SpeciesMetrics:
    species: str
    nc_var: str
    phase: str
    nc_units: str
    n_points: int
    peak_conc: float
    tier: str
    abs_max: float
    abs_mean: float
    rel_max: float
    rel_mean: float
    rel_p95: float
    devmax2: float
    trace_pass: bool | None = None


@dataclass
class ComparisonResult:
    label: str
    full_nc: str
    reduced_nc: str
    ctrl_file: str
    peak_threshold: float = 0.0
    atol: float = 0.0
    resolved: list[SpeciesMetrics] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def compute_species_metrics(
    species: str,
    nc_var: str,
    phase: str,
    nc_units: str,
    t: np.ndarray,
    full: np.ndarray,
    reduced: np.ndarray,
    eps: float,
    peak_threshold: float,
    atol: float,
) -> SpeciesMetrics:
    peak = peak_concentration(full, reduced)
    tier = classify_tier(peak, peak_threshold)
    abs_diff = np.abs(full - reduced)
    abs_max = float(np.max(abs_diff))
    abs_mean = float(np.mean(abs_diff))

    if tier == "active":
        rel = safe_rel_err(full, reduced, eps)
        rel_max = float(np.max(rel))
        rel_mean = float(np.mean(rel))
        rel_p95 = float(np.percentile(rel, 95))
        dev = devmax2(full, reduced, t)
        trace_pass = None
    else:
        rel_max = rel_mean = rel_p95 = dev = float("nan")
        trace_pass = abs_max <= atol if tier == "trace" else None

    return SpeciesMetrics(
        species=species,
        nc_var=nc_var,
        phase=phase,
        nc_units=nc_units,
        n_points=int(t.size),
        peak_conc=peak,
        tier=tier,
        abs_max=abs_max,
        abs_mean=abs_mean,
        rel_max=rel_max,
        rel_mean=rel_mean,
        rel_p95=rel_p95,
        devmax2=dev,
        trace_pass=trace_pass,
    )


def compare_datasets(
    full_nc: Path,
    reduced_nc: Path,
    species_list: Iterable[str],
    eps: float = 1e-12,
    peak_threshold: float = 1e-7,
    atol: float = 1e-7,
) -> tuple[list[SpeciesMetrics], list[str], np.ndarray]:
    """
    Load both NetCDF files, resolve species, align times, return metrics + skipped + time [h].
    """
    skipped: list[str] = []
    metrics: list[SpeciesMetrics] = []
    time_h: np.ndarray | None = None

    with nc.Dataset(full_nc) as ds_full, nc.Dataset(reduced_nc) as ds_red:
        keys_full = set(ds_full.variables.keys())
        keys_red = set(ds_red.variables.keys())
        t_full = load_time_seconds(ds_full)
        t_red = load_time_seconds(ds_red)
        n_full = t_full.size
        n_red = t_red.size

        for sp in species_list:
            var = resolve_nc_variable(sp, keys_full, keys_red)
            if var is None:
                skipped.append(sp)
                continue
            try:
                y_full = load_series(ds_full, var, n_full)
                y_red = load_series(ds_red, var, n_red)
                t_aligned, yf, yr = align_time_series(t_full, y_full, t_red, y_red)
                if time_h is None:
                    time_h = t_aligned / 3600.0
                nc_units = read_nc_var_units(ds_full, var)
                phase = classify_nc_phase(nc_units, var)
                m = compute_species_metrics(
                    sp, var, phase, nc_units, t_aligned, yf, yr, eps, peak_threshold, atol
                )
                metrics.append(m)
            except (KeyError, ValueError) as exc:
                skipped.append(f"{sp} ({exc})")

    if time_h is None:
        raise ValueError("No resolvable TARGET_SPC species in both NetCDF files")

    return metrics, skipped, time_h


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def write_summary_csv(path: Path, metrics: list[SpeciesMetrics]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(asdict(metrics[0]).keys()) if metrics else []
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for m in metrics:
            w.writerow(asdict(m))


def write_summary_json(path: Path, result: ComparisonResult) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    active = [m for m in result.resolved if m.tier == "active"]
    trace = [m for m in result.resolved if m.tier == "trace"]
    negligible = [m for m in result.resolved if m.tier == "negligible"]
    payload = {
        "label": result.label,
        "full_nc": result.full_nc,
        "reduced_nc": result.reduced_nc,
        "ctrl_file": result.ctrl_file,
        "peak_threshold": result.peak_threshold,
        "atol": result.atol,
        "n_resolved": len(result.resolved),
        "n_active": len(active),
        "n_trace": len(trace),
        "n_negligible": len(negligible),
        "n_skipped": len(result.skipped),
        "skipped": result.skipped,
        "warnings": result.warnings,
        "species": [asdict(m) for m in result.resolved],
    }
    if active:
        worst_dev = max(active, key=lambda m: m.devmax2)
        worst_rel = max(active, key=lambda m: m.rel_max)
        payload["worst_devmax2"] = {"species": worst_dev.species, "value": worst_dev.devmax2}
        payload["worst_rel_max"] = {"species": worst_rel.species, "value": worst_rel.rel_max}
    if trace:
        payload["trace_failures"] = [
            {"species": m.species, "abs_max": m.abs_max}
            for m in trace
            if m.trace_pass is False
        ]
    path.write_text(json.dumps(payload, indent=2))


def print_summary_table(metrics: list[SpeciesMetrics], top_n: int = 15) -> None:
    if not metrics:
        print("[issa-verify] no resolved species")
        return

    active = [m for m in metrics if m.tier == "active"]
    trace = [m for m in metrics if m.tier == "trace"]
    negligible = [m for m in metrics if m.tier == "negligible"]

    if active:
        ranked = sorted(active, key=lambda m: m.devmax2, reverse=True)
        print(f"\nActive species (Mauersberger metrics), n={len(active)}")
        print(f"{'Species':<16} {'NC var':<20} {'devmax2':>10} {'rel_max':>10} {'abs_max':>12}")
        print("-" * 72)
        for m in ranked[:top_n]:
            print(
                f"{m.species:<16} {m.nc_var:<20} {m.devmax2:10.4f} {m.rel_max:10.4f} {m.abs_max:12.4e}"
            )
        if len(ranked) > top_n:
            print(f"  ... and {len(ranked) - top_n} more active species")

    if trace:
        ranked_trace = sorted(trace, key=lambda m: m.abs_max, reverse=True)
        print(f"\nTrace species (abs_max vs atol), n={len(trace)}")
        print(f"{'Species':<16} {'peak_conc':>12} {'abs_max':>12} {'pass':>6}")
        print("-" * 50)
        for m in ranked_trace[:top_n]:
            passed = "yes" if m.trace_pass else "no"
            print(f"{m.species:<16} {m.peak_conc:12.4e} {m.abs_max:12.4e} {passed:>6}")
        if len(ranked_trace) > top_n:
            print(f"  ... and {len(ranked_trace) - top_n} more trace species")

    if negligible:
        print(f"\nNegligible species (peak <= {NCDF_Y_MIN:.0e}), n={len(negligible)}")
        names = ", ".join(m.species for m in negligible[:20])
        print(f"  {names}")
        if len(negligible) > 20:
            print(f"  ... and {len(negligible) - 20} more")
