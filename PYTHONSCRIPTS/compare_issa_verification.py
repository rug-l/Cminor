"""
Compare full vs ISSA-reduced Cminor NetCDF outputs.

Usage:
  python PYTHONSCRIPTS/compare_issa_verification.py \\
    --full-run    RUN/TESTRUN/MCM+CAPRAM/MCM+CAPRAM_issa.run \\
    --reduced-run RUN/TESTRUN/MCM+CAPRAM/MCM+CAPRAM_issa_reduced.run \\
    --label       MCM32_CAPRAM40_issa \\
    --out-dir     RUN/TESTRUN/diagnostics
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

from matplotlib import transforms as mtransforms

from issa_compare_utils import (
    ComparisonResult,
    NCDF_Y_MIN,
    RobustAxisLimits,
    align_time_series,
    check_run_parity,
    classify_nc_phase,
    compare_datasets,
    comparison_threshold,
    format_unit_for_axis,
    load_series,
    load_time_seconds,
    parse_run_file,
    parse_species_families,
    parse_target_species,
    parse_target_species_groups,
    print_summary_table,
    read_nc_var_units,
    resolve_nc_variable,
    robust_axis_limits,
    run_comparison_atol,
    write_summary_csv,
    write_summary_json,
)

# resolved[species] = (full, reduced, nc_var, phase, nc_units)
ResolvedSeries = dict[str, tuple[np.ndarray, np.ndarray, str, str, str]]
# panel = (title, species, full[], red[], phases[], unit_labels[])
PlotPanel = tuple[str, list[str], list[np.ndarray], list[np.ndarray], list[str], list[str]]

MAX_SPECIES_PER_PANEL = 2
MIN_PLOTTABLE_FRAC = 0.05
PLOT_ROBUST_LOW_PCT = 2.0
PLOT_ROBUST_HIGH_PCT = 98.0


def _resolve_path(base: Path, p: str | None) -> Path | None:
    if p is None:
        return None
    path = Path(p)
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def load_species_series(
    full_nc: Path,
    reduced_nc: Path,
    species_list: list[str],
) -> tuple[list[str], list[np.ndarray], list[np.ndarray], np.ndarray, ResolvedSeries]:
    """Return time [h] and resolved series dict with phase/units metadata."""
    time_h: np.ndarray | None = None
    resolved: ResolvedSeries = {}

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
                continue
            y_full = load_series(ds_full, var, n_full)
            y_red = load_series(ds_red, var, n_red)
            t_aligned, yf, yr = align_time_series(t_full, y_full, t_red, y_red)
            if time_h is None:
                time_h = t_aligned / 3600.0
            nc_units = read_nc_var_units(ds_full, var)
            phase = classify_nc_phase(nc_units, var)
            resolved[sp] = (yf, yr, var, phase, nc_units)

    if time_h is None:
        raise ValueError("No plottable species resolved in both NetCDF files")
    return list(resolved.keys()), time_h, resolved


def _plottable_fraction(y: np.ndarray, floor: float = NCDF_Y_MIN) -> float:
    if y.size == 0:
        return 0.0
    return float(np.mean(y > floor * 10))


def _is_plottable(yf: np.ndarray, yr: np.ndarray) -> bool:
    return (
        _plottable_fraction(yf) >= MIN_PLOTTABLE_FRAC
        or _plottable_fraction(yr) >= MIN_PLOTTABLE_FRAC
    )


def _prepare_plot_series(y: np.ndarray, floor: float = NCDF_Y_MIN) -> np.ndarray:
    """Mask NetCDF floor and isolated startup/end spikes (plot-only)."""
    out = y.astype(float).copy()
    out[out <= floor * 10] = np.nan
    if out.size > 3:
        interior = out[1:-1]
        med = np.nanmedian(interior)
        if np.isfinite(med) and med > 0:
            for idx in (0, -1):
                if np.isfinite(out[idx]) and out[idx] > 10 * med:
                    out[idx] = np.nan
    return out


def expand_groups_to_panels(
    groups: list[tuple[str, list[str]]],
    resolved: ResolvedSeries,
    max_per_panel: int = MAX_SPECIES_PER_PANEL,
) -> list[tuple[str, list[str]]]:
    """Split by phase and chunk crowded groups for readable subplots."""
    expanded: list[tuple[str, list[str]]] = []
    for group_name, members in groups:
        gas: list[str] = []
        aqua: list[str] = []
        for sp in members:
            if sp not in resolved:
                continue
            yf, yr, _, phase, _ = resolved[sp]
            if not _is_plottable(yf, yr):
                continue
            (gas if phase == "gas" else aqua).append(sp)

        phase_buckets: list[tuple[str, list[str]]] = []
        if gas:
            phase_buckets.append(("gas", gas))
        if aqua:
            phase_buckets.append(("aqua", aqua))
        if not phase_buckets:
            continue

        for phase_name, sps in phase_buckets:
            n_parts = math.ceil(len(sps) / max_per_panel)
            for chunk_i in range(0, len(sps), max_per_panel):
                chunk = sps[chunk_i : chunk_i + max_per_panel]
                part = chunk_i // max_per_panel + 1
                if len(phase_buckets) > 1 or n_parts > 1:
                    if n_parts > 1:
                        title = f"{group_name} ({phase_name} {part}/{n_parts})"
                    else:
                        title = f"{group_name} ({phase_name})"
                else:
                    title = group_name
                expanded.append((title, chunk))
    return expanded


def balanced_subplot_grid(n: int, ncols_override: int | None = None) -> tuple[int, int]:
    """Return (nrows, ncols) with nrows ≈ ncols and nrows * ncols >= n."""
    if n <= 0:
        return 1, 1
    if ncols_override is not None:
        ncols = max(1, min(ncols_override, n))
        return math.ceil(n / ncols), ncols
    ncols = int(math.ceil(math.sqrt(n)))
    nrows = int(math.ceil(n / ncols))
    return nrows, ncols


def _finite_values(series: list[np.ndarray]) -> np.ndarray:
    chunks = [a[np.isfinite(a)] for a in series if a.size]
    return np.concatenate(chunks) if chunks else np.array([], dtype=float)


def _panel_use_log(series: list[np.ndarray]) -> bool:
    vals = _finite_values(series)
    pos = vals[vals > 0]
    if pos.size < 2:
        return False
    return float(np.max(pos) / np.min(pos)) > 10.0


def _apply_axis_limits(
    ax,
    series: list[np.ndarray],
    log_scale: bool,
    robust: bool = True,
    low_pct: float = PLOT_ROBUST_LOW_PCT,
    high_pct: float = PLOT_ROBUST_HIGH_PCT,
) -> None:
    vals = _finite_values(series)
    if vals.size == 0:
        return

    if robust:
        limits = robust_axis_limits(
            vals,
            log_scale=log_scale,
            low_pct=low_pct,
            high_pct=high_pct,
        )
        if limits is None:
            return
        ax.set_ylim(limits.lo, limits.hi)
        _mark_ylim_outliers(ax, limits)
        return

    if log_scale:
        pos = vals[vals > 0]
        if pos.size == 0:
            return
        lo, hi = float(np.min(pos)), float(np.max(pos))
        pad = max(0.08 * (np.log10(hi) - np.log10(lo)), 0.05)
        ax.set_ylim(10.0 ** (np.log10(lo) - pad), 10.0 ** (np.log10(hi) + pad))
    else:
        lo, hi = float(np.min(vals)), float(np.max(vals))
        pad = 0.15 * (hi - lo) if hi > lo else max(abs(hi), 1.0) * 0.15
        ax.set_ylim(lo - pad, hi + pad)


def _mark_ylim_outliers(ax, limits: RobustAxisLimits) -> None:
    """Small axis markers when data extend beyond robust limits (xarray colorbar extend)."""
    trans = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
    if limits.has_lower_outliers:
        ax.plot(
            -0.055,
            limits.lo,
            marker="v",
            markersize=4,
            color="0.45",
            transform=trans,
            clip_on=False,
        )
    if limits.has_upper_outliers:
        ax.plot(
            -0.055,
            limits.hi,
            marker="^",
            markersize=4,
            color="0.45",
            transform=trans,
            clip_on=False,
        )


def _axis_unit_label(unit_labels: list[str]) -> str:
    uniq = {u for u in unit_labels if u}
    if len(uniq) == 1:
        return format_unit_for_axis(next(iter(uniq)))
    return "concentration"


def _plot_species_lines(
    ax,
    time_h: np.ndarray,
    species: list[str],
    full_series: list[np.ndarray],
    red_series: list[np.ndarray],
) -> list[np.ndarray]:
    """Plot reduced (behind) then full; return masked series used."""
    plotted: list[np.ndarray] = []
    for sp, yf, yr in zip(species, full_series, red_series):
        pf = _prepare_plot_series(yf)
        pr = _prepare_plot_series(yr)
        plotted.extend([pf, pr])
        if not np.any(np.isfinite(pf)) and not np.any(np.isfinite(pr)):
            continue
        (line,) = ax.plot(time_h, pf, "-", lw=1.5, zorder=3, label=sp)
        ax.plot(
            time_h,
            pr,
            "--",
            lw=1.1,
            alpha=0.8,
            color=line.get_color(),
            dashes=(4, 2),
            zorder=2,
        )
    if species:
        ax.legend(fontsize=6, loc="best", ncol=1)
    return plotted


def _plot_group_panel(
    ax,
    time_h: np.ndarray,
    species: list[str],
    full_series: list[np.ndarray],
    red_series: list[np.ndarray],
    phases: list[str],
    unit_labels: list[str],
    title: str,
    metrics_by_species: dict | None = None,
    plot_robust: bool = True,
    plot_percentile: tuple[float, float] = (PLOT_ROBUST_LOW_PCT, PLOT_ROBUST_HIGH_PCT),
) -> None:
    if metrics_by_species:
        dev_vals = [
            metrics_by_species[sp].devmax2
            for sp in species
            if sp in metrics_by_species and metrics_by_species[sp].tier == "active"
        ]
        if dev_vals:
            title = f"{title}\nmax devmax2={max(dev_vals):.3f}"

    plotted = _plot_species_lines(ax, time_h, species, full_series, red_series)
    ax.set_ylabel(_axis_unit_label(unit_labels))

    use_log = _panel_use_log(plotted)
    if use_log:
        ax.set_yscale("log")
    _apply_axis_limits(
        ax,
        plotted,
        log_scale=use_log,
        robust=plot_robust,
        low_pct=plot_percentile[0],
        high_pct=plot_percentile[1],
    )

    ax.set_title(title, fontsize=8)
    ax.grid(alpha=0.25, which="both")


def make_grouped_panel_figure(
    panels: list[PlotPanel],
    time_h: np.ndarray,
    metrics_by_species: dict,
    suptitle: str,
    out_png: Path,
    ncols: int | None = None,
    plot_robust: bool = True,
    plot_percentile: tuple[float, float] = (PLOT_ROBUST_LOW_PCT, PLOT_ROBUST_HIGH_PCT),
) -> None:
    n = len(panels)
    if n == 0:
        return

    nrows, ncols = balanced_subplot_grid(n, ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows), dpi=160, squeeze=False)
    fig.suptitle(suptitle, fontsize=12)

    for idx in range(nrows * ncols):
        ax = axes[idx // ncols][idx % ncols]
        if idx >= n:
            ax.axis("off")
            continue
        group_name, species, full_series, red_series, phases, unit_labels = panels[idx]
        _plot_group_panel(
            ax, time_h, species, full_series, red_series, phases, unit_labels,
            group_name, metrics_by_species, plot_robust, plot_percentile,
        )
        if idx // ncols == nrows - 1:
            ax.set_xlabel("Time [h]")

    if n > 0:
        fig.legend(
            [
                plt.Line2D([0], [0], color="k", lw=1.5, linestyle="-"),
                plt.Line2D([0], [0], color="k", lw=1.1, linestyle="--"),
            ],
            ["full", "reduced"],
            loc="upper right",
            fontsize=7,
            framealpha=0.9,
        )

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png)
    plt.close(fig)


def build_group_panels(
    groups: list[tuple[str, list[str]]],
    resolved: ResolvedSeries,
) -> list[PlotPanel]:
    """Expand ctrl groups (phase split, max 2 species/panel); drop unplottable."""
    panels: list[PlotPanel] = []
    for group_name, members in expand_groups_to_panels(groups, resolved):
        species: list[str] = []
        full_series: list[np.ndarray] = []
        red_series: list[np.ndarray] = []
        phases: list[str] = []
        unit_labels: list[str] = []
        for sp in members:
            yf, yr, _var, phase, nc_units = resolved[sp]
            species.append(sp)
            full_series.append(yf)
            red_series.append(yr)
            phases.append(phase)
            unit_labels.append(nc_units)
        if species:
            panels.append((group_name, species, full_series, red_series, phases, unit_labels))
    return panels


def make_quicklook_figure(
    groups: list[tuple[str, list[str]]],
    resolved: ResolvedSeries,
    time_h: np.ndarray,
    metrics_by_species: dict,
    label: str,
    out_png: Path,
    ncols: int | None = None,
    plot_robust: bool = True,
    plot_percentile: tuple[float, float] = (PLOT_ROBUST_LOW_PCT, PLOT_ROBUST_HIGH_PCT),
) -> None:
    panels = build_group_panels(groups, resolved)
    make_grouped_panel_figure(
        panels,
        time_h,
        metrics_by_species,
        f"ISSA verification (TARGET_SPC groups): {label}",
        out_png,
        ncols=ncols,
        plot_robust=plot_robust,
        plot_percentile=plot_percentile,
    )


def _safe_group_filename(group_name: str, index: int) -> str:
    slug = re.sub(r"[^\w\-+().=]+", "_", group_name.strip()).strip("_")
    slug = slug.replace("/", "_") or "group"
    return f"{index:02d}_{slug}"


def make_per_group_figures(
    target_groups: list[tuple[str, list[str]]],
    resolved: ResolvedSeries,
    time_h: np.ndarray,
    metrics_by_species: dict,
    label: str,
    out_dir: Path,
    ncols: int | None = None,
    plot_robust: bool = True,
    plot_percentile: tuple[float, float] = (PLOT_ROBUST_LOW_PCT, PLOT_ROBUST_HIGH_PCT),
) -> list[Path]:
    """One PNG per TARGET_SPC comment group (e.g. S(IV+VI), N(V))."""
    group_dir = out_dir / f"{label}_issa_groups"
    written: list[Path] = []
    for idx, (group_name, members) in enumerate(target_groups, start=1):
        panels = build_group_panels([(group_name, members)], resolved)
        if not panels:
            continue
        out_png = group_dir / f"{_safe_group_filename(group_name, idx)}.png"
        make_grouped_panel_figure(
            panels,
            time_h,
            metrics_by_species,
            f"ISSA verification — {group_name} ({label})",
            out_png,
            ncols=ncols,
            plot_robust=plot_robust,
            plot_percentile=plot_percentile,
        )
        written.append(out_png)
    return written


def make_families_figure(
    families: list[list[str]],
    resolved: ResolvedSeries,
    time_h: np.ndarray,
    metrics_by_species: dict,
    label: str,
    out_png: Path,
    ncols: int | None = None,
    plot_robust: bool = True,
    plot_percentile: tuple[float, float] = (PLOT_ROBUST_LOW_PCT, PLOT_ROBUST_HIGH_PCT),
) -> None:
    # One ctrl family → expand by phase / max species per panel (not one crowded subplot)
    family_groups = [(" / ".join(members), members) for members in families]
    panels = build_group_panels(family_groups, resolved)
    make_grouped_panel_figure(
        panels,
        time_h,
        metrics_by_species,
        f"ISSA verification (SPC_FAMILIES): {label}",
        out_png,
        ncols=ncols,
        plot_robust=plot_robust,
        plot_percentile=plot_percentile,
    )


def main() -> int:
    p = argparse.ArgumentParser(description="ISSA full vs reduced NetCDF verification")
    p.add_argument("--full-run", required=True, help="Full-mechanism .run file")
    p.add_argument("--reduced-run", required=True, help="Reduced-mechanism .run file")
    p.add_argument("--full-nc", help="Override full NetCDF path")
    p.add_argument("--reduced-nc", help="Override reduced NetCDF path")
    p.add_argument("--ctrl", help="Override ISSA control file (TARGET_SPC source)")
    p.add_argument("--label", required=True, help="Output label prefix")
    p.add_argument("--out-dir", default="RUN/TESTRUN/diagnostics", help="Output directory")
    p.add_argument(
        "--ncols",
        type=int,
        default=None,
        help="Quicklook columns (default: ~sqrt(N) for square-ish grid)",
    )
    p.add_argument("--eps", type=float, default=1e-12, help="Relative error floor (active tier only)")
    p.add_argument(
        "--min-peak",
        type=float,
        default=0.0,
        help="Min peak concentration for active tier (threshold = max(atol, min-peak))",
    )
    p.add_argument("--top-n", type=int, default=15, help="Rows printed in stdout table")
    p.add_argument(
        "--no-plot-robust",
        action="store_true",
        help="Use min/max y-limits instead of xarray-style robust percentiles",
    )
    p.add_argument(
        "--plot-percentile",
        type=float,
        nargs=2,
        metavar=("LOW", "HIGH"),
        default=(PLOT_ROBUST_LOW_PCT, PLOT_ROBUST_HIGH_PCT),
        help="Percentiles for robust y-limits (default: 2 98, like xarray robust=True)",
    )
    args = p.parse_args()

    cminor_dir = Path(__file__).parent.parent.resolve()
    full_run_path = _resolve_path(cminor_dir, args.full_run)
    red_run_path = _resolve_path(cminor_dir, args.reduced_run)
    if full_run_path is None or red_run_path is None:
        print("ERROR: invalid run file paths", file=sys.stderr)
        return 2

    full_run = parse_run_file(full_run_path)
    red_run = parse_run_file(red_run_path)

    ctrl_path = _resolve_path(cminor_dir, args.ctrl or full_run.get("RedCtrlFile"))
    if ctrl_path is None or not ctrl_path.is_file():
        print(f"ERROR: RedCtrlFile not found: {ctrl_path}", file=sys.stderr)
        return 2

    full_nc = _resolve_path(cminor_dir, args.full_nc or full_run.get("NetCdfFile"))
    red_nc = _resolve_path(cminor_dir, args.reduced_nc or red_run.get("NetCdfFile"))
    if full_nc is None or not full_nc.is_file():
        print(f"ERROR: full NetCDF not found: {full_nc}", file=sys.stderr)
        return 2
    if red_nc is None or not red_nc.is_file():
        print(f"ERROR: reduced NetCDF not found: {red_nc}", file=sys.stderr)
        return 2

    out_dir = _resolve_path(cminor_dir, args.out_dir) or Path(args.out_dir)
    species_list = parse_target_species(ctrl_path)
    target_groups = parse_target_species_groups(ctrl_path)
    species_families = parse_species_families(ctrl_path)

    parity_warnings = check_run_parity(full_run, red_run)
    if parity_warnings:
        print("[issa-verify] WARNING: full vs reduced run config differs:")
        for w in parity_warnings:
            print(w)

    print(f"[issa-verify] full NC:    {full_nc}")
    print(f"[issa-verify] reduced NC: {red_nc}")
    print(f"[issa-verify] ctrl:       {ctrl_path}")
    print(f"[issa-verify] TARGET_SPC: {len(species_list)} species in {len(target_groups)} groups")
    if species_families:
        print(f"[issa-verify] SPC_FAMILIES: {len(species_families)} families")

    _atol_gas, _atol_aqua, atol = run_comparison_atol(full_run, red_run)
    peak_thresh = comparison_threshold(atol, args.min_peak)
    print(f"[issa-verify] peak threshold: {peak_thresh:.4e}  (atol={atol:.4e}, min-peak={args.min_peak:.4e})")

    metrics, skipped, _time_h = compare_datasets(
        full_nc,
        red_nc,
        species_list,
        eps=args.eps,
        peak_threshold=peak_thresh,
        atol=atol,
    )

    if skipped:
        print(f"[issa-verify] skipped {len(skipped)} species (not in both NC files):")
        for sp in skipped[:20]:
            print(f"  - {sp}")
        if len(skipped) > 20:
            print(f"  ... and {len(skipped) - 20} more")

    if not metrics:
        print("ERROR: no species resolved for comparison", file=sys.stderr)
        return 1

    # Log variable mappings
    for m in metrics:
        print(f"  {m.species} -> {m.nc_var}  [{m.nc_units}, {m.phase}]")

    result = ComparisonResult(
        label=args.label,
        full_nc=str(full_nc),
        reduced_nc=str(red_nc),
        ctrl_file=str(ctrl_path),
        peak_threshold=peak_thresh,
        atol=atol,
        resolved=metrics,
        skipped=skipped,
        warnings=parity_warnings,
    )

    csv_path = out_dir / f"{args.label}_issa_summary.csv"
    json_path = out_dir / f"{args.label}_issa_summary.json"
    png_path = out_dir / f"{args.label}_issa_quicklook.png"
    families_png_path = out_dir / f"{args.label}_issa_families.png"

    write_summary_csv(csv_path, metrics)
    write_summary_json(json_path, result)
    print_summary_table(metrics, top_n=args.top_n)

    plot_species = sorted({sp for group in target_groups for sp in group[1]} | {sp for fam in species_families for sp in fam})
    _sp_plot, time_h, resolved = load_species_series(full_nc, red_nc, plot_species)
    metrics_map = {m.species: m for m in metrics}
    plot_robust = not args.no_plot_robust
    low_pct, high_pct = args.plot_percentile
    if not (0.0 <= low_pct < high_pct <= 100.0):
        print("ERROR: --plot-percentile LOW HIGH must satisfy 0 <= LOW < HIGH <= 100", file=sys.stderr)
        return 2
    plot_pct = (low_pct, high_pct)
    make_quicklook_figure(
        target_groups,
        resolved,
        time_h,
        metrics_map,
        args.label,
        png_path,
        ncols=args.ncols,
        plot_robust=plot_robust,
        plot_percentile=plot_pct,
    )
    group_pngs = make_per_group_figures(
        target_groups,
        resolved,
        time_h,
        metrics_map,
        args.label,
        out_dir,
        ncols=args.ncols,
        plot_robust=plot_robust,
        plot_percentile=plot_pct,
    )
    if species_families:
        make_families_figure(
            species_families,
            resolved,
            time_h,
            metrics_map,
            args.label,
            families_png_path,
            ncols=args.ncols,
            plot_robust=plot_robust,
            plot_percentile=plot_pct,
        )

    print(f"\n[issa-verify] wrote: {csv_path}")
    print(f"[issa-verify] wrote: {json_path}")
    print(f"[issa-verify] wrote: {png_path}")
    if group_pngs:
        print(f"[issa-verify] wrote {len(group_pngs)} group figures under: {group_pngs[0].parent}")
    if species_families:
        print(f"[issa-verify] wrote: {families_png_path}")
    print(f"[issa-verify] resolved {len(metrics)} / {len(species_list)} target species")
    return 0


if __name__ == "__main__":
    sys.exit(main())
