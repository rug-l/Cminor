"""
Create diagnostic PNGs comparing Cminor test vs reference NetCDF outputs.

Usage:
  python PYTHONSCRIPTS/plot_test_diagnosis.py \
    --test RUN/TESTRUN/LLNL_nHeptane/LLNL_nHeptane_test.nc \
    --ref  RUN/TESTRUN/LLNL_nHeptane/LLNL_nHeptane_reference.nc \
    --mode comb \
    --label LLNL_nHeptane \
    --out RUN/TESTRUN/diagnostics/LLNL_nHeptane_diagnostic.png
"""

import argparse
import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt


IGNORE_EXACT = {
    "trajectory",
    "Step_Size",
}
IGNORE_PREFIX = (
    "Number_",
    "dryRadius_",
    "wetRadius_",
    "pH_Value_",
)


def safe_rel_err(test_arr, ref_arr, eps):
    denom = np.maximum(np.abs(ref_arr), eps)
    return np.abs(test_arr - ref_arr) / denom


def is_timeseries(var, n_t):
    # Keep only 1D vars aligned with time axis
    try:
        return len(var.shape) == 1 and var.shape[0] == n_t
    except Exception:
        return False


def should_ignore_key(key):
    if key in IGNORE_EXACT:
        return True
    for pref in IGNORE_PREFIX:
        if key.startswith(pref):
            return True
    return False


def load_time(ds, mode):
    if "time" not in ds.variables:
        raise ValueError("NetCDF missing 'time' variable")

    t = np.array(ds["time"][:], dtype=float)
    if mode == "comb":
        # keep current project convention: combustion plots in ms
        return t * 1000.0, "Time [ms]"
    return t, "Time [h]"


def top_species(test_ds, ref_ds, n_t, top_n=5):
    common = set(test_ds.variables.keys()) & set(ref_ds.variables.keys())

    scored = []
    for key in common:
        if should_ignore_key(key):
            continue

        v_test = test_ds.variables[key]
        v_ref = ref_ds.variables[key]

        if not is_timeseries(v_test, n_t) or not is_timeseries(v_ref, n_t):
            continue

        a = np.array(v_test[:], dtype=float)
        b = np.array(v_ref[:], dtype=float)

        if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
            continue

        abs_max = float(np.max(np.abs(a - b)))
        scored.append((key, abs_max))

    scored.sort(key=lambda x: x[1], reverse=True)
    return [k for k, _ in scored[:top_n]]


def make_plot(test_nc, ref_nc, mode, label, out_png, eps, top_n):
    os.makedirs(os.path.dirname(out_png), exist_ok=True)

    with nc.Dataset(test_nc) as ds_test, nc.Dataset(ref_nc) as ds_ref:
        time, time_label = load_time(ds_test, mode)
        n_t = time.size

        if "Temperature" in ds_test.variables and "Temperature" in ds_ref.variables:
            has_temp = (
                is_timeseries(ds_test.variables["Temperature"], n_t)
                and is_timeseries(ds_ref.variables["Temperature"], n_t)
            )
        else:
            has_temp = False

        keys = top_species(ds_test, ds_ref, n_t, top_n=top_n)
        if "Temperature" in keys:
            keys = [k for k in keys if k != "Temperature"]

        # ---------- Figure ----------
        fig = plt.figure(figsize=(12, 10), dpi=180)

        # Row 1: Temperature + DeltaT
        ax1 = fig.add_subplot(3, 1, 1)
        if has_temp:
            t_test = np.array(ds_test["Temperature"][:], dtype=float)
            t_ref = np.array(ds_ref["Temperature"][:], dtype=float)
            dT = t_test - t_ref

            ax1.plot(time, t_ref, label="Temperature ref", lw=2)
            ax1.plot(time, t_test, label="Temperature test", lw=1.5, alpha=0.85)
            ax1.set_ylabel("Temperature [K]")
            ax1.grid(alpha=0.3)

            ax1b = ax1.twinx()
            ax1b.plot(time, dT, "--", label="DeltaT", alpha=0.8)
            ax1b.set_ylabel("DeltaT [K]")

            lines1, labels1 = ax1.get_legend_handles_labels()
            lines2, labels2 = ax1b.get_legend_handles_labels()
            ax1.legend(lines1 + lines2, labels1 + labels2, loc="best")
        else:
            ax1.text(0.02, 0.5, "Temperature variable not available", transform=ax1.transAxes)
            ax1.set_ylabel("Temperature [K]")
            ax1.grid(alpha=0.3)

        ax1.set_title(f"Cminor regression diagnostics: {label}")
        ax1.set_xlabel(time_label)

        # Row 2: top absolute-error species trajectories
        ax2 = fig.add_subplot(3, 1, 2)
        if len(keys) == 0:
            ax2.text(0.02, 0.5, "No comparable 1D species found", transform=ax2.transAxes)
        else:
            for key in keys:
                a = np.array(ds_test[key][:], dtype=float)
                b = np.array(ds_ref[key][:], dtype=float)
                ax2.plot(time, b, lw=2, label=f"{key} ref")
                ax2.plot(time, a, "--", lw=1.2, alpha=0.9, label=f"{key} test")
        ax2.set_ylabel("Value")
        ax2.set_xlabel(time_label)
        ax2.grid(alpha=0.3)
        ax2.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8)

        # Row 3: relative error of top species
        ax3 = fig.add_subplot(3, 1, 3)
        if len(keys) > 0:
            for key in keys:
                a = np.array(ds_test[key][:], dtype=float)
                b = np.array(ds_ref[key][:], dtype=float)
                rel = safe_rel_err(a, b, eps=eps)
                ax3.plot(time, rel, lw=1.5, label=key)
        ax3.set_yscale("log")
        ax3.set_ylabel("Relative error")
        ax3.set_xlabel(time_label)
        ax3.grid(alpha=0.3, which="both")
        ax3.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8)

        plt.tight_layout()
        fig.savefig(out_png)
        plt.close(fig)

        # Metrics for console
        metrics = {}
        if has_temp:
            t_test = np.array(ds_test["Temperature"][:], dtype=float)
            t_ref = np.array(ds_ref["Temperature"][:], dtype=float)
            metrics["temp_abs_max"] = float(np.max(np.abs(t_test - t_ref)))

        key_metrics = {}
        for key in keys:
            a = np.array(ds_test[key][:], dtype=float)
            b = np.array(ds_ref[key][:], dtype=float)
            rel = safe_rel_err(a, b, eps=eps)
            key_metrics[key] = {
                "abs_max": float(np.max(np.abs(a - b))),
                "rel_p95": float(np.percentile(rel, 95)),
                "rel_max": float(np.max(rel)),
            }
        metrics["species"] = key_metrics

    return metrics


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--test", required=True, help="Path to *_test.nc")
    p.add_argument("--ref", required=True, help="Path to *_reference.nc")
    p.add_argument("--mode", required=True, choices=["atm", "comb"], help="Comparison mode")
    p.add_argument("--label", required=True, help="Label for figure title")
    p.add_argument("--out", required=True, help="Output PNG path")
    p.add_argument("--eps", type=float, default=1e-12, help="Relative error floor")
    p.add_argument("--top-n", type=int, default=5, help="Top species count")
    args = p.parse_args()

    metrics = make_plot(
        test_nc=args.test,
        ref_nc=args.ref,
        mode=args.mode,
        label=args.label,
        out_png=args.out,
        eps=args.eps,
        top_n=args.top_n,
    )

    print(f"[diagnostics] wrote: {args.out}")
    if "temp_abs_max" in metrics:
        print(f"[diagnostics] temp_abs_max: {metrics['temp_abs_max']:.6e}")
    if metrics.get("species"):
        print("[diagnostics] top species metrics:")
        for k, v in metrics["species"].items():
            print(
                f"  {k}: abs_max={v['abs_max']:.6e}, "
                f"rel_p95={v['rel_p95']:.6e}, rel_max={v['rel_max']:.6e}"
            )


if __name__ == "__main__":
    main() 