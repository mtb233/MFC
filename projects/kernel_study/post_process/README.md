# Bubble evolution plotter

`plot_bubble_evolution.py` plots the time history of a single Lagrangian bubble from one
or more MFC cases. A JSON config lists the figures, the bubble property on each one, and
which cases supply curves. MFC writes its output without units; the script converts it
back to physical units before plotting.

Run it from the MFC root directory:

```bash
python projects/kernel_study/post_process/plot_bubble_evolution.py \
    projects/kernel_study/post_process/bubble_plots.example.json
```

Requires `jsonschema`, `pandas`, and `matplotlib`.

## Input data

With `lag_params%write_bubbles = 'T'`, every MPI rank writes `D/lag_bubble_evol_<rank>.dat`
in the case directory. Each file has one header line and whitespace-separated columns:

| Column | Meaning | Dimension |
|---|---|---|
| `currentTime` | simulation time | time |
| `particleID` | bubble ID | — |
| `x`, `y`, `z` | bubble position | length |
| `coreVaporMass` | vapor mass in the bubble core | mass |
| `coreVaporConcentration` | vapor mass fraction in the core | dimensionless |
| `radius` | bubble radius | length |
| `interfaceVelocity` | dR/dt | velocity |
| `corePressure` | pressure inside the bubble | pressure |

**This tool only supports cases with one bubble.** Ranks that don't own the bubble write
header-only files. For each case, the script uses the one file with data rows. It stops
with an error in these situations:

- no file has data (the case hasn't run, or `write_bubbles` is off)
- more than one file has data (for example, the bubble moved between ranks)
- the chosen file contains more than one `particleID`

## Converting to physical units

MFC data has no units. It is scaled by reference values that default to those in
`projects/kernel_study/generate_cases.py`:

| Reference | Default | Meaning |
|---|---|---|
| `x_0` | `1.0e-3` m | length scale |
| `rho_0` | `1.0e3` kg/m³ | density scale |
| `c_0` | `1475.0` m/s | velocity scale |

These give the conversion factor (in SI units) for each dimension:

| Dimension | Scale |
|---|---|
| time | `x_0 / c_0` |
| length | `x_0` |
| velocity | `c_0` |
| pressure | `rho_0 * c_0**2` |
| mass | `rho_0 * x_0**3` |
| dimensionless | 1 |

If a case set uses different reference values, override them in the config's `reference`
block. The script can't import them from `generate_cases.py`, because that file generates
cases as soon as it is imported. Keep the two in sync by hand.

### Display units

Each axis takes a `unit` that must fit its property's dimension:

| Dimension | Units (the first is the default) |
|---|---|
| time | `us`, `s`, `ms`, `ns` |
| length | `um`, `m`, `mm` |
| velocity | `m/s` |
| pressure | `kPa`, `Pa`, `MPa`, `atm` |
| mass | `ng`, `kg`, `pg` |
| dimensionless | `-` |

Any property can use `"unit": "normalized"` instead. This divides the curve by its first
value (for example, `radius` → R/R₀). A unit that doesn't fit the property, such as `kPa`
for `radius`, stops the script with an error that lists the valid choices.

## Config file

The config is validated against `bubble_plots.schema.json`. Unknown keys are rejected,
so a typo such as `ylimit` gets an error instead of being silently ignored. Relative
paths (`output_dir`, `case`) are resolved from the current working directory, which
should be the MFC root (e.g. `projects/kernel_study/cases/mc_large`).

### Top level

| Key | Required | Default | Description |
|---|---|---|---|
| `plots` | yes | — | list of plot objects |
| `reference` | no | see above | `{x_0, rho_0, c_0}` overrides (any subset) |
| `output_dir` | no | `"projects/kernel_study/results"` | where figures are written |
| `format` | no | `"png"` | `png`, `pdf`, `svg`, or `eps` |
| `dpi` | no | `300` | resolution of the saved figure |

### Plot

| Key | Required | Default | Description |
|---|---|---|---|
| `name` | yes | — | output file name without extension (`[A-Za-z0-9_.-]`) |
| `y` | yes | — | axis object for the y-axis |
| `curves` | yes | — | list of curve objects |
| `x` | no | `currentTime` in `us` | axis object for the x-axis |
| `title` | no | none | figure title |
| `figsize` | no | `[5.5, 4.0]` | figure size in inches |
| `xlim`, `ylim` | no | automatic | `[min, max]` in display units |
| `legend` | no | `true` | show the legend |

### Axis

| Key | Required | Default | Description |
|---|---|---|---|
| `property` | yes | — | one of the column names above, except `particleID` |
| `unit` | no | the dimension's default | display unit or `"normalized"` |
| `label` | no | generated | axis label (matplotlib mathtext is allowed) |
| `scale` | no | `"linear"` | `linear` or `log` |

Generated labels look like `Radius (μm)` or, for normalized axes, `$R/R_0$`.

### Curve

| Key | Required | Default | Description |
|---|---|---|---|
| `case` | yes | — | case directory (the one that contains `D/`) |
| `label` | no | case directory name | legend entry |
| `style` | no | `{}` | keyword arguments passed to matplotlib's `Axes.plot` (`color`, `linestyle`, `linewidth`, `marker`, ...) |

## Example

`bubble_plots.example.json` draws two figures for `mc_large`: R/R₀ and core pressure in
kPa. Below is a config that compares kernels on the same grid:

```json
{
  "output_dir": "projects/kernel_study/results",
  "plots": [
    {
      "name": "radius_large",
      "y": {"property": "radius", "unit": "normalized"},
      "ylim": [0.5, 2.0],
      "curves": [
        {"case": "projects/kernel_study/cases/df_large", "label": "Default"},
        {"case": "projects/kernel_study/cases/mc_large", "label": "MC"},
        {"case": "projects/kernel_study/cases/k1_large", "label": "K1", "style": {"linestyle": "--"}},
        {"case": "projects/kernel_study/cases/k2_large", "label": "K2", "style": {"linestyle": ":"}}
      ]
    },
    {
      "name": "wall_velocity_large",
      "x": {"property": "currentTime", "unit": "us"},
      "y": {"property": "interfaceVelocity", "unit": "m/s", "label": "$\\dot{R}$ (m/s)"},
      "curves": [
        {"case": "projects/kernel_study/cases/mc_large", "label": "MC"}
      ]
    }
  ]
}
```

A case that is used in several plots is only read once.

## Other files in this directory

- `plot_comparison.py`: the older single-purpose script (R/R₀ against golden and analytical data).
- `D_golden/`: reference output from an 8-rank run (only rank 7 holds the bubble).
- `analytical_solution.csv`: reference R/R₀ curve, with time in units of acoustic periods.
