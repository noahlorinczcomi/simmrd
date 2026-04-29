# simmrd CLI

A command-line interface for running [`simmrd`](https://github.com/noahlorinczcomi/simmrd) simulations from a YAML parameter file, with built-in support for multi-iteration Monte Carlo runs.

The CLI is managed by [pixi](https://pixi.sh), which handles the R and Python environment automatically.

---

## Setup (once)

> [!NOTE]
> Install [pixi](https://pixi.sh/latest/#installation) if you don't already have it: `curl -fsSL https://pixi.sh/install.sh | sh`

```bash
cd cli/
pixi install       # resolve and download the R + Python environment
pixi run setup     # install simmrd from the parent source tree
```

---

## Usage

```
python simmrd_cli.py --params FILE --output FILE [--iterations INT] [--seed INT]
```

Or through pixi:

```bash
pixi run simulate --params FILE --output FILE [--iterations INT] [--seed INT]
```

### Flags

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--params` | path | required | YAML file of simulation parameters |
| `--output` | path | required | Output file (saved as `.rds`) |
| `--iterations` | int | `1` | Number of simulation replicates |
| `--seed` | int | none | Master random seed for reproducibility |

---

## Examples

### Single replicate

```bash
pixi run simulate \
  --params params/example_summary.yaml \
  --output results.rds
```

### 500 replicates, reproducible

```bash
pixi run simulate \
  --params params/example_summary.yaml \
  --output results.rds \
  --iterations 500 \
  --seed 42
```

### Individual-level data, 100 replicates

```bash
pixi run simulate \
  --params params/example_individual.yaml \
  --output results_individual.rds \
  --iterations 100 \
  --seed 1
```

---

## Parameter files

YAML files map directly to [`set_params()`](https://github.com/noahlorinczcomi/simmrd) arguments.
Two annotated examples are in `params/`:

| File | Generator called |
|------|-----------------|
| `params/example_summary.yaml` | `generate_summary()` |
| `params/example_individual.yaml` | `generate_individual()` |

The `type` field at the top of the YAML controls which generator is used (`summary` or `individual`). All other keys are passed directly to `set_params()`. Unrecognised keys are ignored with a warning.

### Minimal summary-data YAML

```yaml
type: summary
number_of_exposures:        2
true_causal_effects:        [0.3, 0.1]
prop_gwas_overlap_Xs_and_Y: 0
number_of_CHP_causal_SNPs:  20
ratio_of_CHP_variance:      0.25
CHP_correlation:            -0.5
```

---

## Output

The output is an R list saved as an `.rds` file.

- **Single replicate** (`--iterations 1`): a named list with elements `bx`, `bxse`, `by`, `byse`, `RhoME`, `LDMatrix`, `LDhatMatrix`, `theta`, `IVtype`, etc.
- **Multiple replicates** (`--iterations > 1`): a length-`n` list where each element is one replicate's named list.

Read the output back in R:

```r
results <- readRDS("results.rds")

# single replicate
head(results$bx)

# multiple replicates — extract by from each
by_list <- lapply(results, `[[`, "by")
```

> **Format TBD** — the `.rds` format preserves all output without data loss while the final output format is decided.

---

## Reproducibility

| Scenario | Behaviour |
|----------|-----------|
| `--seed` given | Master seed set once; per-iteration seeds derived deterministically. Identical output on every re-run. |
| `--seed` absent | Fresh random seed sequence each run. |

The per-iteration seeding strategy prevents `generate_individual()`'s internal `set.seed(1)` default from making every replicate identical when looping.
