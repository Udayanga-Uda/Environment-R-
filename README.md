# Residual emissions in mitigation scenarios

An R workflow for comparing scenario data from several modelling teams, focused on residual gross emissions: the emissions that remain when net emissions reach zero, and which sectors they come from.

**The data are simulated.** Three fictitious modelling teams (Model A, B, C) submit gross emissions by sector and carbon dioxide removal (CDR) for three scenarios (1.5C, 2C, 2C delayed) from 2020 to 2100. The numbers are illustrative and are not results from real integrated assessment models.

## What the script does

1. Simulates the submissions and stores them in a SQLite database (sectors, emissions, CDR, milestones).
2. Reads the data back with SQL joins.
3. Checks the submissions: missing rows, unit errors (one series is reported in Gt instead of Mt on purpose) and negative values. Problems are logged, fixed and saved to a QA table.
4. Computes gross, CDR and net emissions, the year net emissions reach zero, and residual gross emissions at that point and in 2050 and 2100.
5. Compares the three models (median and range) and looks at the agriculture share of what remains.
6. Keeps a small deliverable tracker with status and overdue flags.
7. Saves figures and CSV tables.

## How to run

You need R 4.1 or newer (RStudio is optional).

```r
source("analysis.R")
```

Missing packages are installed automatically: tidyverse, DBI, RSQLite, patchwork. The seed is fixed, so results are reproducible.

## Output

- `data/scenarios.sqlite`: the simulated database
- `outputs/figures/`: emissions pathways, residual emissions by sector group, agriculture share, net-zero years, and a combined panel
- `outputs/tables/`: QA log, totals, net-zero years and residuals, cross-model spread, agriculture share, tracker

## Limitations

- The scenario shapes come from simple decay curves with a floor for each sector, not from an energy-economy or land-use model.
- Emissions are aggregated into CO2-equivalents, so the choice of metric is not explored.
- CDR is a smooth curve split into three types, without feasibility or sustainability limits.
- The tracker entries and dates are made up to show the structure.

## Next steps

Replace the simulation (section 2 of the script) with real scenario data, for example the IPCC AR6 Scenario Database, which is available through the IIASA Scenario Explorer and the `pyam` / `ixmp4` tools, and keep the checks and figures.

## Author

Prabhath Samarahewa
