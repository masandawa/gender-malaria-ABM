# gender-malaria-ABM

MATLAB code for the manuscript *A calibrated agent-based model of how sex-specific biology, gendered behaviour and bed-net access shape malaria disparities and net prioritisation in Tanzania* (Masandawa et al.).

The code reproduces every number, table and figure in the paper. It runs in MATLAB R2019b or later without toolboxes, and in GNU Octave 6 or later.

## Quick start

```matlab
% from the repository folder
run_all            % MODE = 'paper' by default; set MODE = 'quick' at the top for a fast check
```

| Mode | Purpose | Approximate run time |
|---|---|---|
| `quick` | Small populations and few realisations; checks that the whole pipeline runs and produces all CSVs and figures | about 6 min in Octave (one core); faster in MATLAB |
| `paper` | Settings used for the published results | about 3 h in Octave on one core; faster in MATLAB |

The experiments are independent once calibration is done. To save time, set any `RUN.Ax` flag in `run_all.m` to `false` to skip an experiment whose results already exist, or run experiments in separate MATLAB sessions (see below).

To redraw tables and figures from existing results without re-running any simulation:

```matlab
addpath model plotting config analyses
make_tables; make_figures
```

## Repository layout

```
run_all.m                 single entry point (choose MODE, switch experiments on/off)
config/settings.m         run sizes for 'quick' and 'paper'
model/                    the agent-based model
  abm_params.m            all parameters with sources (Table 2 of the paper)
  abm_init.m              initial population and mosquitoes
  abm_run.m               daily simulator (demography, pregnancy, hourly exposure,
                          ITN use, infection, treatment, attribution, mosquitoes)
  activity_rows.m         hourly probability of being outdoors (gender contrast)
  allocate.m              ITN allocation and selection rules
  scenario.m              mechanism blocks [B A G S] (behaviour, programme access,
                          pregnancy biology, sex-specific clearance)
  summarize_run.m         outcome measures
  concentration.m         event concentration and clinical-transmission overlap
  rbinom.m rpois.m set_seed.m base_params.m addrow.m struct2csv.m csv2struct.m
analyses/
  A0_calibration.m        single-target calibration (Morogoro, TDHS-MIS 2022) and profile
  A1_behaviour.m          female-male activity contrast; event concentration
  A2_factorial.m          2^4 factorial experiment (Shapley decomposition in make_tables)
  A3_pregnancy_sensitivity.m
  A4_itn_use.m            habitual vs nightly ITN use; standardised infection risk
  A5_allocation.m         fixed-resource ITN recipient selection (7 rules)
  A6_biting_ratio.m       outdoor:indoor biting ratio, fixed and re-calibrated density
  A7_uncertainty.m        Latin hypercube joint parameter uncertainty
  A8_verification.m       verification tests
  A9_mechanism_checks.m   clearance pathway, alternative immunity structures, ITN protection components
plotting/
  make_tables.m           summary CSVs used in Tables 3-5 and the text
  make_figures.m          Figures 2-6 (17 cm wide, 600-dpi PNG and vector PDF)
  fig_new.m fig_save.m panel.m mlegend.m eb.m okabe.m
results/                  CSV outputs (created by the run)
figures/                  figures (created by the run)
```

## Which file produces which output

| Paper item | Produced by | Result files |
|---|---|---|
| Calibration, profile (Fig. 2a-b) | `A0_calibration.m` | `A0_*.csv` |
| Table 3, Fig. 2c-d (empirical evaluation) | `A2_factorial.m` (reference configuration) | `summary_evaluation.csv` |
| Fig. 3 (activity, concentration) | `A1_behaviour.m` | `A1_*.csv`, `summary_A1*.csv` |
| Table 4, Fig. 4 (mechanisms, pregnancy) | `A2_factorial.m`, `A3_pregnancy_sensitivity.m` | `summary_A2_*.csv`, `summary_A3.csv` |
| Fig. 5a-b (ITN protection) | `A2_factorial.m`, `A4_itn_use.m` | `summary_A4.csv` |
| Table 5, Fig. 5c-d (selection rules) | `A5_allocation.m` | `summary_A5.csv` |
| Fig. 6 (robustness) | `A6_biting_ratio.m`, `A7_uncertainty.m` | `summary_A6.csv`, `summary_A7.csv` |
| Verification tests | `A8_verification.m` | `A8_verification.txt` |
| Section 3.4 and Table 4 lower part (clearance pathway, immunity structures); Section 3.5 (randomised ITN protection components) | `A9_mechanism_checks.m` | `summary_A9_immunity.csv`, `summary_A9_itn.csv` |

## Running experiments in parallel

After calibration (`A0_calibration`), each experiment needs only the settings `X` and the calibrated ratio `m`. Use separate MATLAB sessions to run them at the same time:

```matlab
addpath model analyses config plotting
X = settings('paper'); m = 19.30;     % m* from results/A0_calibration_summary.csv
A5_allocation(X, m)                   % e.g. in one session
A7_uncertainty(X, m, 1:20)            % the Latin hypercube can be split across sessions;
A7_uncertainty(X, m, 21:40)           % make_tables merges the parts automatically
```

## Reproducibility

Every realisation uses a fixed seed (`set_seed.m`). Scenarios compared within a realisation share a common burn-in and common random numbers. MATLAB and Octave use different random-number streams, so individual realisations differ between them, but summary results agree within Monte Carlo error.

## Licence

MIT (add a LICENSE file before release if required by your institution).
