# The First Kalman Filter to Fly — Apollo free-return LKF reconstruction

A learning project: the linearized Kalman filter of NASA TR R-135 (Smith,
Schmidt & McGee, 1962) — the first practical application of the Kalman
filter — reimplemented in MATLAB and flown over a complete Earth → Moon →
Earth ballistic free-return mission. Companion report:
`Apollo_LKF_Report.pdf` (the full story, mathematics, design decisions and
original NASA sources).

<img width="800" height="450" alt="Apollo_FreeReturn_LKF-ezgif com-video-to-gif-converter" src="https://github.com/user-attachments/assets/0d48449c-3fa5-4c8d-a2f2-8a116486d8dc" />

Author: Sreeram Anil.

## Files

| File | Purpose |
|---|---|
| `apollo_lkf_mission.m`   | Everything: free-return trajectory targeting, truth simulation, R-135 optical measurements, both filter variants (LKF / estimate-linearized), all six figures. Saves `apollo_lkf_results.mat`. |
| `apollo_lkf_animation.m` | Renders the mission animation (MP4 via VideoWriter) from the saved results. |
| `Apollo_LKF_Report.pdf`  | The report. |
| `Apollo_FreeReturn_LKF.mp4` | Pre-rendered animation. |

## Run (MATLAB, no toolboxes required)

```matlab
R = apollo_lkf_mission();     % ~1–3 min; prints targeting iterations and
                              % filter results; writes fig1..fig6 PNGs
apollo_lkf_animation();       % writes Apollo_FreeReturn_LKF.mp4 (~1–5 min)
```

Notes
- Trajectory targeting warm-starts from the recorded corridor solution;
  set `cfg.do_scan = true` inside `apollo_lkf_mission.m` to redo the full
  corridor scan from scratch (~2 min extra).
- Every model constant and design decision is annotated in the source with
  the NASA document it came from (R-135 appendices/equations, Bellcomm
  TR-66-310-4, Apollo 13 Mission Report Suppl. 1, NASA TM-86847, MIT R-700).
- Also runs under GNU Octave ≥ 8 (used for headless validation); without a
  video backend the animation falls back to PNG frames plus an ffmpeg
  one-liner it prints.

## Headline result

Uncorrected, the true trajectory (R-135's own sampled injection error,
scaled, plus a venting-class disturbance) misses Earth entry by ~31,000 km.
On identical sextant-grade angle sightings (σ = 10 arcsec), the filter
linearized about the nominal diverges after the lunar flyby — while the
estimate-linearized variant (the modification R-135 p. 15 already
recommends, later named the extended Kalman filter) tracks the whole
5.72-day mission and arrives at entry with ~1.4 km error, having also
identified the unmodeled venting acceleration and the 8-km Earth-horizon
bias taken from Apollo 13's actual P23 sightings.
