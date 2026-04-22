# Miscellaneous

Commands that don't fit cleanly into any of the other command categories.

---

## `-findblends`

```
-findblends
    matchrad ["radec"] ["xycol" xcol ycol]
    < "fix" period | "list" ["column" col] | "fixcolumn" <colname | colnum> >
    ["starlist" starlistfile] ["zeromag" zeromagval] ["nofluxconvert"]
    ["Nharm" Nharm] ["omatches" outputmatchfile]
```

Determine whether a detected periodic signal is likely due to contamination (blending) from a nearby variable star. For each potential variable, the routine measures the flux amplitude of all nearby light curves and reports the one with the highest amplitude.

A light curve list (`-l`) is required with x and y coordinates as additional columns.

**Parameters**

- `matchrad` — Matching radius for identifying potentially blended stars. In arcseconds if `"radec"` is given; in pixel units otherwise.
- `"radec"` — Treat x and y coordinates as RA and Dec (in degrees) and use `matchrad` in arcseconds.
- `"xycol" xcol ycol` — Columns in the input list for x and y coordinates (default: next two available columns).
- Period source: `"fix" period`, `"list"`, or `"fixcolumn"`.
- `"starlist" starlistfile` — Match the input list against this external catalog instead of itself. Format: `lcname x y`.
- `"zeromag" zeromagval` — Zero-point magnitude for converting magnitudes to fluxes (default: 25.0).
- `"nofluxconvert"` — Skip the magnitude-to-flux conversion (use if input is already in flux units).
- `"Nharm" Nharm` — Number of harmonics for the Fourier series amplitude measurement (default: 2; use 0 for a pure sinusoid).
- `"omatches" outputmatchfile` — Output names and flux amplitudes of all stars matching each target to this file.

---

