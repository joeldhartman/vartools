"""PerLC: per-light-curve parameter values for batch pipeline runs."""


class PerLC:
    """Wrap a sequence of per-light-curve values for use as a command parameter.

    Pass a PerLC as a command parameter value to supply a different scalar for
    each light curve in a batch run.  numpy arrays and pd.Series are also
    accepted directly by commands and auto-wrapped.  Plain Python lists require
    explicit PerLC([...]) to avoid ambiguity with fixed multi-valued parameters
    such as ld_coeffs.

    Examples
    --------
    Explicit:
        cmd.LS(minper=PerLC([0.1, 0.2, 0.15]), maxper=10.0)

    Via numpy (auto-detected):
        cmd.LS(minper=np.array([0.1, 0.2, 0.15]), maxper=10.0)

    Via DataFrame column (auto-detected):
        cmd.LS(minper=catalog_df["min_period"], maxper=10.0)
    """

    def __init__(self, values):
        import numpy as np
        try:
            import pandas as pd
            if isinstance(values, pd.Series):
                values = values.to_numpy()
        except ImportError:
            pass
        self._values = np.asarray(values, dtype=float)
        if self._values.ndim != 1:
            raise ValueError("PerLC values must be 1-dimensional")

    @property
    def values(self):
        return self._values

    def __len__(self):
        return len(self._values)

    def __getitem__(self, i):
        return float(self._values[i])

    def __repr__(self):
        return f"PerLC(n={len(self._values)})"
