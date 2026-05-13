"""PerLC: per-light-curve parameter values for batch pipeline runs."""


class PerLC:
    """Wrap a sequence of per-light-curve values for use as a command parameter.

    Pass a PerLC as a command parameter value to supply a different value
    for each light curve in a batch run.  numpy arrays and pd.Series are
    also accepted directly by commands and auto-wrapped.  Plain Python
    lists require explicit PerLC([...]) to avoid ambiguity with fixed
    multi-valued parameters such as ld_coeffs.

    Numeric values, strings, and any other 1-D iterable are all accepted;
    the wrapper preserves element type rather than forcing a numeric cast,
    so PerLC(["a", "b", "c"]) can drive string-valued parameters such as
    cmd.o(outname=PerLC([...])).

    Examples
    --------
    Numeric, explicit:
        cmd.LS(minper=PerLC([0.1, 0.2, 0.15]), maxper=10.0)

    Numeric, via numpy (auto-detected):
        cmd.LS(minper=np.array([0.1, 0.2, 0.15]), maxper=10.0)

    String-valued (per-LC output filename basenames):
        cmd.o(outdir="results", outname=PerLC(["a", "b", "c"]))
    """

    def __init__(self, values):
        import numpy as np
        try:
            import pandas as pd
            if isinstance(values, pd.Series):
                values = values.to_numpy()
        except ImportError:
            pass
        # Preserve dtype: numeric callers see float arrays; string callers
        # see object/str arrays.  np.asarray with no dtype= chooses based
        # on contents.
        self._values = np.asarray(values)
        if self._values.ndim != 1:
            raise ValueError("PerLC values must be 1-dimensional")

    @property
    def values(self):
        return self._values

    def __len__(self):
        return len(self._values)

    def __getitem__(self, i):
        item = self._values[i]
        # Unwrap numpy scalars to Python built-ins so callers see plain
        # int / float / str rather than numpy wrapper types.
        if hasattr(item, "item"):
            return item.item()
        return item

    def __repr__(self):
        return f"PerLC(n={len(self._values)})"
