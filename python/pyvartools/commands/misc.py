"""Miscellaneous vartools command wrappers."""

from typing import List, Optional, Union

from pyvartools._command import VartoolsCommand
from pyvartools.userlib import UserCommand  # noqa: F401 — re-exported
from ._helpers import _flag, _bool, _outtoken, _pval, _varexpr


class addfitskeyword(VartoolsCommand):
    """Add a FITS keyword to the output statistics table (-addfitskeyword).

    Parameters
    ----------
    keyword : str
        FITS keyword name (max 8 characters).
    dtype : str
        Data type: ``"TDOUBLE"``, ``"TINT"``, ``"TLONG"``, or ``"TSTRING"``.
    value : str
        Value specification.  Either a bare Python scalar (``"fix"`` is
        prepended automatically) or a full vartools token string such as
        ``"fix 3.14"``, ``"var myvar"``.
    comment : str, optional
        FITS comment string.
    hdu : str, optional
        ``"primary"`` or ``"extension"`` (vartools default is primary).
    mode : str, optional
        ``"append"`` or ``"update"``.
    combinelc : str, optional
        Variable name holding the LC number for combined-LC mode.
    """

    _vt_name = "addfitskeyword"

    def __init__(
        self,
        keyword: str,
        dtype: str,
        value,
        comment: Optional[str] = None,
        hdu: Optional[str] = None,
        mode: Optional[str] = None,
        combinelc: Optional[str] = None,
    ) -> None:
        self.keyword = keyword
        self.dtype = dtype
        self.value = value
        self.comment = comment
        self.hdu = hdu
        self.mode = mode
        self.combinelc = combinelc

    def _to_cli_args(self) -> List[str]:
        args = ["-addfitskeyword", str(self.keyword)]
        if self.combinelc is not None:
            args += ["combinelc", str(self.combinelc)]
        args += [str(self.dtype)]
        # Auto-prefix bare scalars with "fix"
        val = self.value
        if isinstance(val, (int, float)):
            args += ["fix", str(val)]
        else:
            args += str(val).split()
        if self.comment is not None:
            args += ["comment", str(self.comment)]
        if self.hdu is not None:
            args += [str(self.hdu)]
        if self.mode is not None:
            args += [str(self.mode)]
        return args

    def _output_file_specs(self) -> dict:
        return {}


class converttime(VartoolsCommand):
    """Convert the time column between time systems (-converttime).

    Parameters
    ----------
    input_format : str
        Input time format: ``"mjd"``, ``"jd"``, ``"hjd"``, or ``"bjd"``.
    output_format : str
        Output time format: ``"mjd"``, ``"jd"``, ``"hjd"``, or ``"bjd"``.
    ra : float or str, optional
        Right ascension for HJD/BJD conversion.  A float is treated as
        degrees and passed as ``"fix ra dec"`` (requires ``dec`` too).
        A string is split and forwarded verbatim (e.g. ``"list"``).
    dec : float, optional
        Declination in degrees (used together with a float ``ra``).
    input_subtract : float, optional
        Subtract this constant from the input times before conversion.
    output_subtract : float, optional
        Subtract this constant from the output times after conversion.
    input_sys : str, optional
        Input time system: ``"tdb"`` or ``"utc"``.
    output_sys : str, optional
        Output time system: ``"tdb"`` or ``"utc"``.
    ephemfile : str, optional
        Path to a JPL ephemeris file.
    leapsecfile : str, optional
        Path to a leap-second table.
    """

    _vt_name = "converttime"

    def __init__(
        self,
        input_format: str,
        output_format: str,
        ra=None,
        dec: Optional[float] = None,
        input_subtract: Optional[float] = None,
        output_subtract: Optional[float] = None,
        input_sys: Optional[str] = None,
        output_sys: Optional[str] = None,
        ephemfile: Optional[str] = None,
        leapsecfile: Optional[str] = None,
    ) -> None:
        self.input_format = input_format
        self.output_format = output_format
        self.ra = ra
        self.dec = dec
        self.input_subtract = input_subtract
        self.output_subtract = output_subtract
        self.input_sys = input_sys
        self.output_sys = output_sys
        self.ephemfile = ephemfile
        self.leapsecfile = leapsecfile

    def _to_cli_args(self) -> List[str]:
        args = ["-converttime",
                "input", str(self.input_format)]
        if self.input_subtract is not None:
            args += ["inputsubtract", str(self.input_subtract)]
        if self.input_sys is not None:
            args += [f"inputsys-{self.input_sys}"]
        args += ["output", str(self.output_format)]
        if self.output_subtract is not None:
            args += ["outputsubtract", str(self.output_subtract)]
        if self.output_sys is not None:
            args += [f"outputsys-{self.output_sys}"]
        if self.ra is not None:
            if isinstance(self.ra, (int, float)):
                args += ["radec", "fix", str(self.ra), str(self.dec)]
            else:
                args += ["radec"] + str(self.ra).split()
        if self.ephemfile is not None:
            args += ["ephemfile", str(self.ephemfile)]
        if self.leapsecfile is not None:
            args += ["leapsecfile", str(self.leapsecfile)]
        return args

    def _output_file_specs(self) -> dict:
        return {}


class R(VartoolsCommand):
    """Run R code on each light curve (-R).

    Parameters
    ----------
    command : str
        Either an inline R code string or, if ``fromfile=True``, the path
        to an R script file.
    fromfile : bool
        If True, ``command`` is treated as a file path (passes
        ``"fromfile"`` keyword).  Default False (inline string).
    init : str, optional
        Initialisation R code (inline string) or file path run once before
        processing starts.
    init_fromfile : bool
        If True, ``init`` is a file path (passes ``"file"`` keyword).
    vars : str, optional
        Comma-separated list of vartools variables to pass in and out.
    invars : str, optional
        Variables to pass into R (alternative to ``vars``).
    outvars : str, optional
        Variables to receive back from R (alternative to ``vars``).
    outputcolumns : str, optional
        Variables to append as extra output columns.
    process_all_lcs : bool
        Pass ``"process_all_lcs"`` flag.
    verbose : bool
        Pass ``"verbose"`` flag.
    continueprocess : int, optional
        Prior R command number to share a subprocess with.
    """

    _vt_name = "R"

    def __init__(
        self,
        command: str,
        fromfile: bool = False,
        init: Optional[str] = None,
        init_fromfile: bool = False,
        vars: Optional[str] = None,
        invars: Optional[str] = None,
        outvars: Optional[str] = None,
        outputcolumns: Optional[str] = None,
        process_all_lcs: bool = False,
        verbose: bool = False,
        continueprocess: Optional[int] = None,
    ) -> None:
        self.command = command
        self.fromfile = fromfile
        self.init = init
        self.init_fromfile = init_fromfile
        self.vars = vars
        self.invars = invars
        self.outvars = outvars
        self.outputcolumns = outputcolumns
        self.process_all_lcs = process_all_lcs
        self.verbose = verbose
        self.continueprocess = continueprocess

    def _to_cli_args(self) -> List[str]:
        args = ["-R"]
        if self.fromfile:
            args += ["fromfile", str(self.command)]
        else:
            args += [str(self.command)]
        if self.continueprocess is not None:
            args += ["continueprocess", str(self.continueprocess)]
        if self.init is not None:
            if self.init_fromfile:
                args += ["init", "file", str(self.init)]
            else:
                args += ["init", str(self.init)]
        if self.vars is not None:
            args += ["vars", str(self.vars)]
        else:
            if self.invars is not None:
                args += ["invars", str(self.invars)]
            if self.outvars is not None:
                args += ["outvars", str(self.outvars)]
        if self.outputcolumns is not None:
            args += ["outputcolumns", str(self.outputcolumns)]
        if self.process_all_lcs:
            args += ["process_all_lcs"]
        if self.verbose:
            args += ["verbose"]
        return args

    def _output_file_specs(self) -> dict:
        return {}


class match(VartoolsCommand):
    """Match the light curve against an external catalog file (-match).

    Parameters
    ----------
    catalog : str
        Path to the catalog file.  When ``source="inlist"``, use
        ``inlist_column`` instead to specify the column holding the catalog
        path, and ``catalog`` is ignored.
    source : str
        ``"file"`` (default) or ``"inlist"``.
    matchcolumn : str
        Column specification used for matching, e.g. ``"t:1"``
        (vartools-variable:column-number) or just ``"1"``.
    addcolumns : str
        Comma-separated column specs to add to the light curve, e.g.
        ``"ra:2,dec:3"`` or ``"ra:2:TDOUBLE,dec:3:TDOUBLE"``.
    missing : str
        How to handle unmatched rows: ``"cullmissing"``, ``"nanmissing"``,
        or ``"missingval <value>"``.
    inlist_column : str or int, optional
        Column number/name in the input list that holds per-LC catalog
        file paths.  Required when ``source="inlist"``; ``catalog`` is
        then ignored.
    skipnum : int, optional
        Number of header lines to skip in the catalog.
    skipchar : str, optional
        Comment character(s) to skip (comma-separated if multiple).
    delimiter : str, optional
        Column delimiter character.
    opencommand : str, optional
        Shell command used to open the catalog file (e.g. for compressed files).
    """

    _vt_name = "match"

    def __init__(
        self,
        catalog: str,
        matchcolumn: str,
        addcolumns: str,
        missing: str = "nanmissing",
        source: str = "file",
        inlist_column: Optional[Union[str, int]] = None,
        skipnum: Optional[int] = None,
        skipchar: Optional[str] = None,
        delimiter: Optional[str] = None,
        opencommand: Optional[str] = None,
    ) -> None:
        self.catalog = catalog
        self.matchcolumn = matchcolumn
        self.addcolumns = addcolumns
        self.missing = missing
        self.source = source
        self.inlist_column = inlist_column
        self.skipnum = skipnum
        self.skipchar = skipchar
        self.delimiter = delimiter
        self.opencommand = opencommand

    def _to_cli_args(self) -> List[str]:
        if self.source == "inlist":
            if self.inlist_column is None:
                raise ValueError(
                    "cmd.match: inlist_column is required when source='inlist'"
                )
            args = ["-match", "inlist", str(self.inlist_column)]
        else:
            args = ["-match", str(self.source), str(self.catalog)]
        if self.opencommand is not None:
            args += ["opencommand", str(self.opencommand)]
        if self.skipnum is not None:
            args += ["skipnum", str(self.skipnum)]
        if self.skipchar is not None:
            args += ["skipchar", str(self.skipchar)]
        if self.delimiter is not None:
            args += ["delimiter", str(self.delimiter)]
        args += ["matchcolumn", str(self.matchcolumn)]
        args += ["addcolumns", str(self.addcolumns)]
        args += str(self.missing).split()
        return args

    def _output_file_specs(self) -> dict:
        return {}


class o(VartoolsCommand):
    """Output the current light curve to a file (-o).

    Parameters
    ----------
    filename : str, optional
        Output file path (or directory in batch mode).  Required unless
        ``capture=True``.
    nameformat : str, optional
        Format string for constructing per-LC output filenames in batch
        mode, e.g. ``"file_%s_%05d.txt"`` (``%s`` = LC basename, ``%d``
        = sequence number).
    columnformat : str, optional
        Output column specification, e.g.
        ``"t:%17.9f,mag:%9.5f,err:%9.5f"``.
    fits : bool
        Write output in FITS format.
    noclobber : bool
        Do not overwrite existing output files.
    copyheader : bool
        Copy the input FITS header to the output file.
    namecommand : str, optional
        Shell command used to generate the output filename.
    namefromlist : bool or str, optional
        Derive output filename from the input list.  Pass ``True`` to use
        the default column, or a column number/name string to use a
        specific column (emits ``"namefromlist" "column" col``).
    delimiter : str, optional
        Column delimiter character for the output file.
    logcommandline : bool
        Write the vartools command line to the output file header.
    capture : bool
        If ``True``, capture the written light curve and return it in
        ``result.files[key]``.  For single-LC runs, ``result.files[key]``
        is a ``LightCurve``; for batch runs it is a list of
        ``LightCurve`` objects (one per input LC, ``None`` if missing).

        When ``filename`` is ``None`` the output is written to a
        temporary file/directory that is cleaned up automatically after
        the run.  When ``filename`` is also supplied the file is written
        to disk *and* captured.  Default ``False``.
    key : str
        Key under which the captured LC(s) are stored in
        ``result.files``.  Default ``"o"``.  Use a unique key when the
        pipeline contains more than one ``cmd.o(capture=True)`` command.
    """

    _vt_name = "o"

    def __init__(
        self,
        filename: Optional[str] = None,
        nameformat: Optional[str] = None,
        columnformat: Optional[str] = None,
        allcols: bool = False,
        fits: bool = False,
        noclobber: bool = False,
        copyheader: bool = False,
        namecommand: Optional[str] = None,
        namefromlist: Union[bool, str, None] = None,
        delimiter: Optional[str] = None,
        logcommandline: bool = False,
        capture: bool = False,
        key: str = "o",
    ) -> None:
        if filename is None and not capture:
            raise ValueError(
                "cmd.o() requires either a filename or capture=True"
            )
        if allcols and columnformat is not None:
            raise ValueError(
                "cmd.o(): 'allcols' and 'columnformat' are mutually exclusive"
            )
        self.filename = filename
        self.nameformat = nameformat
        self.columnformat = columnformat
        self.allcols = allcols
        self.fits = fits
        self.noclobber = noclobber
        self.copyheader = copyheader
        self.namecommand = namecommand
        self.namefromlist = namefromlist
        self.delimiter = delimiter
        self.logcommandline = logcommandline
        self.capture = capture
        self.key = key
        # Injected by Pipeline before _to_cli_args() is called when
        # capture=True and filename is None.
        self._capture_path: Optional[str] = None

    def _to_cli_args(self) -> List[str]:
        path = self.filename if self.filename is not None else self._capture_path
        if path is None:
            raise RuntimeError(
                "cmd.o with capture=True must be run through a Pipeline "
                "(capture path has not been assigned yet)"
            )
        args = ["-o", str(path)]
        if self.nameformat is not None:
            args += ["nameformat", str(self.nameformat)]
        if self.columnformat is not None:
            args += ["columnformat", str(self.columnformat)]
        elif self.allcols or self.capture:
            # When capturing, default to `allcols` so the captured DataFrame
            # contains every LC-vector variable defined by earlier commands —
            # matching the library-mode fast path.  The explicit allcols flag
            # also takes this branch for non-capturing callers.
            args += ["allcols"]
        if self.fits:
            args += ["fits"]
        if self.noclobber:
            args += ["noclobber"]
        if self.copyheader:
            args += ["copyheader"]
        if self.namecommand is not None:
            args += ["namecommand", str(self.namecommand)]
        if self.namefromlist is not None and self.namefromlist is not False:
            if self.namefromlist is True:
                args += ["namefromlist"]
            else:
                args += ["namefromlist", "column", str(self.namefromlist)]
        if self.delimiter is not None:
            args += ["delimiter", str(self.delimiter)]
        if self.logcommandline:
            args += ["logcommandline"]
        return args

    def _output_file_specs(self) -> dict:
        return {}

    # ------------------------------------------------------------------
    # Capture helpers
    # ------------------------------------------------------------------

    def _columnformat_names(self):
        """Return the list of column names declared in ``self.columnformat``.

        Returns ``None`` when no ``columnformat`` is set.  Each entry in
        the comma-separated spec is ``name[:printf_format]`` — strip the
        format suffix and keep the bare name, so the captured ASCII
        DataFrame can be renamed from auto-generated ``col4``/``col5``/…
        to the user-declared names.
        """
        if not self.columnformat:
            return None
        names = []
        for entry in str(self.columnformat).split(","):
            entry = entry.strip()
            if not entry:
                continue
            names.append(entry.split(":", 1)[0])
        return names or None


class ifcmd(VartoolsCommand):
    """Open a conditional block (-if).

    Parameters
    ----------
    condition : str
        The vartools condition expression (passed verbatim).

    See Also
    --------
    elifcmd : ``-elif`` branch.
    elsecmd : ``-else`` branch.
    ficmd : ``-fi`` — closes an `ifcmd` / `elifcmd` / `elsecmd` block.
    """

    _vt_name = "if"

    def __init__(self, condition: str) -> None:
        self.condition = condition

    def _to_cli_args(self) -> List[str]:
        # vartools -if takes the condition as a single token; do not split
        # on whitespace.
        return ["-if", str(self.condition)]

    def _output_file_specs(self) -> dict:
        return {}


class elifcmd(VartoolsCommand):
    """Open a conditional `-elif` branch.

    Must be preceded by a matching :class:`ifcmd` and closed by a
    :class:`ficmd`.

    Parameters
    ----------
    condition : str
        The vartools condition expression (passed verbatim).
    """

    _vt_name = "elif"

    def __init__(self, condition: str) -> None:
        self.condition = condition

    def _to_cli_args(self) -> List[str]:
        return ["-elif", str(self.condition)]

    def _output_file_specs(self) -> dict:
        return {}


class elsecmd(VartoolsCommand):
    """Open a conditional `-else` branch.

    Must be preceded by a matching :class:`ifcmd` (or :class:`elifcmd`) and
    closed by a :class:`ficmd`.  Takes no parameters.
    """

    _vt_name = "else"

    def __init__(self) -> None:
        pass

    def _to_cli_args(self) -> List[str]:
        return ["-else"]

    def _output_file_specs(self) -> dict:
        return {}


class ficmd(VartoolsCommand):
    """Close a conditional block (`-fi`).

    Must follow a matching :class:`ifcmd` / :class:`elifcmd` / :class:`elsecmd`
    sequence.  Takes no parameters.
    """

    _vt_name = "fi"

    def __init__(self) -> None:
        pass

    def _to_cli_args(self) -> List[str]:
        return ["-fi"]

    def _output_file_specs(self) -> dict:
        return {}


class binlc(VartoolsCommand):
    """Bin the light curve in time (-binlc).

    Parameters
    ----------
    method : str
        Binning method: ``"average"`` (default), ``"median"``, or
        ``"weightedaverage"``.
    binsize : float, optional
        Bin size in the same units as the time column.  Either ``binsize``
        or ``nbins`` must be given.
    nbins : int, optional
        Number of bins (alternative to ``binsize``).
    time_output : str
        How to set the output time for each bin: ``"tcenter"`` (default),
        ``"taverage"``, ``"tmedian"``, or ``"tnoshrink"``.
    bincolumns : str, optional
        Extra columns to bin, e.g. ``"col1,col2:median"``.
    bincolumnsonly : bool
        Only output the ``bincolumns`` variables (with ``"tnoshrink"``).
    T0 : float or str, optional
        Reference time for bin-edge alignment.  A float is passed as
        ``"fix T0"``.
    firstbinshift : float, optional
        Shift the first bin edge by this amount.
    maskpoints : str, optional
    """

    _vt_name = "binlc"

    def __init__(
        self,
        method: str = "average",
        binsize: Optional[float] = None,
        nbins: Optional[int] = None,
        time_output: str = "tcenter",
        bincolumns: Optional[str] = None,
        bincolumnsonly: bool = False,
        T0=None,
        firstbinshift: Optional[float] = None,
        maskpoints: Optional[str] = None,
    ) -> None:
        if binsize is None and nbins is None:
            raise ValueError("binlc requires either binsize or nbins")
        self.method = method
        self.binsize = binsize
        self.nbins = nbins
        self.time_output = time_output
        self.bincolumns = bincolumns
        self.bincolumnsonly = bincolumnsonly
        self.T0 = T0
        self.firstbinshift = firstbinshift
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        args = ["-binlc", str(self.method)]
        if self.binsize is not None:
            args += ["binsize"] + _varexpr(self.binsize)
        else:
            args += ["nbins"] + _varexpr(self.nbins)
        if self.bincolumns is not None:
            args += ["bincolumns", str(self.bincolumns)]
        if self.T0 is not None:
            if isinstance(self.T0, (int, float)):
                args += ["T0", "fix", str(self.T0)]
            else:
                args += ["T0"] + str(self.T0).split()
        if self.firstbinshift is not None:
            args += ["firstbinshift"] + _varexpr(self.firstbinshift)
        to = str(self.time_output)
        args += [to]
        if to == "tnoshrink" and self.bincolumnsonly:
            args += ["bincolumnsonly"]
        if self.maskpoints is not None:
            args += ["maskpoints", str(self.maskpoints)]
        return args

    def _output_file_specs(self) -> dict:
        return {}


class columnsuffix(VartoolsCommand):
    """Set the column-name suffix for all subsequent commands (-columnsuffix).

    vartools appends a numeric suffix (the command's 0-based index) to output
    variable names by default, e.g. ``LS_Period_1_0``.  Use this command to
    replace that suffix with a meaningful string so the key is predictable
    regardless of pipeline position.

    Parameters
    ----------
    suffix : str
        Suffix to use for all commands that follow in the pipeline.

    Examples
    --------
    ::

        pipe = vt.Pipeline([
            cmd.clip(sigclip=5.0),
            cmd.columnsuffix("ls"),
            cmd.LS(0.5, 10.0, 1e-3),
        ])
        result = pipe.run(lc)
        best_period = float(result.vars["LS_Period_1_ls"])
    """

    _vt_name = "columnsuffix"

    def __init__(self, suffix: str) -> None:
        self.suffix = suffix

    def _to_cli_args(self) -> List[str]:
        return ["-columnsuffix", str(self.suffix)]

    def _output_file_specs(self) -> dict:
        return {}


class Raw(VartoolsCommand):
    """Pass arbitrary raw CLI tokens directly to vartools.

    Use this as an escape hatch for commands not yet wrapped or for
    experimental options.

    Parameters
    ----------
    args : list of str or single str
        Raw tokens to inject verbatim.  A plain string is split on
        whitespace; a list is used as-is.

    Examples
    --------
    ::

        Raw("-LS 0.1 10.0 0.1 5 1 0")
        Raw(["-LS", "0.1", "10.0", "0.1", "5", "1", "0"])
    """

    _vt_name = ""  # no canonical name; skip help population

    def __init__(self, args: Union[str, List[str]]) -> None:
        if isinstance(args, str):
            self._raw_args = args.split()
        else:
            self._raw_args = list(args)

    def _to_cli_args(self) -> List[str]:
        return self._raw_args

    def _output_file_specs(self) -> dict:
        return {}

    @classmethod
    def _populate_help(cls) -> None:
        # Raw has no -help entry; set a short doc manually.
        cls.__doc__ = (
            "Pass arbitrary raw CLI tokens directly to vartools.\n\n"
            "Use as an escape hatch for commands not yet wrapped."
        )
