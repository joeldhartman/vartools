"""
Command wrapper classes for pyvartools.

Each class corresponds to one vartools command.  All are subclasses of
VartoolsCommand defined in pyvartools._command.
"""

from pyvartools._command import VartoolsCommand

# Periodicity / period-search commands
from .periodicity import (
    LS,
    aov,
    aov_harm,
    PDM,
    FTP,
    BLS,
    BLSFixPer,
    BLSFixDurTc,
    BLSFixPerDurTc,
    autocorrelation,
    dftclean,
    wwz,
    GetLSAmpThresh,
    Phase,
)

# Light-curve manipulation / statistics commands
from .manipulation import (
    clip,
    rms,
    rmsbin,
    chi2,
    chi2bin,
    alarm,
    vonNeumann,
    percentileratios,
    beyondNsigma,
    slopestats,
    CodyM,
    CodyQ,
    structurefunction,
    drwfit,
    rescalesig,
    ensemblerescalesig,
    stats,
    harmonicfilter,
    fourierfilter,
    Killharm,
    linfit,
    Injectharm,
    Injecttransit,
    sortlc,
    restricttimes,
    restoretimes,
    savelc,
    restorelc,
    difffluxtomag,
    fluxtomag,
    magtoflux,
    changeerror,
    changevariable,
    copylc,
    medianfilter,
    expr,
    print_cols,
    FFT,
    IFFT,
    resample,
    decorr,
    Jstet,
)

# Model-fitting commands
from .fitting import (
    TFA,
    TFA_SR,
    SYSREM,
    MandelAgolTransit,
    SoftenedTransit,
    Starspot,
    microlens,
    nonlinfit,
    addnoise,
    findblends,
    MatchedFilter,
)

# Miscellaneous / utility commands
from .misc import (
    addfitskeyword,
    converttime,
    R,
    python,
    match,
    o,
    ifcmd,
    elifcmd,
    elsecmd,
    ficmd,
    binlc,
    columnsuffix,
    Raw,
    UserCommand,
)

# Typed wrappers for USERLIB extension commands
from .userlibs import (
    fastchi2,
    ftuneven,
    hatpiflag,
    jktebop,
    macula,
    magadd,
    splinedetrend,
    stitch,
)

__all__ = [
    "VartoolsCommand",
    # periodicity
    "LS", "aov", "aov_harm", "PDM", "FTP", "BLS", "BLSFixPer", "BLSFixDurTc",
    "BLSFixPerDurTc", "autocorrelation", "dftclean", "wwz",
    "GetLSAmpThresh", "Phase",
    # manipulation
    "clip", "rms", "rmsbin", "chi2", "chi2bin", "alarm", "vonNeumann",
    "percentileratios", "beyondNsigma", "slopestats", "CodyM", "CodyQ",
    "structurefunction", "drwfit",
    "rescalesig", "ensemblerescalesig", "stats", "harmonicfilter",
    "fourierfilter", "Killharm", "linfit",
    "Injectharm", "Injecttransit", "sortlc", "restricttimes", "restoretimes",
    "savelc", "restorelc", "difffluxtomag", "fluxtomag", "magtoflux", "changeerror",
    "changevariable", "copylc", "medianfilter", "expr", "print_cols",
    "FFT", "IFFT", "resample", "decorr", "Jstet",
    # fitting
    "TFA", "TFA_SR", "SYSREM", "MandelAgolTransit", "SoftenedTransit",
    "Starspot", "microlens", "nonlinfit", "addnoise", "findblends",
    "MatchedFilter",
    # misc
    "addfitskeyword", "converttime", "R", "python", "match", "o",
    "ifcmd", "elifcmd", "elsecmd", "ficmd",
    "binlc", "columnsuffix", "Raw",
    # user extensions
    "UserCommand",
    # typed USERLIB wrappers
    "fastchi2", "ftuneven", "hatpiflag", "jktebop", "macula",
    "magadd", "splinedetrend", "stitch",
]
