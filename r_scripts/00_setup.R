# 00_setup.R
# One-off environment checks. Run this once per machine before anything else.
# Nothing here writes to data/ or outputs/. It exists only to confirm that the
# four things the rest of the pipeline depends on are actually working.


# Section 1: Initialize environment
# renv::init() creates the project library and the lockfile. The console must be
# restarted after it finishes, otherwise the new library is not on the search path.
renv::init()

# renv::status() compares the lockfile against the installed library. A clean
# report here means the project can be rebuilt from renv.lock alone.
renv::status()


# Section 2: Establish project folder path
# here() should return the project root. Every path in scripts 01 to 04 is built
# from here(), so if this points somewhere unexpected the rest of the pipeline fails.
here::here()


# Section 3: Confirm Java is available to R
# MaxEnt runs on Java through rJava. A Java version mismatch does not surface until
# maxent() is first called, deep inside the cross-validation loops, so it is checked
# up front instead.
rJava::.jvmState()


# Section 4: Confirm predicts can find maxent.jar
# Called with no arguments, MaxEnt() reports whether the jar is installed. If it is
# not, the regularised GLM in 03_analysis.R still runs but MaxEnt does not.
predicts::MaxEnt()
