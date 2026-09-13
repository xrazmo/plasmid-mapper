class PipelineError(Exception):
    """User-facing error with a plain-English remediation message.

    Raised instead of letting subprocess/IO errors propagate as raw
    tracebacks, since the target audience is microbiologists running this
    from the command line, not developers debugging a stack trace.
    """
