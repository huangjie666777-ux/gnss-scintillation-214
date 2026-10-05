"""Shared exception types with file-position context."""


class GnssError(Exception):
    """Base error carrying a human-readable file location."""

    def __init__(self, message: str, file: str = "", line: int = 0):
        self.file = file
        self.line = line
        loc = f"{file}:{line}" if file else "input"
        super().__init__(f"{loc}: {message}")


class RinexParseError(GnssError):
    pass


class Sp3ParseError(GnssError):
    pass


class RejectedContentError(GnssError):
    """Raised for explicitly unsupported content (events, corrections, etc.)."""
