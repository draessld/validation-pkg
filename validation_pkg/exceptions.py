from __future__ import annotations
from typing import Optional




class ValidationError(Exception):
    """Base exception for the EFSA project.


    Includes an optional error code and a causal exception for improved debugging.
    """


    def __init__(self, message: str = "", *, code: Optional[str] = None, cause: Optional[BaseException] = None):
        super().__init__(message)
        self.code = code
        self.cause = cause


    def __str__(self) -> str: # pragma: no cover - trivial
        base = super().__str__()
        parts = [base]
        if self.code:
            parts.append(f"[code={self.code}]")
        if self.cause:
            parts.append(f"(caused by {self.cause.__class__.__name__}: {self.cause})")
        return " ".join(parts)




class ConfigError(ValidationError):
    """Raised for problems loading or validating configuration."""




class CoordinatorError(ValidationError):
    """Raised by the coordinator on orchestration/runtime errors."""




class ParseError(ValidationError):
    """Raised when parsing input files fails."""




class ValidationError(ValidationError):
    """Raised by validators when biological or structural validation fails."""




class EditError(ValidationError):
    """Raised by editors when an edit operation cannot be applied."""




__all__ = [
    "ValidationError",
    "ConfigError",
    "CoordinatorError",
    "ParseError",
    "ValidationError",
    "EditError",
]