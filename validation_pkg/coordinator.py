from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

from .exceptions import ConfigError, CoordinatorError


SUPPORTED_FORMATS = {"fasta", "genbank"}
EXT_TO_FORMAT = {
    ".fa": "fasta",
    ".fasta": "fasta",
    ".fna": "fasta",
    ".gb": "genbank",
    ".gbk": "genbank",
    ".genbank": "genbank",
}


@dataclass
class Config:
    input_path: Path
    output_dir: Path
    input_format: str
    validators: Dict[str, Any] = field(default_factory=dict)
    editors: Iterable[Dict[str, Any]] = field(default_factory=list)
    options: Dict[str, Any] = field(default_factory=dict)


class Coordinator:
    """Coordinates the validation pipeline.

    This first version focuses on *config management* and *exception integration*.
    The actual domain steps (parse → validate → edit → write) are intentionally
    stubbed and can be implemented later.
    """

    def __init__(self, config: Config):
        self.config = config

    # ---------- Public API ----------
    @classmethod
    def load(cls, config_input: Any) -> "Coordinator":
        """Load and validate configuration from a path, file-like, dict or JSON string.

        Parameters
        ----------
        config_input: Any
            One of:
            - str/Path to a JSON file
            - dict with configuration values
            - JSON string

        Returns
        -------
        Coordinator
            Coordinator instance with normalized Config.
        """
        raw = _read_config_any(config_input)
        cfg = _validate_and_normalize(raw)
        return cls(cfg)

    def summary(self) -> Dict[str, Any]:
        """Return a serializable summary of the normalized configuration."""
        return {
            "input_path": str(self.config.input_path),
            "output_dir": str(self.config.output_dir),
            "input_format": self.config.input_format,
            "validators": self.config.validators,
            "editors": list(self.config.editors),
            "options": self.config.options,
        }

    def run(self) -> Dict[str, Any]:
        """Stub run method.

        For now, we only demonstrate orchestration decisions based on config.
        Actual calls to parsers/validators/editors should be added later.
        """
        try:
            actions = {
                "parse_with": _choose_parser(self.config.input_format),
                "validators": _enabled_validators(self.config.validators),
                "planned_edits": list(self.config.editors),
                "outputs": str(self.config.output_dir),
            }
            # Here is where you'd call the real pipeline.
            # e.g., parser = FastaParser(...) / GenbankParser(...)
            #       records = parser.parse(self.config.input_path)
            #       ... then validators & editors ...
            return actions
        except KeyError as exc:  # defensive — should not be hit if validation ran
            raise CoordinatorError("Coordinator run failed due to missing key.", code="COORD/KEY", cause=exc) from exc
        except Exception as exc:
            raise CoordinatorError("Unexpected error during coordination.", code="COORD/UNEXPECTED", cause=exc) from exc


# ---------- Helper functions (module-private) ----------

def _read_config_any(config_input: Any) -> Dict[str, Any]:
    """Read raw config into a dict, raising ConfigError on problems."""
    try:
        if isinstance(config_input, dict):
            return dict(config_input)  # shallow copy
        if isinstance(config_input, (str, Path)):
            p = Path(config_input)
            if p.exists():
                try:
                    text = p.read_text(encoding="utf-8")
                except FileNotFoundError as exc:
                    raise ConfigError(f"Config file not found: {p}", code="CFG/NOT_FOUND", cause=exc) from exc
                except OSError as exc:
                    raise ConfigError(f"Cannot read config file: {p}", code="CFG/READ", cause=exc) from exc
            else:
                # If it's not a path, treat as JSON string
                text = str(config_input)
            try:
                return json.loads(text)
            except json.JSONDecodeError as exc:
                raise ConfigError("Invalid JSON for configuration.", code="CFG/JSON", cause=exc) from exc
        # Fallback: attempt to serialize unknown object via __dict__
        if hasattr(config_input, "__dict__"):
            return dict(getattr(config_input, "__dict__"))
    except ConfigError:
        raise
    except Exception as exc:
        raise ConfigError("Failed to read configuration.", code="CFG/READ_ANY", cause=exc) from exc

    raise ConfigError("Unsupported config input type.", code="CFG/TYPE")


def _validate_and_normalize(raw: Dict[str, Any]) -> Config:
    """Validate keys & types and build a Config dataclass."""
    if not isinstance(raw, dict):
        raise ConfigError("Configuration root must be an object/dict.", code="CFG/ROOT_TYPE")

    required = ["input_path", "output_dir"]
    missing = [k for k in required if k not in raw]
    if missing:
        raise ConfigError(f"Missing required config key(s): {', '.join(missing)}", code="CFG/MISSING")

    input_path = Path(raw["input_path"]).expanduser()
    output_dir = Path(raw["output_dir"]).expanduser()

    if not input_path.exists():
        raise ConfigError(f"Input file does not exist: {input_path}", code="CFG/INPUT_MISSING")

    # Determine input format
    input_format: Optional[str] = raw.get("input_format")
    if input_format is None:
        input_format = _guess_format_from_extension(input_path)
    input_format = str(input_format).lower().strip()
    if input_format not in SUPPORTED_FORMATS:
        raise ConfigError(
            f"Unsupported input_format '{input_format}'. Supported: {sorted(SUPPORTED_FORMATS)}",
            code="CFG/FORMAT",
        )

    validators = raw.get("validators", {})
    if validators is None:
        validators = {}
    if not isinstance(validators, dict):
        raise ConfigError("'validators' must be an object/dict.", code="CFG/VALIDATORS_TYPE")

    editors = raw.get("editors", [])
    if editors is None:
        editors = []
    if not isinstance(editors, (list, tuple)):
        raise ConfigError("'editors' must be an array/list.", code="CFG/EDITORS_TYPE")

    options = raw.get("options", {})
    if options is None:
        options = {}
    if not isinstance(options, dict):
        raise ConfigError("'options' must be an object/dict.", code="CFG/OPTIONS_TYPE")

    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except Exception as exc:
        raise ConfigError(f"Cannot create or access output_dir: {output_dir}", code="CFG/OUTPUT_DIR", cause=exc) from exc

    return Config(
        input_path=input_path,
        output_dir=output_dir,
        input_format=input_format,
        validators=validators,
        editors=list(editors),
        options=options,
    )


def _guess_format_from_extension(path: Path) -> str:
    ext = path.suffix.lower()
    if ext in EXT_TO_FORMAT:
        return EXT_TO_FORMAT[ext]
    # Second chance: double extension like .fa.gz → use the first suffix
    if len(path.suffixes) >= 2:
        first_ext = path.suffixes[-2].lower()
        if first_ext in EXT_TO_FORMAT:
            return EXT_TO_FORMAT[first_ext]
    # As a safe default, prefer fasta if the filename hints so
    name = path.name.lower()
    if any(tok in name for tok in ("fasta", "fna", "fa")):
        return "fasta"
    if any(tok in name for tok in ("gb", "gbk", "genbank")):
        return "genbank"
    raise ConfigError(
        f"Cannot infer input format from extension of '{path.name}'. Provide 'input_format'.",
        code="CFG/FORMAT_GUESS",
    )


def _choose_parser(fmt: str) -> str:
    """Return the dotted path of the parser class to use (deferred import).

    Keeping this as a string avoids importing heavy deps at import time and keeps
    tests fast. The actual import/usage can be done later.
    """
    if fmt == "fasta":
        return "parsers.fasta_parser.FastaParser"
    if fmt == "genbank":
        return "parsers.genbank_parser.GenbankParser"
    # Should be unreachable due to validation
    raise CoordinatorError(f"Unknown format '{fmt}'", code="COORD/PARSER")


def _enabled_validators(validators_cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Filter/normalize validators configuration to only enabled ones.

    Example input:
        {"read": true, "genome": {"threshold": 0.95}, "feature": false}
    Output keeps falsy entries out and normalizes bool → {} for convenience.
    """
    enabled: Dict[str, Any] = {}
    for name, cfg in validators_cfg.items():
        if isinstance(cfg, bool):
            if cfg:
                enabled[name] = {}
        elif isinstance(cfg, dict):
            enabled[name] = cfg
        else:
            raise ConfigError(f"Validator '{name}' must be bool or object.", code="CFG/VALIDATOR_ITEM")
    return enabled
