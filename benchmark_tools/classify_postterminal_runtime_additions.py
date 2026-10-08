"""Classify exact postterminal inventory additions; never authorize an attempt."""

from datetime import datetime


def require(condition, message):
    if not condition:
        raise ValueError(message)


def classify(comparison, allowed_additions, native_end, prior_inventory_latest_at, installed_at):
    require(all(isinstance(value, datetime) and value.utcoffset() is not None for value in
                (native_end, prior_inventory_latest_at, installed_at)),
            "Require explicit timezone-aware chronology")
    require(native_end <= prior_inventory_latest_at < installed_at,
            "Additions are not established strictly after native terminal inventory")
    require(comparison.get("equal") is False and comparison.get("missing") == []
        and comparison.get("changed") == [] and comparison.get("metadata_differences") == {},
        "Require addition-only difference; no changed, missing or metadata entries")
    require(isinstance(allowed_additions, list) and bool(allowed_additions)
        and comparison.get("added") == allowed_additions,
        "Added entries differ from the exact prospective inventory boundary")
    paths = []
    for row in allowed_additions:
        require(set(row) == {"path", "mode", "kind", "bytes", "sha256"}
            and row["kind"] == "file" and isinstance(row["path"], str) and row["path"].startswith("/")
            and type(row["mode"]) is int and 0 <= row["mode"] <= 0o7777
            and type(row["bytes"]) is int and row["bytes"] > 0
            and isinstance(row["sha256"], str) and len(row["sha256"]) == 64
            and all(char in "0123456789abcdef" for char in row["sha256"]),
            "Require exact regular-file additions with mode, bytes and SHA256")
        paths.append(row["path"])
    require(paths == sorted(set(paths)), "Require sorted unique additional paths")
    require(type(comparison.get("expected_records")) is int and comparison["expected_records"] > 0
        and type(comparison.get("observed_records")) is int
        and comparison["observed_records"] == comparison["expected_records"] + len(allowed_additions),
        "Inventory record totals do not recover exact additions")
    return dict(status="postterminal_additions_classified_not_admitted",
        original_inventory_equality=False, original_entries_changed_or_missing=False,
        additional_entries=allowed_additions, native_end=native_end.isoformat(),
        prior_inventory_latest_at=prior_inventory_latest_at.isoformat(), installed_at=installed_at.isoformat(),
        historical_failure_cause_established=False, terminal_reviewed=False,
        next_identity_authorized=False, accuracy_evaluated=False, automatic_retry=False,
        publication_ready=False, limitations=[
            "Classification uses supplied retained comparison and timestamps; callers must bind their provenance.",
            "This is not a replacement for the original exact runtime verifier or full-review admission.",
            "Before/after inventory and installation records do not establish continuous runtime integrity.",
            "Any prospective history/conversion or runtime adoption needs a separate explicit contract."])
