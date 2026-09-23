"""Reconcile full-table coverage against independently admitted replay batches.

This is a final-admission component, not an execution or database admission.
Callers must bind all inputs to freshly verified provenance records.
"""


def unique_ids(values, label):
    values = list(values)
    if any(not isinstance(value, str) or not value for value in values):
        raise ValueError("Invalid identifier in " + label)
    result = set(values)
    if len(result) != len(values):
        raise ValueError("Duplicate identifier in " + label)
    return result


def reconcile(genes, retained_queries, retained_rows, batches, content):
    universe = unique_ids(genes, "input")
    retained = unique_ids(retained_queries, "retained prefix")
    if not universe or type(retained_rows) is not int or retained_rows < len(retained):
        raise ValueError("Invalid retained row count or empty input")
    replay, absent, failed = set(), set(), set()
    replay_rows = 0
    for batch in batches:
        queries = unique_ids(batch["query_ids"], "batch queries")
        hits = unique_ids((b["query"] for b in batch["query_blocks"]), "batch hit blocks")
        no_hits = unique_ids(batch["coverage"]["queries_without_hits"], "batch no-hit queries")
        failures = unique_ids(batch["coverage"]["failed_queries"], "batch failed queries")
        diagnostics = batch["diagnostics"]
        if (not set(diagnostics) <= queries
                or any(type(d["query_failed"]) is not bool for d in diagnostics.values())
                or {g for g, d in diagnostics.items() if d["query_failed"]} != failures):
            raise ValueError("Batch failure diagnostics disagree")
        if replay & queries or not hits <= queries or no_hits != queries - hits:
            raise ValueError("Inconsistent or overlapping replay coverage")
        if not failures <= no_hits:
            raise ValueError("Failed query outside no-hit replay disposition")
        for block in batch["query_blocks"]:
            if type(block["rows"]) is not int or block["rows"] <= 0:
                raise ValueError("Invalid replay block row count")
            replay_rows += block["rows"]
        replay.update(queries)
        absent.update(no_hits)
        failed.update(failures)
    if retained & replay or retained | replay != universe:
        raise ValueError("Retained and replay queries do not partition input")
    actual_absent = unique_ids(content["query_ids_without_hits"], "full-table no-hit queries")
    diagnostics = content["diagnostics"]
    diagnostic_ids = unique_ids((d["gene"] for d in diagnostics), "full-table diagnostics")
    if not diagnostic_ids <= universe:
        raise ValueError("Unknown full-table diagnostic query")
    if any(type(d["query_failed"]) is not bool for d in diagnostics):
        raise ValueError("Invalid diagnostic failure flag")
    actual_failed = {d["gene"] for d in diagnostics if d["query_failed"]}
    if actual_absent != absent or actual_failed != failed:
        raise ValueError("Full-table query dispositions disagree with replay evidence")
    for diagnostic in diagnostics:
        if diagnostic["has_query_hits"] is not (diagnostic["gene"] not in absent):
            raise ValueError("Diagnostic outgoing-hit flag disagrees")
    expected = dict(input_proteins=len(universe), hsp_rows=retained_rows + replay_rows,
        queries_with_hits=len(universe - absent),
        queries_without_hits=len(absent), failed_queries=len(failed),
        queries_without_hits_and_without_logged_failure=len(absent - failed),
        failed_queries_with_query_hits=0, hsp_rows_above_1e_minus_5=0)
    for name, value in expected.items():
        if type(content[name]) is not int or content[name] != value:
            raise ValueError("Full-table count disagrees: " + name)
    return dict(status="recovered_query_dispositions_reconciled", **expected,
        retained_queries=len(retained), replay_queries=len(replay), replay_hsp_rows=replay_rows,
        search_admitted=False, publication_ready=False)
