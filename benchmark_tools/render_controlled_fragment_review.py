"""Apply the established profile-table print fix to the fresh fragment review."""

import argparse
import json
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools import render_profile_stratum_review as profile


def run(root, original_assets, output, assets):
    base = Path(root).resolve() / "benchmark_tools/results"
    original_assets = Path(original_assets).resolve()
    output, assets = Path(output).absolute(), Path(assets).absolute()
    profile.require(output != assets and output.parent.resolve() == assets.parent.resolve() == base,
                    "Preserve manuscript-relative assets")
    for path in (output, assets):
        if path.exists() or path.is_symlink():
            raise FileExistsError(path)
    receipt_ref = record(original_assets)
    receipt = json.loads(original_assets.read_text())
    profile.require(receipt["status"] == "manuscript_review_rendered"
                    and receipt["publication_ready"] is False
                    and Path(receipt["html"]["path"]).parent == base, "Wrong original review")
    source = record(__file__)
    styling_source = record(profile.__file__)
    checked = [receipt_ref, receipt["html"], *receipt["sources"], *receipt["targets"], source, styling_source]
    for ref in checked:
        check(ref)
    html = profile.styled(Path(receipt["html"]["path"]).read_text())
    with output.open("x") as stream:
        stream.write(html)
    revised = dict(receipt, html=record(output),
        sources=[*receipt["sources"], receipt_ref, receipt["html"], source, styling_source],
        render_command_scope="Inherited executed Pandoc invocation; only established profile-table print CSS inserted",
        rendering_recovery=dict(operation="insert_scoped_print_css_only", original_assets=receipt_ref,
            original_html=receipt["html"], source=source, styling_source=styling_source,
            scientific_html_unchanged=True, profile_rows=23, all_columns_required=6, selector=profile.SELECTOR,
            fragment_tables_changed=False))
    for ref in checked:
        check(ref)
    with assets.open("x") as stream:
        json.dump(revised, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return revised


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "original-assets", "output", "assets"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.root, args.original_assets, args.output, args.assets)))
