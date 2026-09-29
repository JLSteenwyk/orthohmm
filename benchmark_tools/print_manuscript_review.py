"""Print a verified HTML review into a new, non-overwriting attempt directory."""

import argparse
import json
from pathlib import Path
import subprocess
import time

import fitz

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def print_review(assets, output, browser, no_sandbox=False):
    assets, browser = assets.resolve(), browser.resolve()
    # Do not resolve the final component: a dangling symlink is occupied too.
    output = output.parent.resolve() / output.name
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    receipt_record = record(assets)
    receipt = json.loads(assets.read_text())
    if receipt.get("status") != "manuscript_review_rendered":
        raise ValueError("Require a successful HTML render receipt")
    refs = [receipt_record, receipt["html"], *receipt["sources"],
            *receipt["targets"], record(browser), record(Path(__file__))]
    for ref in refs:
        check(ref)
    html = Path(receipt["html"]["path"])
    if html.suffix.lower() != ".html":
        raise ValueError("Require an HTML source")
    output.mkdir(parents=True, exist_ok=False)
    pdf = output / "document.pdf"
    command = [str(browser), "--headless", "--disable-gpu",
               "--no-pdf-header-footer", f"--user-data-dir={output / 'browser_profile'}",
               f"--print-to-pdf={pdf}"]
    if no_sandbox:
        command.append("--no-sandbox")
    command.append(html.as_uri())
    result = dict(status="print_started", command=command, checked_records=refs,
                  publication_ready=False, visual_review_complete=False,
                  limitations=[
                      "Printing and PDF parsing do not establish scientific correctness or visual fidelity.",
                      "Direct render inputs are checked; transitive links and network assets are not audited.",
                      "Browser launcher identity is recorded, not its complete runtime dependency closure."])
    started = time.monotonic()
    try:
        process = subprocess.run(command, capture_output=True, text=True, timeout=120)
        result.update(returncode=process.returncode, stdout=process.stdout, stderr=process.stderr)
        process.check_returncode()
        for ref in refs:
            check(ref)
        with fitz.open(pdf) as document:
            if not document.is_pdf or document.page_count < 1:
                raise ValueError("Browser did not produce a nonempty PDF")
            result["page_count"] = document.page_count
        result.update(status="verified_html_printed", pdf=record(pdf))
    except Exception as error:
        result.update(status="print_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        result["elapsed_seconds"] = time.monotonic() - started
        with (output / "print.json").open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("assets", "output", "browser"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--no-sandbox", action="store_true",
                        help="Explicitly disable the browser sandbox for trusted local input")
    args = parser.parse_args()
    result = print_review(args.assets, args.output, args.browser, args.no_sandbox)
    print(json.dumps({key: result[key] for key in ("status", "page_count", "pdf")}))
