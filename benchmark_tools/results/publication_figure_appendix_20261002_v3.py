import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import fitz

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / 'benchmark_tools/results'
parser = argparse.ArgumentParser(description='Combine retained text/figures for working review, not publication readiness.')
parser.add_argument('--output', type=Path, default=RESULTS / 'publication_main_with_figures_20261002_v3')
OUTPUT = parser.parse_args().output.absolute()
MAIN = RESULTS / 'publication_main_print_20261002_v3/document.pdf'
MAIN_SHA = '65e271daee8bd8047543d614d937d26f285a212ab8e3823eea06a5a6fce87b53'
FIGURES = [
    ('figures_publication_method_20260916/publication_method.pdf', 'HMM-centered workflow',
     'Initial search, profile expansion, candidate families and phylogenetic refinement have distinct outputs. This schematic describes the frozen workflow, not proof of an HMM accuracy advantage.'),
    ('figures_current_accuracy_overview_20261002_v2/current_accuracy_overview.pdf', 'Current eight-method accuracy overview',
     'OrthoBench precision/recall and Three Kingdoms pair F1 use different reference scopes. These descriptive comparisons do not form a cross-benchmark ranking; Three Kingdoms is supplementary conserved-family evidence.'),
    ('figures_corrected_qfo_endpoints_20260926/corrected_qfo_endpoints.pdf', 'Individual corrected QfO endpoints',
     'GO and EC are similarity scores; VGNC, SwissTrees and TreeFam-A are F1 endpoints, and FAS measures feature-architecture similarity. Eligibility and denominators differ; this plot is descriptive, not paired uncertainty.'),
    ('figures_ob_complete_uncertainty_20260928/ob_complete_uncertainty.pdf', 'OrthoBench paired differences',
     'All 21 method-metric contrasts retain their multiplicity-adjusted RefOG intervals for F1, precision and recall. Development exposure and family-exchangeability assumptions limit inference. An interval containing zero establishes neither superiority nor equivalence.'),
    ('qfo_corrected_factorial_figures_20260919/qfo_factorial_swiss.pdf', 'HMM, expansion and reconciliation factorial',
     'The corrected SwissTrees factorial retains all cells and component contrasts. Candidate expansion interacts with reconciliation; profile-refinement F1 intervals include zero. These conditional effects do not establish universal improvement.'),
    ('figures_qfo_scored_pair_decomposition_20260930/qfo_pair_decomposition.pdf', 'GO and EC scored-pair composition',
     'Shared-denominator and exclusive-pair terms reproduce the aggregate differences for all seven comparators against full OrthoFinder. This is arithmetic decomposition of eligible scored sets, not proteome-wide coverage, causal attribution or confidence intervals.'),
    ('matched_graph_figure_v2_20260926/matched_graph.pdf', 'Matched-recall HMM search control',
     'The fixed-graph diagnostic compares HMM-derived search evidence with DIAMOND on 35 simulation datasets. It is development-exposed, without profile expansion or phylogeny, and does not match computational effort or compare full OrthoFinder.'),
    ('figures_frozen_null_scores_20261002_v2/frozen_null_scores.pdf', 'Synthetic null-score tails',
     'Ninety endpoints arise from 90,000 independent sequence pairs and 180,000 dependent full/banded evaluations. Composition and length change approximate tail behavior. This is not calibrated orthology confidence, a biological false-positive rate or fitted parameters.'),
    ('figures_simulation_fixed_native_20260916/simulation_evidence.pdf', 'Fixed-length evolutionary stress',
     'Native failures remain in the fixed-length panel. No OrthoFinder output passes the stated native-completion/finite-graph gate, so its 14 planned comparisons are unavailable rather than OrthoHMM wins.'),
    ('figures_simulation_variable_native_v2_20260916/simulation_evidence.pdf', 'Heterogeneous-length evolutionary stress',
     'The heterogeneous-length panel retains admitted and failed runs. Full OrthoFinder has higher paired mean F1 than both OrthoHMM modes in all seven conditions. Effects are conditioned on paired success, not failure-adjusted population estimates.'),
    ('figures_simulation_tree_robustness_20260917/simulation_tree_robustness.pdf', 'Species-tree controls and perturbations',
     'Generating and perturbed trees are diagnostic interventions, not independent inferred pipelines or empirical tree posteriors. All 560 outcomes are retained. Two upstream OrthoFinder discrepancies limit strict tree-only causal attribution.'),
    ('qfo_parameter_complete_export_20261001/qfo_parameter_neighborhood.pdf', 'Complete local parameter neighborhood',
     'All seven arms, including separately admitted private-runtime recovery, remain in the 18-endpoint adjustment denominator. Every adjusted interval includes zero; no default or scientific setting was promoted.'),
    ('figures_ygob_overlap_20260928/ygob_overlap_strata.pdf', 'YGOB transfer and overlap strata',
     'The frozen separate-clade test supplies bounded transfer evidence. Its subsequent overlap partition is descriptive, preserves false-positive allocations and does not establish family-disjoint confirmation or a second independent test.'),
    ('figures_wgd_application_20260917/wgd_application.pdf', 'Whole-genome-duplicate biological application',
     'Anchor separation is distinguished from supported separation and reference-homolog coverage. Coverage can be high for merged anchors and is not orthology recall. Prespecified losses and negative OrthoHMM comparisons are retained.')]


def record(path):
    path = Path(path).resolve()
    content = path.read_bytes()
    return dict(path=str(path), bytes=len(content), sha256=hashlib.sha256(content).hexdigest())


def save(path, data):
    with path.open('x') as stream:
        json.dump(data, stream, indent=2, sort_keys=True)
        stream.write('\n')


def box(page, rect, text, size, font='helv'):
    assert page.insert_textbox(rect, text, fontsize=size, fontname=font) >= 0, text


assert not OUTPUT.exists()
fitz.TOOLS.mupdf_warnings(reset=True)
commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
driver_ref = record(Path(__file__))
main_ref = record(MAIN)
assert main_ref['sha256'] == MAIN_SHA
inputs, figures = [main_ref, driver_ref], []
selected_sources = [RESULTS / relative for relative, _, _ in FIGURES]
selected_sources += [path.parent / 'manifest.json' for path in selected_sources]
tracked = set(subprocess.check_output(['git', 'ls-files', '--',
    *[str(path.relative_to(ROOT)) for path in selected_sources]], cwd=ROOT, text=True).splitlines())
local_only = {str(path.relative_to(ROOT)) for path in selected_sources} - tracked
assert local_only == {
    'benchmark_tools/results/figures_simulation_fixed_native_20260916/simulation_evidence.pdf',
    'benchmark_tools/results/figures_simulation_fixed_native_20260916/manifest.json',
    'benchmark_tools/results/qfo_parameter_complete_export_20261001/manifest.json'}
for index, (relative, title, caption) in enumerate(FIGURES, 1):
    path = RESULTS / relative
    ref = record(path)
    manifest_path = path.parent / 'manifest.json'
    manifest_ref = record(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    candidates = [row for row in manifest['outputs'] if Path(row['path']).name == path.name]
    assert len(candidates) == 1
    assert (candidates[0]['bytes'], candidates[0]['sha256']) == (ref['bytes'], ref['sha256'])
    for source in (path, manifest_path):
        relative_source = str(source.relative_to(ROOT))
        if relative_source in tracked:
            blob = subprocess.check_output(['git', 'show', commit + ':' + relative_source], cwd=ROOT)
            assert source.read_bytes() == blob
    if manifest.get('source_results'):
        source_result = manifest['source_results']
        assert record(source_result['path']) == source_result
        inputs.append(source_result)
    if relative.startswith('qfo_parameter_complete_export_20261001/'):
        readback_path = path.parent / 'readback.json'
        readback = json.loads(readback_path.read_text())
        assert readback['manifest'] == manifest_ref and ref in readback['outputs']
        assert readback_path.read_bytes() == subprocess.check_output(
            ['git', 'show', commit + ':' + str(readback_path.relative_to(ROOT))], cwd=ROOT)
        inputs.append(record(readback_path))
    inputs.extend([ref, manifest_ref])
    figures.append(dict(number=index, title=title, caption=caption, pdf=ref, producer_manifest=manifest_ref))
assert MAIN.read_bytes() == subprocess.check_output(['git', 'show', commit + ':' + str(MAIN.relative_to(ROOT))], cwd=ROOT)
OUTPUT.mkdir()
document = fitz.open()
with fitz.open(MAIN) as source:
    assert len(source) == 9
    document.insert_pdf(source)
main_pages = len(document)
guide_pages = (len(figures) + 3) // 4
first_figure = main_pages + guide_pages
guide_links = []
for group in range(guide_pages):
    page = document.new_page(width=595.276, height=841.89)
    box(page, fitz.Rect(40, 38, 555, 72), 'Figure Appendix - Working Review', 18, 'hebo')
    box(page, fitz.Rect(40, 78, 555, 120),
        'Original main text: pages 1-9. Figures retain their original vector pages and dimensions. '
        'Controlled-resource evidence is still missing; no timing figure or publication-readiness claim is added.', 10)
    for slot, item in enumerate(figures[group * 4:(group + 1) * 4]):
        top = 135 + slot * 160
        title = 'A%d. %s' % (item['number'], item['title'])
        box(page, fitz.Rect(40, top, 555, top + 30), title, 12, 'hebo')
        box(page, fitz.Rect(40, top + 35, 555, top + 110), item['caption'], 10.5)
        box(page, fitz.Rect(40, top + 115, 555, top + 145),
            'Source: ' + str(Path(item['pdf']['path']).relative_to(RESULTS)), 8)
        target = first_figure + item['number'] - 1
        guide_links.append((page.number, fitz.Rect(40, top, 555, top + 145), target))
    box(page, fitz.Rect(40, 790, 555, 817), 'Guide %d of %d - click an entry or use PDF bookmarks.' % (group + 1, guide_pages), 9)
for item in figures:
    item['combined_page'] = len(document) + 1
    with fitz.open(item['pdf']['path']) as source:
        assert len(source) == 1
        document.insert_pdf(source)
for guide_index, rectangle, target in guide_links:
    document[guide_index].insert_link({'kind': fitz.LINK_GOTO, 'from': rectangle, 'page': target})
destinations = {item['pdf']['path']: item['combined_page'] - 1 for item in figures}
restored_file_actions = []
with fitz.open(MAIN) as source_main:
    for index in range(main_pages):
        page = document[index]
        for original_link in source_main[index].get_links():
            if original_link.get('file') in destinations or original_link['kind'] != fitz.LINK_LAUNCH:
                continue
            uri_type, uri = source_main.xref_get_key(original_link['xref'], 'A/URI')
            assert uri_type == 'string' and uri.startswith('file:')
            copied = [link for link in page.get_links() if link['from'] == original_link['from']]
            assert len(copied) == 1
            page.update_link({'xref': copied[0]['xref'], 'kind': fitz.LINK_URI,
                              'from': original_link['from'], 'uri': uri})
            restored_file_actions.append(dict(main_page=index + 1, uri=uri,
                original_kind=original_link['kind'], copied_kind=copied[0]['kind']))
assert restored_file_actions
redirected = []
for index in range(main_pages):
    page = document[index]
    for link in page.get_links():
        filename = link.get('file')
        if filename in destinations:
            old = dict(link)
            rectangle = link['from']
            page.update_link({'xref': link['xref'], 'kind': fitz.LINK_GOTO, 'from': rectangle,
                              'page': destinations[filename]})
            redirected.append(dict(main_page=index + 1, source=filename,
                                   combined_page=destinations[filename] + 1, original_kind=old['kind']))
assert len(redirected) == 6
document.set_toc([[1, 'Main Text', 1], [1, 'Figure Guide', main_pages + 1],
                 *[[1, 'A%d. %s' % (item['number'], item['title']), item['combined_page']] for item in figures]])
document.set_metadata(dict(title='OrthoHMM Working Manuscript With Figure Appendix',
                           subject='Review copy; controlled timing and publication readiness incomplete'))
output_pdf = OUTPUT / 'document.pdf'
document.save(output_pdf, garbage=3, deflate=True)
document.close()
page_checks, renders = [], []
with fitz.open(output_pdf) as combined:
    assert len(combined) == main_pages + guide_pages + len(figures)
    source_pages = [(MAIN, page, page) for page in range(main_pages)] + [
        (Path(item['pdf']['path']), 0, item['combined_page'] - 1) for item in figures]
    for path, original_index, combined_index in source_pages:
        with fitz.open(path) as original:
            before, after = original[original_index], combined[combined_index]
            assert before.rect == after.rect
            assert before.get_text('words') == after.get_text('words')
            a, b = before.get_pixmap(alpha=False), after.get_pixmap(alpha=False)
            assert (a.width, a.height, a.n) == (b.width, b.height, b.n)
            assert a.samples == b.samples
            page_checks.append(dict(source=record(path), original_page=original_index + 1,
                combined_page=combined_index + 1, pixel_sha256=hashlib.sha256(a.samples).hexdigest(),
                exact_text_geometry_and_pixels=True))
    for index in range(main_pages, len(combined)):
        image = OUTPUT / ('page_%02d.png' % (index + 1))
        combined[index].get_pixmap(matrix=fitz.Matrix(1.25, 1.25), alpha=False).save(image)
        renders.append(record(image))
    links = [[dict(kind=link['kind'], page=link.get('page')) for link in page.get_links()
              if link['kind'] == fitz.LINK_GOTO] for page in combined]
    assert sum(len(links[index]) for index in range(main_pages, first_figure)) == len(figures)
    assert all(0 <= link['page'] < len(combined) for row in links for link in row)
    outlines = combined.get_toc()
    assert len(outlines) == len(figures) + 2
for ref in inputs:
    assert record(ref['path']) == ref
backend_warnings = fitz.TOOLS.mupdf_warnings(reset=True)
assert not backend_warnings, backend_warnings
report = dict(status='working_manuscript_and_existing_figures_combined', source_commit=commit,
    source=driver_ref, pdf=record(output_pdf), inputs=inputs, figures=figures,
    main_pages=main_pages, guide_pages=guide_pages, figure_pages=len(figures), total_pages=27,
    redirected_main_figure_links=redirected, restored_original_file_uri_actions=restored_file_actions, bookmarks=outlines,
    preserved_source_pages=page_checks, rendered_pages=renders, pymupdf_version=fitz.__version__,
    local_only_direct_sources=sorted(local_only),
    initial_attempt_failure='Pre-assembly Git check failed on the retained local-only fixed-length figure. No output directory existed. Final checks distinguish exactly three local-only sources and bind producer/output/result or committed compact-readback evidence.',
    mupdf_warnings=backend_warnings,
    previous_assembly=record(RESULTS / 'publication_main_with_figures_20261002_v2/assembly.json'),
    previous_version_not_accepted='V1 guide A2 incorrectly named an absent QfO mean and emitted six transient xref warnings. V2 corrected those but its separate readback caught lossy local-file link actions from page copying. Preserve both attempts; V3 restores original file URI actions and redirects only the six included figure links.',
    original_text_and_figure_pages_pixel_identical=True, visual_review_complete=False,
    scientific_settings_or_results_changed=False, scientific_timings_admitted=False, publication_ready=False,
    limitations=[
        'Existing vector PDFs combined, not scientific plotting, statistics or native inference recomputed.',
        'All direct producer PDF hashes checked; 25 of 28 figure/manifest files and main PDF match current Git. Three explicit local-only sources use retained producer/result or committed compact readback bindings.',
        'Available direct source-result file hashes checked, not historical transitive provenance or statistics recomputed.',
        'Figure pages retain differing original sizes; this is a working review, not journal typesetting.',
        'Audit/data/code links still need external files; self-contained for included text/figures, not executable release.',
        'Controlled timing, final manuscript reconciliation and complete publication requirements remain missing.'])
save(OUTPUT / 'assembly.json', report)
print(json.dumps(dict(pdf=report['pdf'], assembly=record(OUTPUT / 'assembly.json'),
                     pages=27, preserved_pages=len(page_checks), redirected_links=len(redirected)), indent=2))
