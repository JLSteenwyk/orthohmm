"""Rebuild the review presentation from committed PDFs, not scientific analyses."""
import argparse
import hashlib
import json
from pathlib import Path

import fitz

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / 'benchmark_tools/results'
parser = argparse.ArgumentParser(description='Rebuild the review PDF from committed inputs; not inference or publication readiness.')
parser.add_argument('--output', type=Path, default=RESULTS / 'publication_review_rebuild')
OUTPUT = parser.parse_args().output.absolute()
MAIN = RESULTS / 'publication_main_print_20261002_v3/document.pdf'
MAIN_SHA = '65e271daee8bd8047543d614d937d26f285a212ab8e3823eea06a5a6fce87b53'
ASSEMBLY = RESULTS / 'publication_main_with_figures_20261002_v3/assembly.json'
ASSEMBLY_SHA = '861b869a4e359106843692c50a3e43e7805c177b32cce151fcf9cc2ab271176e'


def relocated(ref):
    marker = '/benchmark_tools/results/'
    parts = ref['path'].split(marker)
    if len(parts) != 2:
        raise ValueError('Expected a pinned results path: ' + ref['path'])
    relative = Path(parts[1])
    if relative.is_absolute() or '..' in relative.parts:
        raise ValueError('Unsafe pinned path: ' + ref['path'])
    path = RESULTS / relative
    current = record(path)
    if (current['bytes'], current['sha256']) != (ref['bytes'], ref['sha256']):
        raise ValueError('Input checksum mismatch: ' + str(relative))
    return current


def record(path):
    path = Path(path).resolve()
    content = path.read_bytes()
    return dict(path=str(path), bytes=len(content), sha256=hashlib.sha256(content).hexdigest())


def save(path, data):
    with path.open('x') as stream:
        json.dump(data, stream, indent=2, sort_keys=True)
        stream.write('\n')


def require(condition, message='Presentation check failed'):
    if not condition:
        raise ValueError(message)


def box(page, rect, text, size, font='helv'):
    require(page.insert_textbox(rect, text, fontsize=size, fontname=font) >= 0, text)


if OUTPUT.exists():
    raise FileExistsError('Refusing existing output directory: ' + str(OUTPUT))
fitz.TOOLS.mupdf_warnings(reset=True)
driver_ref = record(Path(__file__))
assembly_ref = record(ASSEMBLY)
if assembly_ref['sha256'] != ASSEMBLY_SHA:
    raise ValueError('Accepted assembly receipt checksum mismatch')
accepted = json.loads(ASSEMBLY.read_text())
main_ref = relocated(accepted['preserved_source_pages'][0]['source'])
if main_ref['sha256'] != MAIN_SHA or Path(main_ref['path']) != MAIN:
    raise ValueError('Unexpected manuscript source')
inputs, figures = [main_ref, driver_ref, assembly_ref], []
for original in accepted['figures']:
    ref = relocated(original['pdf'])
    inputs.append(ref)
    figures.append(dict(number=original['number'], title=original['title'],
        caption=original['caption'], pdf=ref, original_pdf_path=original['pdf']['path']))
if len(figures) != 14 or [item['number'] for item in figures] != list(range(1, 15)):
    raise ValueError('Unexpected figure inventory')
OUTPUT.mkdir()
document = fitz.open()
with fitz.open(MAIN) as source:
    require(len(source) == 9)
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
        require(len(source) == 1)
        document.insert_pdf(source)
for guide_index, rectangle, target in guide_links:
    document[guide_index].insert_link({'kind': fitz.LINK_GOTO, 'from': rectangle, 'page': target})
# The original manuscript embeds historical absolute paths. Match those too.
destinations = {path: item['combined_page'] - 1 for item in figures
                for path in (item['pdf']['path'], item['original_pdf_path'])}
restored_file_actions = []
with fitz.open(MAIN) as source_main:
    for index in range(main_pages):
        page = document[index]
        for original_link in source_main[index].get_links():
            if original_link.get('file') in destinations or original_link['kind'] != fitz.LINK_LAUNCH:
                continue
            uri_type, uri = source_main.xref_get_key(original_link['xref'], 'A/URI')
            require(uri_type == 'string' and uri.startswith('file:'))
            copied = [link for link in page.get_links() if link['from'] == original_link['from']]
            require(len(copied) == 1)
            page.update_link({'xref': copied[0]['xref'], 'kind': fitz.LINK_URI,
                              'from': original_link['from'], 'uri': uri})
            restored_file_actions.append(dict(main_page=index + 1, uri=uri,
                original_kind=original_link['kind'], copied_kind=copied[0]['kind']))
require(restored_file_actions)
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
require(len(redirected) == 6)
document.set_toc([[1, 'Main Text', 1], [1, 'Figure Guide', main_pages + 1],
                 *[[1, 'A%d. %s' % (item['number'], item['title']), item['combined_page']] for item in figures]])
document.set_metadata(dict(title='OrthoHMM Working Manuscript With Figure Appendix',
                           subject='Review copy; controlled timing and publication readiness incomplete'))
output_pdf = OUTPUT / 'document.pdf'
document.save(output_pdf, garbage=3, deflate=True)
document.close()
page_checks, renders = [], []
with fitz.open(output_pdf) as combined:
    require(len(combined) == main_pages + guide_pages + len(figures))
    source_pages = [(MAIN, page, page) for page in range(main_pages)] + [
        (Path(item['pdf']['path']), 0, item['combined_page'] - 1) for item in figures]
    for path, original_index, combined_index in source_pages:
        with fitz.open(path) as original:
            before, after = original[original_index], combined[combined_index]
            require(before.rect == after.rect)
            require(before.get_text('words') == after.get_text('words'))
            a, b = before.get_pixmap(alpha=False), after.get_pixmap(alpha=False)
            require((a.width, a.height, a.n) == (b.width, b.height, b.n))
            require(a.samples == b.samples)
            page_checks.append(dict(source=record(path), original_page=original_index + 1,
                combined_page=combined_index + 1, pixel_sha256=hashlib.sha256(a.samples).hexdigest(),
                exact_text_geometry_and_pixels=True))
    for index in range(main_pages, len(combined)):
        image = OUTPUT / ('page_%02d.png' % (index + 1))
        combined[index].get_pixmap(matrix=fitz.Matrix(1.25, 1.25), alpha=False).save(image)
        renders.append(record(image))
    links = [[dict(kind=link['kind'], page=link.get('page')) for link in page.get_links()
              if link['kind'] == fitz.LINK_GOTO] for page in combined]
    require(sum(len(links[index]) for index in range(main_pages, first_figure)) == len(figures))
    require(all(0 <= link['page'] < len(combined) for row in links for link in row))
    outlines = combined.get_toc()
    require(len(outlines) == len(figures) + 2)
for ref in inputs:
    require(record(ref['path']) == ref)
backend_warnings = fitz.TOOLS.mupdf_warnings(reset=True)
require(not backend_warnings, backend_warnings)
report = dict(status='committed_input_presentation_rebuild', source=driver_ref,
    accepted_assembly=assembly_ref, pdf=record(output_pdf), inputs=inputs, figures=figures,
    main_pages=main_pages, guide_pages=guide_pages, figure_pages=len(figures), total_pages=27,
    redirected_main_figure_links=redirected, restored_original_file_uri_actions=restored_file_actions,
    bookmarks=outlines, preserved_source_pages=page_checks, rendered_pages=renders,
    pymupdf_version=fitz.__version__, mupdf_warnings=backend_warnings,
    original_text_and_figure_pages_pixel_identical=True, visual_review_complete=False,
    scientific_settings_or_results_changed=False, scientific_timings_admitted=False,
    publication_ready=False, limitations=[
        'Presentation replay only: producer manifests, statistical inputs and inference are not replayed.',
        'Hash-pinned committed PDFs and accepted assembly provide captions and presentation inventory.',
        'Original external audit/data/code links retain historical locations and are not portable.',
        'Figure page sizes differ; no journal typesetting or new manual visual acceptance is claimed.',
        'Controlled resources, final reconciliation and complete release remain unfinished.'])

save(OUTPUT / 'assembly.json', report)
print(json.dumps(dict(pdf=report['pdf'], assembly=record(OUTPUT / 'assembly.json'),
                     pages=27, preserved_pages=len(page_checks), redirected_links=len(redirected)), indent=2))
