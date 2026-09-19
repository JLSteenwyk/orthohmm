"""Exploratory native-output comparison; neither inferred tree is truth."""

import argparse
import hashlib
import json
from pathlib import Path


def load(path, expected_sha):
    data = Path(path).read_bytes()
    if hashlib.sha256(data).hexdigest() != expected_sha:
        raise ValueError('Input checksum differs: ' + str(path))
    return json.loads(data)


def compare_trees(left, right):
    import dendropy
    from dendropy.calculate import treecompare

    namespace = dendropy.TaxonNamespace()
    trees = [dendropy.Tree.get(path=str(path), schema='newick', taxon_namespace=namespace,
                              preserve_underscores=True, rooting='force-rooted') for path in (left, right)]
    labels = [{leaf.taxon.label for leaf in tree.leaf_node_iter()} for tree in trees]
    if labels[0] != labels[1]:
        raise ValueError('Species universes differ')
    for tree in trees:
        tree.encode_bipartitions()
    clades = [{frozenset(leaf.taxon.label for leaf in node.leaf_iter())
               for node in tree.postorder_node_iter()
               if node.parent_node is not None and not node.is_leaf()} for tree in trees]
    distance = int(treecompare.symmetric_difference(*trees))
    if distance != len(clades[0] ^ clades[1]):
        raise ValueError('Library distance differs from explicit rooted clade comparison')
    return dict(species=len(labels[0]), rooted_symmetric_difference=distance,
                left_clades=len(clades[0]), right_clades=len(clades[1]),
                shared_clades=len(clades[0] & clades[1]),
                left_only_clades=len(clades[0] - clades[1]), right_only_clades=len(clades[1] - clades[0]),
                dendropy_version=dendropy.__version__, either_tree_is_ground_truth=False)


def compare(left, right):
    reports = [load(*item) for item in (left, right)]
    for report in reports:
        if report['status'] != 'corrected_qfo_native_pair_output_verified':
            raise ValueError('Require admitted corrected native output')
    for key in ('candidate_admission', 'environment', 'prepared'):
        if reports[0][key]['sha256'] != reports[1][key]['sha256']:
            raise ValueError('Unmatched native provenance: ' + key)
    rows, tree_paths = [], []
    for report in reports:
        integrity = report['native_group_integrity']
        metrics = load(integrity['native_metrics']['path'], integrity['native_metrics']['sha256'])
        if (metrics['counts']['ortholog_pairs'] != report['native_pair_count']
                or metrics['counts']['candidate_families'] != integrity['partition']['candidate_families']
                or metrics['counts']['root_hogs'] != integrity['partition']['root_hogs']):
            raise ValueError('Native metric counts differ from admission')
        tree = integrity['species_tree']
        data = Path(tree['path']).read_bytes()
        if len(data) != tree['bytes'] or hashlib.sha256(data).hexdigest() != tree['sha256']:
            raise ValueError('Native tree changed')
        tree_paths.append(tree['path'])
        rows.append(dict(cell=report['cell'], native_pairs=report['native_pair_count'],
                         partition=integrity['partition'], membership=integrity['membership'],
                         counts=metrics['counts'], species_tree=tree, native_metrics=integrity['native_metrics']))
    return dict(status='exploratory_corrected_native_comparison',
                inputs=[dict(path=str(p), sha256=s) for p,s in (left,right)], rows=rows,
                tree_comparison=compare_trees(*tree_paths), accuracy_evaluated=False,
                publication_ready=False,
                limitations=['Post-outcome descriptive diagnostic, not a prespecified accuracy endpoint.',
                    'Different inferred trees do not establish which is correct or explain an accuracy change causally.',
                    'Candidate expansion contrast includes downstream tree estimation; it is not a fixed-tree isolation.',
                    'Native pair/group/event counts are not biological accuracy or reference coverage.'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('left', 'right'):
        parser.add_argument('--' + name, nargs=2, required=True, metavar=('ADMISSION', 'SHA256'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = compare(args.left, args.right)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write('\n')
