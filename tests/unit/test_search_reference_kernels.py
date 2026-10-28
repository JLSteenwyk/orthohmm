"""Compare interpreted and compiled search kernels with independent oracles."""

from collections import Counter, defaultdict
import math

import numpy as np
import pytest

from orthohmm.search import evalue, prefilter, profile_spans, viterbi


def interpreted(function):
    return getattr(function, "py_func", function)


def packed(sequences):
    lengths = np.array([len(s) for s in sequences], dtype=np.int32)
    offsets = np.zeros(len(sequences), dtype=np.int64)
    offsets[1:] = np.cumsum(lengths[:-1])
    flat = np.array([x for s in sequences for x in s], dtype=np.uint8)
    return flat, offsets, lengths


@pytest.mark.parametrize("score", [-5, 0, 1, 50, 1000, 10000])
@pytest.mark.parametrize("query_length,total_residues", [(1, 1), (100, 10000)])
def test_evalue_python_and_jit_match_analytical_formula(score, query_length, total_residues):
    expected = (1e10 if score <= 0 else
                .041 * query_length * total_residues * math.exp(-.267 * score))
    arguments = score, query_length, total_residues, .267, .041
    assert interpreted(evalue.estimate_evalue)(*arguments) == pytest.approx(expected, rel=1e-12, abs=0)
    assert evalue.estimate_evalue(*arguments) == pytest.approx(expected, rel=1e-12, abs=0)


@pytest.mark.parametrize("scores", [[], [-1, 0, 1, 50, 1000]])
def test_batch_evalues_match_formula_without_mutating_inputs(scores):
    values = np.array(scores, dtype=np.int32)
    lengths = np.arange(1, len(scores) + 1, dtype=np.int32)
    before = values.tobytes(), lengths.tobytes()
    expected = np.array([1e10 if score <= 0 else .041 * length * 1000 * math.exp(-.267 * score)
                         for score, length in zip(scores, lengths)], dtype=np.float64)
    for function in (interpreted(evalue.batch_estimate_evalues), evalue.batch_estimate_evalues):
        result = function(values, lengths, 1000, .267, .041)
        np.testing.assert_allclose(result, expected, rtol=1e-12, atol=0)
        assert result.dtype == np.float64
    assert before == (values.tobytes(), lengths.tobytes())


def full_matrix_score(emissions, insertions, transitions, target, band):
    """Full storage oracle: each cell reads its explicit three predecessors."""
    length, size = len(emissions), len(target)
    states = np.full((3, length + 1, size + 1), viterbi.NEG_INF, dtype=np.int64)
    width = band if band > 0 and max(length, size) > 50 else max(length, size)
    best = 0
    for i in range(1, length + 1):
        for j in range(1, size + 1):
            if abs(j - i * size // length) > width:
                continue
            aa = int(target[j - 1])
            emission = int(emissions[i - 1, aa]) if aa < 20 else 0
            insertion = int(insertions[aa]) if aa < 20 else -1
            states[0, i, j] = emission + max(0,
                states[0, i - 1, j - 1] + transitions[0],
                states[1, i - 1, j - 1] + transitions[3],
                states[2, i - 1, j - 1] + transitions[5])
            previous = max(states[0, i, j - 1] + transitions[1],
                           states[1, i, j - 1] + transitions[4])
            states[1, i, j] = previous + insertion if previous > viterbi.NEG_INF // 2 else viterbi.NEG_INF
            states[2, i, j] = max(states[0, i - 1, j] + transitions[2],
                                 states[2, i - 1, j] + transitions[6])
            best = max(best, int(states[0, i, j]))
    return best


def model(symbols):
    emissions = np.full((len(symbols), 20), -9, dtype=np.int8)
    for i, aa in enumerate(symbols):
        emissions[i, aa] = 5
    return (emissions, np.full(20, -1, dtype=np.int8),
            np.array([0, -1, -1, 0, 0, 0, 0], dtype=np.int32))


@pytest.mark.parametrize("symbols,target", [
    ([], []), ([], [0]), ([0], []), ([0], [20]), ([0], [255]),
    ([0, 1], [0, 1]), ([0, 1], [0, 2, 2, 2, 1]),
    ([0, 2, 2, 2, 1], [0, 1]), ([0, 1, 0], [2, 0, 1, 0, 2]),
    ([0, 1, 0], [0, 20, 0]), ([0] * 51, [0] * 60),
])
@pytest.mark.parametrize("band", [0, 1, 8])
def test_viterbi_and_span_python_jit_match_full_matrix_oracle(symbols, target, band):
    emissions, insertions, transitions = model(symbols)
    sequence = np.array(target, dtype=np.uint8)
    before = [a.tobytes() for a in (emissions, insertions, transitions, sequence)]
    expected = full_matrix_score(emissions, insertions, transitions, sequence, band)
    arguments = emissions, insertions, transitions, sequence, len(target), len(symbols), band
    assert interpreted(viterbi._viterbi_score_one)(*arguments) == expected
    assert viterbi._viterbi_score_one(*arguments) == expected
    span_arguments = emissions, insertions, transitions, sequence, band
    python_span = interpreted(profile_spans._span_one)(*span_arguments)
    native_span = profile_spans._span_one(*span_arguments)
    np.testing.assert_array_equal(python_span, native_span)
    assert python_span[0] == expected
    assert [a.tobytes() for a in (emissions, insertions, transitions, sequence)] == before


def test_viterbi_batch_offsets_order_duplicates_and_empty_pairs(monkeypatch):
    models = [model([])[0], model([0, 1])[0], model([1, 0, 1, 0])[0]]
    lengths = np.array([len(m) for m in models], dtype=np.int32)
    offsets = np.array([0, 0, 2], dtype=np.int64)
    emissions = np.concatenate(models)
    _, insertions, transitions = model([])
    sequences = [[], [0], [0, 1, 0], [1, 0, 1, 0, 1]]
    flat, target_offsets, target_lengths = packed(sequences)
    for rows in ([], [[2, 3], [0, 1], [1, 2], [2, 3], [1, 0]]):
        pairs = np.array(rows, dtype=np.int32).reshape(-1, 2)
        expected = [full_matrix_score(models[q], insertions, transitions, sequences[t], 0) for q, t in rows]
        arguments = (emissions, insertions, transitions, offsets, lengths,
                     flat, target_offsets, target_lengths, pairs, 0)
        with monkeypatch.context() as patch:
            patch.setattr(viterbi, "_viterbi_score_one", interpreted(viterbi._viterbi_score_one))
            result = interpreted(viterbi.batch_viterbi_score)(*arguments)
        np.testing.assert_array_equal(result, expected)
        np.testing.assert_array_equal(viterbi.batch_viterbi_score(*arguments), expected)
    indices = np.array([3, 1, 3, 0], dtype=np.int32)
    arguments = (models[1], insertions, transitions, flat, target_offsets, target_lengths, indices, 2, 0)
    expected = [full_matrix_score(models[1], insertions, transitions, sequences[t], 0) for t in indices]
    with monkeypatch.context() as patch:
        patch.setattr(viterbi, "_viterbi_score_one", interpreted(viterbi._viterbi_score_one))
        result = interpreted(viterbi.viterbi_all_targets)(*arguments)
    np.testing.assert_array_equal(result, expected)
    np.testing.assert_array_equal(viterbi.viterbi_all_targets(*arguments), expected)


TARGETS = [[], [0], [0, 1, 0, 1], [1, 0, 1, 0], [2, 2, 2], [0, 1, 20, 0, 1]]


def postings(sequences, k=2, alphabet=3):
    result = defaultdict(list)
    for index, sequence in enumerate(sequences):
        for position in range(len(sequence) - k + 1):
            word = sequence[position:position + k]
            if all(aa < alphabet for aa in word):
                code = 0
                for aa in word:
                    code = code * alphabet + aa
                result[tuple(word)].append((index, position, code))
    return result


def test_python_and_jit_kmer_index_match_independent_postings(monkeypatch):
    flat, offsets, lengths = packed(TARGETS)
    expected = postings(TARGETS)
    for python in (True, False):
        with monkeypatch.context() as patch:
            if python:
                patch.setattr(prefilter, "_count_kmers", interpreted(prefilter._count_kmers))
                patch.setattr(prefilter, "_fill_entries", interpreted(prefilter._fill_entries))
            pointers, entries, frequencies = prefilter.build_kmer_index(flat, offsets, lengths, k=2, alpha_size=3)
        assert pointers.shape == (10,) and frequencies.shape == (9,)
        assert pointers[0] == 0 and pointers[-1] == len(entries)
        for code in range(9):
            wanted = [(index << 32) | pos for values in expected.values()
                      for index, pos, value in values if value == code]
            assert entries[pointers[code]:pointers[code + 1]].tolist() == wanted
            assert frequencies[code] == len(wanted)


def candidate_oracle(query, targets, frequency, minimum, diagonal, bin_width):
    index = postings(targets)
    totals, bins = Counter(), defaultdict(Counter)
    for pos in range(len(query) - 1):
        values = index.get(tuple(query[pos:pos + 2]), [])
        if len(values) > frequency:
            continue
        for target, target_pos, _ in values:
            totals[target] += 1
            bins[target][(target_pos - pos + len(query)) // bin_width] += 1
    result = {}
    for target, count in sorted(totals.items()):
        if count >= minimum:
            score = min(count, 32767) if diagonal <= 1 else max(bins[target].values())
            if diagonal <= 1 or score >= diagonal:
                result[target] = score
    return result


@pytest.mark.parametrize("query", [[], [0], [0, 1, 0, 1], [20, 0, 1, 20, 0, 1], [2, 2, 2]])
@pytest.mark.parametrize("frequency,minimum,diagonal,width", [(0, 1, 1, 1), (3, 1, 1, 1),
    (100, 1, 1, 1), (100, 3, 2, 1), (100, 1, 3, 2), (100, 100, 2, 1)])
def test_two_stage_python_and_jit_candidates_match_diagonal_oracle(query, frequency, minimum, diagonal, width):
    flat, offsets, lengths = packed(TARGETS)
    index = prefilter.build_kmer_index(flat, offsets, lengths, k=2, alpha_size=3)
    arguments = (np.array(query, dtype=np.uint8), len(query), 2, 3, *index,
                 frequency, len(TARGETS), minimum, diagonal, width)
    expected = candidate_oracle(query, TARGETS, frequency, minimum, diagonal, width)
    for function in (interpreted(prefilter._score_query_two_stage), prefilter._score_query_two_stage):
        ids, scores = function(*arguments)
        assert dict(zip(ids.tolist(), scores.tolist())) == expected
        assert ids.dtype == scores.dtype == np.int32


def test_repeated_kmer_counts_saturate_without_integer_wraparound():
    sequences = [[0] * 183]
    flat, offsets, lengths = packed(sequences)
    index = prefilter.build_kmer_index(flat, offsets, lengths, k=2, alpha_size=3)
    arguments = flat, len(flat), 2, 3, *index, 1000, 1, 1, 1, 1
    for function in (interpreted(prefilter._score_query_two_stage), prefilter._score_query_two_stage):
        ids, scores = function(*arguments)
        assert ids.tolist() == [0] and scores.tolist() == [32767]


@pytest.mark.parametrize("cap", [0, 1, 4])
def test_batch_prefilter_csr_matches_per_query_oracle(monkeypatch, cap):
    queries = [[], [0], [0, 1, 0, 1], [2, 2, 2], [20, 20]]
    query_flat, query_offsets, query_lengths = packed(queries)
    flat, offsets, lengths = packed(TARGETS)
    index = prefilter.build_kmer_index(flat, offsets, lengths, k=2, alpha_size=3)
    arguments = (query_flat, query_offsets, query_lengths, 2, 3, *index,
                 100, len(TARGETS), 1, 2, 1, cap)
    expected_ids, expected_offsets = [], [0]
    for query in queries:
        expected_ids.extend(list(candidate_oracle(query, TARGETS, 100, 1, 2, 1))[:cap])
        expected_offsets.append(len(expected_ids))
    with monkeypatch.context() as patch:
        patch.setattr(prefilter, "_score_query_two_stage", interpreted(prefilter._score_query_two_stage))
        ids, pointers = interpreted(prefilter.batch_prefilter)(*arguments)
    np.testing.assert_array_equal(ids, expected_ids)
    np.testing.assert_array_equal(pointers, expected_offsets)
    compiled_ids, compiled_offsets = prefilter.batch_prefilter(*arguments)
    np.testing.assert_array_equal(compiled_ids, expected_ids)
    np.testing.assert_array_equal(compiled_offsets, expected_offsets)


def test_reduced_alphabet_python_and_jit_preserve_unknown_symbols():
    sequence = np.array([0, 1, 19, 20, 255], dtype=np.uint8)
    expected = [int(prefilter.REDUCED_ALPHA[i]) if i < 20 else 20 for i in sequence]
    for function in (interpreted(prefilter._remap_flat), prefilter._remap_flat):
        np.testing.assert_array_equal(function(sequence, prefilter.REDUCED_ALPHA, len(sequence)), expected)
