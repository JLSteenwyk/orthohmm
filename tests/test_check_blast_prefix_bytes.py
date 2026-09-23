import hashlib
import io

import pytest

from benchmark_tools.check_blast_prefix_bytes import scan


def fixture():
    data = b"a\trow\na\trow2\nb\trow\nBROKEN\0"
    cut = len(b"a\trow\na\trow2\n")
    end = cut+len(b"b\trow\n")
    blocks = [dict(query=q, start=a, end=b, rows=n, input_ordinal_0based=i,
                   sha256=hashlib.sha256(data[a:b]).hexdigest(),
                   final_observed_query=final, reuse_authorized=False)
              for q,a,b,n,i,final in [("a",0,cut,2,0,False),("b",cut,end,1,4,True)]]
    boundary = dict(block_count=2,last_query="b",last_query_start=cut,prefix_end=end,excluded_tail_bytes=7)
    return data, blocks, boundary


def test_full_bytes_and_excluded_query():
    data, blocks, boundary = fixture()
    result = scan(io.BytesIO(data), blocks, boundary, hashlib.sha256(data).hexdigest())
    assert result["retained_blocks"] == 1 and result["retained_rows"] == 2
    assert result["excluded_final_complete_rows"] == 1 and result["full_file_bytes"] == len(data)


@pytest.mark.parametrize("problem", ["hash", "gap", "rows", "order", "final", "tail", "truncated", "nul", "boundary"])
def test_invalid_scan(problem):
    data, blocks, boundary = fixture()
    digest = hashlib.sha256(data).hexdigest()
    if problem == "hash": blocks[0]["sha256"] = "0"*64
    elif problem == "gap": blocks[1]["start"] += 1
    elif problem == "rows": blocks[0]["rows"] += 1
    elif problem == "order": blocks[1]["input_ordinal_0based"] = 0
    elif problem == "final": blocks[0]["final_observed_query"] = True
    elif problem == "tail": data += b"x"
    elif problem == "truncated": data = data[:5]
    elif problem == "nul": data = b"\0"+data[1:]
    else: boundary["last_query_start"] -= 1
    with pytest.raises(ValueError):
        scan(io.BytesIO(data), blocks, boundary, digest)
