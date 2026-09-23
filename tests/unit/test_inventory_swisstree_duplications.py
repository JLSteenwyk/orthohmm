import pytest

from benchmark_tools.inventory_swisstree_duplications import parse, FAMILIES


def table():
    return ("<table><tr><th>ID</th><th>Family</th><th>Duplications</th></tr>" +
        "".join(f"<tr><td>{name}</td><td><a>Family*</a></td><td>2</td></tr>" for name in FAMILIES) + "</table>")


def test_complete_inventory_preserves_updated_flags_and_excluded_family():
    result = parse(table())
    assert len(result) == 19
    assert result["ST012"]["qfo_family"] is None
    assert result["ST001"]["updated_marker"] is True
    assert result["ST001"]["published_duplications"] == 2


@pytest.mark.parametrize("problem", ["duplicate", "missing", "unknown", "negative", "absent", "header"])
def test_invalid_source_rejected(problem):
    text = table()
    if problem == "duplicate":
        text = text.replace("ST002", "ST001")
    elif problem == "missing":
        text = text.replace("<td>ST002</td>", "<td>ignored</td>")
    elif problem == "unknown":
        text = text.replace("ST002", "ST020")
    elif problem in {"negative", "absent"}:
        text = text.replace("<td>2</td>", "<td>-1</td>" if problem == "negative" else "<td>-</td>")
    else:
        text = text.replace("Duplications", "Other")
    with pytest.raises(ValueError):
        parse(text)
