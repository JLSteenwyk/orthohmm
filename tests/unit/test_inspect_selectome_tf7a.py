import pytest

pytest.importorskip("sqlglot")

from benchmark_tools.inspect_selectome_tf7a import literal_rows


def test_literal_rows_preserve_mysql_escaped_text():
    sql = r"INSERT INTO `t` VALUES ('a\'b',1,'(a,b);'),('c',2,'x\\y');"
    assert list(literal_rows(sql, "t")) == [["a'b", 1, "(a,b);"], ["c", 2, "x\\y"]]


@pytest.mark.parametrize("sql", ["SELECT 1", "INSERT INTO other VALUES (1)",
    "INSERT INTO t VALUES (SLEEP(1))", "INSERT INTO t SELECT 1",
    "INSERT INTO t VALUES (1); DROP TABLE t", "INSERT INTO t(a) VALUES (1)"])
def test_nonliteral_or_wrong_statement_rejected(sql):
    with pytest.raises(ValueError):
        list(literal_rows(sql, "t"))
