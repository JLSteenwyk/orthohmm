import pytest

from benchmark_tools.export_qfo_threshold_endpoints import scores, AXES


def assessment():
    records, endpoints = [], {}
    for challenge, axes in {**AXES, "SwissTrees-family": ("TPR", "PPV")}.items():
        x, y = (10, .8) if axes[0] == "NR_ORTHOLOGS" else (.5, .8)
        for metric, value in zip(axes, (x, y)):
            records.append(dict(type="assessment", community_id="QfO", participant_id="fixture",
                _id=challenge+metric, challenge_id=challenge,
                metrics=dict(metric_id=metric, value=value, stderr=0)))
        if challenge in AXES:
            endpoints[challenge] = dict(native_participant=dict(participant_id="fixture", metric_x=x, metric_y=y),
                axes=dict(x_axis=axes[0], y_axis=axes[1]), score=y if x == 10 else 2*x*y/(x+y))
    return dict(participant="fixture", swiss_reference_families=["family"], native_assessments=records,
        endpoints=endpoints, secondary_six_metric_mean=sum(r["score"] for r in endpoints.values())/6)


def test_distinct_endpoint_semantics():
    result = scores(assessment())
    assert result["GO"] == result["EC"] == result["FAS"] == .8
    assert result["SwissTrees"] == pytest.approx(2*.5*.8/1.3)


@pytest.mark.parametrize("change", ["participant", "axis", "score", "mean", "missing", "native", "count"])
def test_corrupted_endpoint_rejected(change):
    value = assessment()
    if change == "participant":
        value["endpoints"]["GO"]["native_participant"]["participant_id"] = "other"
    elif change == "axis":
        value["endpoints"]["GO"]["axes"]["x_axis"] = "TPR"
    elif change == "score":
        value["endpoints"]["GO"]["score"] = .4
    elif change == "mean":
        value["secondary_six_metric_mean"] = .4
    elif change == "missing":
        value["endpoints"].pop("EC")
    elif change == "native":
        value["native_assessments"].pop()
    elif change == "count":
        value["endpoints"]["GO"]["native_participant"]["metric_x"] = 11
    with pytest.raises(ValueError):
        scores(value)
