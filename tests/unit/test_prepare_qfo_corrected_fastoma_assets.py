import pytest

from benchmark_tools.prepare_qfo_corrected_fastoma_assets import IMAGE_ID, IMAGE_DIGEST, image_identity, validate_nextflow


def image():
    return {"Id": IMAGE_ID, "RepoDigests": [IMAGE_DIGEST], "Architecture": "amd64", "Os": "linux",
            "Config": {"Labels": {"org.opencontainers.image.version": "0.3.5"}}, "RootFS": {"Layers": []}}


def test_image_and_nextflow_identity():
    assert image_identity([image()])["Id"] == IMAGE_ID
    validate_nextflow("N E X T F L O W\n version 22.10.8 build 5860\n")


@pytest.mark.parametrize("key,value", [("Id", "other"), ("RepoDigests", []), ("Architecture", "arm64"), ("Os", "other")])
def test_changed_image_rejected(key, value):
    row = image()
    row[key] = value
    with pytest.raises(ValueError, match="identity"):
        image_identity([row])


@pytest.mark.parametrize("rows", [[], [image(), image()]])
def test_ambiguous_image_rejected(rows):
    with pytest.raises(ValueError, match="one"):
        image_identity(rows)


def test_wrong_nextflow_rejected():
    with pytest.raises(ValueError):
        validate_nextflow("version 26.04.6")
