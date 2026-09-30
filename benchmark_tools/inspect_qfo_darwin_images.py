"""Inspect selected public QfO image build-context layers, never run/extract them."""

import argparse
import hashlib
import io
import json
from pathlib import Path
import re
import tarfile
from urllib.error import HTTPError, URLError
from urllib.request import HTTPRedirectHandler, Request, build_opener


REPOSITORY = "qfobenchmark/darwin"
REGISTRY = f"https://registry-1.docker.io/v2/{REPOSITORY}"
TAGS_URL = f"https://hub.docker.com/v2/repositories/{REPOSITORY}/tags?page_size=100"
ACCEPT = ", ".join((
    "application/vnd.oci.image.index.v1+json",
    "application/vnd.oci.image.manifest.v1+json",
    "application/vnd.docker.distribution.manifest.list.v2+json",
    "application/vnd.docker.distribution.manifest.v2+json",
))


class PublicRedirect(HTTPRedirectHandler):
    def redirect_request(self, request, *args):
        redirected = super().redirect_request(request, *args)
        if redirected is not None:
            redirected.remove_header("Authorization")
        return redirected


def digest(data):
    return "sha256:" + hashlib.sha256(data).hexdigest()


def verified(data, descriptor):
    if len(data) != descriptor["size"] or digest(data) != descriptor["digest"]:
        raise ValueError("Registry blob size/digest mismatch")
    return data


def cached_layer(descriptor, directories):
    if (not re.fullmatch(r"sha256:[0-9a-f]{64}", descriptor["digest"])
            or type(descriptor["size"]) is not int or not 0 <= descriptor["size"] <= 32_000_000):
        raise ValueError("Invalid or oversized context descriptor")
    filename = descriptor["digest"].split(":", 1)[1] + ".tar.gz"
    for directory in directories:
        path = Path(directory) / filename
        if path.is_symlink():
            raise ValueError("Cached context must be a direct file")
        if path.exists():
            if path.stat().st_size != descriptor["size"]:
                raise ValueError("Cached context size differs from registry descriptor")
            data = verified(path.read_bytes(), descriptor)
            return data, dict(path=str(path.resolve()), bytes=len(data), digest=digest(data))
    return None, None


def context_layers(manifest, config):
    history = [item for item in config["history"] if not item.get("empty_layer", False)]
    layers = manifest["layers"]
    if len(history) != len(layers):
        raise ValueError("Nonempty history does not correspond to manifest layers")
    return [(index, layer, item["created_by"])
            for index, (item, layer) in enumerate(zip(history, layers))
            if re.search(r"\b(?:COPY|ADD)\b", item["created_by"])
            and "/benchmark" in item["created_by"]]


def inventory(data):
    members = []
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as archive:
        for member in archive:
            members.append({"name": member.name, "bytes": member.size,
                            "type": member.type.decode("ascii", errors="backslashreplace"),
                            "linkname": member.linkname})
    candidates = [item for item in members if
                  "treefam2reference" in item["name"].lower()
                  or item["name"].lower().endswith((".nhx", ".nhx.gz", ".nhx.bz2"))]
    references = [item for item in members if "treefam" in item["name"].lower()]
    return {"member_count": len(members), "members": members,
            "original_filename_candidates": candidates, "treefam_named_members": references}


def inspect(destination, tags, cache_dirs=()):
    if (not tags or len(set(tags)) != len(tags)
            or any(not isinstance(tag, str) or not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}", tag)
                   for tag in tags)):
        raise ValueError("Require unique safe tag names")
    destination.mkdir(parents=True, exist_ok=False)
    opener = build_opener(PublicRedirect())
    evidence = []

    def fetch(url, name, limit, token=None, expected=None, manifest=False):
        headers = {"User-Agent": "OrthoHMM-public-reference-audit", "Accept-Encoding": "identity"}
        if token:
            headers["Authorization"] = "Bearer " + token
        if manifest:
            headers["Accept"] = ACCEPT
        with opener.open(Request(url, headers=headers), timeout=30) as response:
            data = response.read(limit + 1)
            if len(data) > limit:
                raise ValueError("Response exceeds prospective byte limit")
            if expected is not None:
                if isinstance(expected, str):
                    if digest(data) != expected:
                        raise ValueError("Manifest digest mismatch")
                else:
                    verified(data, expected)
            path = destination / name
            path.write_bytes(data)
            evidence.append({"file": name, "url": url, "bytes": len(data),
                             "digest": digest(data), "http_status": response.status})
            return data

    tag_data = json.loads(fetch(TAGS_URL, "tags.json", 2_000_000))
    if tag_data["next"] is not None:
        raise ValueError("Tag listing requires pagination; this bounded investigation stops")
    auth_url = ("https://auth.docker.io/token?service=registry.docker.io"
                f"&scope=repository:{REPOSITORY}:pull")
    with opener.open(auth_url, timeout=30) as response:
        token = json.load(response)["token"]
    rows = []
    unique_layers = {}
    total, reused = 0, []
    for tag in tags:
        row = {"tag": tag, "layers": []}
        rows.append(row)
        try:
            raw = fetch(f"{REGISTRY}/manifests/{tag}", f"{tag}.manifest.json",
                        2_000_000, token=token, manifest=True)
            row["tag_manifest_digest"] = digest(raw)
            manifest = json.loads(raw)
            if "manifests" in manifest:
                platforms = [item for item in manifest["manifests"]
                             if item.get("platform", {}).get("os") == "linux"
                             and item.get("platform", {}).get("architecture") == "amd64"]
                if len(platforms) != 1:
                    raise ValueError("Expected exactly one linux/amd64 platform")
                platform = platforms[0]
                raw = fetch(f"{REGISTRY}/manifests/{platform['digest']}",
                            f"{tag}.amd64.manifest.json", 2_000_000,
                            token=token, expected=platform, manifest=True)
                manifest = json.loads(raw)
            row["image_manifest_digest"] = digest(raw)
            config = json.loads(fetch(f"{REGISTRY}/blobs/{manifest['config']['digest']}",
                                      f"{tag}.config.json", 2_000_000, token=token,
                                      expected=manifest["config"]))
            row["total_image_layers"] = len(manifest["layers"])
            selected = context_layers(manifest, config)
            if not selected:
                raise ValueError("No build-context layer recognized; no absence inference")
            for index, layer, created_by in selected:
                layer_digest = layer["digest"]
                if layer_digest not in unique_layers:
                    data, cache_ref = cached_layer(layer, cache_dirs)
                    filename = layer_digest.split(":", 1)[1] + ".tar.gz"
                    if data is None:
                        if total + layer["size"] > 96_000_000:
                            raise ValueError("Build-context download exceeds prospective byte budget")
                        data = fetch(f"{REGISTRY}/blobs/{layer_digest}", filename,
                                     32_000_000, token=token, expected=layer)
                        total += len(data)
                    else:
                        reused.append(cache_ref)
                        (destination / filename).write_bytes(data)
                    result = inventory(data)
                    unique_layers[layer_digest] = result
                    inv_name = layer_digest.split(":", 1)[1] + ".inventory.json"
                    (destination / inv_name).write_text(json.dumps(result, indent=2) + "\n")
                row["layers"].append({"index": index, **layer, "created_by": created_by,
                                      **{key: value for key, value in unique_layers[layer_digest].items()
                                         if key != "members"}})
            row["status"] = "selected_context_layers_inspected"
        except (URLError, OSError, ValueError, KeyError, tarfile.TarError) as error:
            row["status"] = "unresolved"
            row["error"] = str(error)
            if isinstance(error, HTTPError) and error.code == 429:
                print(tag, row["status"], flush=True)
                rows.extend(dict(tag=remaining, layers=[], status="not_attempted_rate_limit")
                            for remaining in tags[len(rows):])
                break
        print(tag, row["status"], flush=True)
    report = {"repository": REPOSITORY, "selected_tags": tags,
              "listed_tag_count": tag_data["count"], "tag_list_complete": True,
              "downloaded_unique_context_bytes": total, "reused_contexts": reused,
              "results": rows, "evidence": evidence,
              "containers_executed": False, "archive_members_extracted": False,
              "scope": "Selected build-context layers only; not all layers, tags, or registries"}
    (destination / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--tags", nargs="+", default=["2020.1", "2020.2", "2022.1"])
    parser.add_argument("--cache-dir", type=Path, action="append", default=[])
    args = parser.parse_args()
    inspect(args.destination, args.tags, args.cache_dir)
