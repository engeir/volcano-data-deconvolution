"""Test for JSON schema conformance.

Based on [this PR](https://github.com/ProjectDrawdown/solutions/pull/165), but updated
with the latest schema from Zenodo and release-please. This is also referred to in the
[GitHub section](https://developers.zenodo.org/#github).
"""

import json
import pathlib
import urllib.request

import jsonschema

parentdir = pathlib.Path(__file__).parents[1]


def test_zenodo_json_schema() -> None:
    """Test the .zenodo.json file against the Zenodo schema."""
    with urllib.request.urlopen(
        "https://github.com/zenodo/zenodo/raw/master/zenodo/modules/deposit/jsonschemas/deposits/records/legacyrecord.json"
    ) as s:
        schema = json.loads(s.read())
    with open(
        str(parentdir.joinpath(".zenodo.json")),
        encoding="locale",
    ) as z:
        instance = json.loads(z.read())
    jsonschema.validate(instance=instance, schema=schema)


def test_release_please_json_schema() -> None:
    """Test the release-please-config.json file against the release-please schema."""
    with urllib.request.urlopen(
        "https://github.com/googleapis/release-please/raw/main/schemas/config.json"
    ) as s:
        schema = json.loads(s.read())
    with open(
        str(parentdir.joinpath("./.github/release-please/release-please-config.json")),
        encoding="locale",
    ) as z:
        instance = json.loads(z.read())
    jsonschema.validate(instance=instance, schema=schema)
