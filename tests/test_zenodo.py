"""Test for zenodo.json schema conformance.

Based on [this PR](https://github.com/ProjectDrawdown/solutions/pull/165), but updated
with the [latest schema from
Zenodo](https://github.com/zenodo/zenodo/blob/master/zenodo/modules/deposit/jsonschemas/deposits/records/legacyrecord.json).
This is also referred to in the [GitHub section](https://developers.zenodo.org/#github).
"""

import json
import pathlib

import jsonschema

thisdir = pathlib.Path(__file__).parents[0]
parentdir = pathlib.Path(__file__).parents[1]


def test_zenodo_json_schema():
    """Test the .zenodo.json file against the Zenodo schema."""
    with open(str(thisdir.joinpath("zenodo_v1.schema.json")), encoding="locale") as s:
        schema = json.loads(s.read())
    with open(str(parentdir.joinpath(".zenodo.json")), encoding="locale") as z:
        instance = json.loads(z.read())
    jsonschema.validate(instance=instance, schema=schema)
