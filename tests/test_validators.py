"""Tests for validators."""

from bioclaw.utils.validators import validate_file_path, validate_smiles


def test_validate_smiles_valid():
    assert validate_smiles("CCO") is True
    assert validate_smiles("CC(=O)Oc1ccccc1C(=O)O") is True
    assert validate_smiles("c1ccccc1") is True


def test_validate_smiles_invalid():
    assert validate_smiles("") is False
    assert validate_smiles(None) is False
    assert validate_smiles("not a smiles!") is False


def test_validate_file_path():
    assert validate_file_path("test.csv") is True
    assert validate_file_path("") is False
    assert validate_file_path("../etc/passwd") is False


def test_validate_file_path_extensions():
    assert validate_file_path("test.csv", [".csv", ".tsv"]) is True
    assert validate_file_path("test.pdf", [".csv", ".tsv"]) is False
