import os

import pytest
from checkv import contamination, completeness, repeats, quality_summary


def test_database_exists(database):
    if not os.path.exists(database):
        assert 0


def test_contamination_module(database, threads):
    args = {
        "input": "test/test_sequences.fna",
        "output": "test/output_files",
        "db": database,
        "threads": threads,
        "restart": False,
        "quiet": False,
        "exclude": None,
    }
    contamination.main(args)
    with open("test/output_files/contamination.tsv") as test_output:
        with open("test/ground_truth/contamination.tsv") as ground_truth_output:
            for l1, l2 in zip(test_output, ground_truth_output):
                assert l1 == l2
    with open("test/output_files/viruses.fna") as test_output:
        with open("test/ground_truth/viruses.fna") as ground_truth_output:
            for l1, l2 in zip(test_output, ground_truth_output):
                assert l1 == l2
    with open("test/output_files/proviruses.fna") as test_output:
        with open("test/ground_truth/proviruses.fna") as ground_truth_output:
            for l1, l2 in zip(test_output, ground_truth_output):
                assert l1 == l2


def test_completeness_module(database, threads):
    args = {
        "input": "test/test_sequences.fna",
        "output": "test/output_files",
        "db": database,
        "threads": threads,
        "restart": False,
        "percent_of_top_hit": 50.0,
        "max_aai": None,
        "exclude_identical": False,
        "quiet": False,
        "exclude_list": None,
    }
    completeness.main(args)
    with open("test/output_files/completeness.tsv") as test_output:
        with open("test/ground_truth/completeness.tsv") as ground_truth_output:
            for l1, l2 in zip(test_output, ground_truth_output):
                assert l1 == l2


def test_repeats_module():
    args = {
        "input": "test/test_sequences.fna",
        "output": "test/output_files",
        "min_tr_len": 20,
        "max_tr_count": 5,
        "max_tr_dust": 20.0,
        "quiet": False,
    }
    repeats.main(args)
    with open("test/output_files/repeats.tsv") as test_output:
        with open("test/ground_truth/repeats.tsv") as ground_truth_output:
            for l1, l2 in zip(test_output, ground_truth_output):
                assert l1 == l2


def test_quality_summary_module():
    args = {
        "input": "test/test_sequences.fna",
        "output": "test/output_files",
        "quiet": False,
    }
    quality_summary.main(args)
    with open("test/output_files/quality_summary.tsv") as test_output:
        with open("test/ground_truth/quality_summary.tsv") as ground_truth_output:
            for l1, l2 in zip(test_output, ground_truth_output):
                assert l1 == l2
