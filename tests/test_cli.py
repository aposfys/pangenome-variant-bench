"""The CLI surface. Added when the repo gained a console script it had never declared."""

from __future__ import annotations

import pytest

from panbench.cli import build_parser, main


def test_experiment_without_data_names_the_fetch_command(tmp_path):
    with pytest.raises(SystemExit) as excinfo:
        main(["--data-dir", str(tmp_path), "--results-dir", str(tmp_path), "experiment"])
    message = str(excinfo.value)
    assert "no sliced data" in message
    assert "panbench fetch" in message


def test_report_without_findings_says_what_to_run(tmp_path):
    with pytest.raises(SystemExit) as excinfo:
        main(["--results-dir", str(tmp_path), "report"])
    assert "panbench experiment" in str(excinfo.value)


def test_only_hg002_is_wired_up(tmp_path):
    with pytest.raises(SystemExit) as excinfo:
        main(["--data-dir", str(tmp_path), "fetch", "--sample", "HG003"])
    assert "HG002" in str(excinfo.value)


def test_subcommands_parse():
    parser = build_parser()
    assert parser.parse_args(["fetch"]).region == "chr20"
    assert parser.parse_args(["experiment"]).command == "experiment"
