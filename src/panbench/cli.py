"""Command line entry point: ``python -m panbench.cli`` or ``panbench``.

The Nextflow workflow is the caller comparison and needs Docker. Everything reachable from
here is the evaluation layer, which does not.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from panbench import __version__


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="panbench",
        description="Where does a pangenome reference actually help?",
    )
    parser.add_argument("--version", action="version", version=f"panbench {__version__}")
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    sub = parser.add_subparsers(dest="command", required=True)

    fetch = sub.add_parser("fetch", help="GIAB truth, confident regions and stratifications")
    fetch.add_argument("--sample", default="HG002", help="only HG002 is wired up")
    fetch.add_argument("--region", default="chr20")

    experiment = sub.add_parser(
        "experiment", help="the real stratification, and the evaluation-procedure simulation"
    )
    experiment.add_argument("--region", default="chr20")
    experiment.add_argument("--seed", type=int, default=0)

    sub.add_parser("report", help="print the summary from an existing findings.json")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return _dispatch(args)
    except RuntimeError as exc:
        raise SystemExit(str(exc)) from exc
    except OSError as exc:
        raise SystemExit(f"could not reach the GIAB FTP: {exc}") from exc


def _dispatch(args: argparse.Namespace) -> int:
    if args.command == "fetch":
        from panbench.fetch import fetch

        if args.sample != "HG002":
            raise SystemExit(f"only HG002 is wired up, not {args.sample!r}")
        result = fetch(args.data_dir, region=args.region)
        print(
            f"\n{result.n_truth_variants:,} variants, "
            f"{result.confident_bases:,} confident bases"
        )
        return 0

    if args.command == "experiment":
        from panbench.experiment import run

        sliced = args.data_dir / args.region
        truth = sliced / f"truth.{args.region}.vcf"
        if not truth.exists():
            # A missing download used to surface as FileNotFoundError from deep inside
            # the reader; say which command produces it instead.
            raise SystemExit(
                f"no sliced data at {sliced}.\n"
                f"Fetch it first:  panbench fetch --region {args.region}"
            )
        run(args.data_dir, args.results_dir, region=args.region, seed=args.seed)
        print(f"\nwrote {args.results_dir / 'findings.json'}")
        return 0

    if args.command == "report":
        import json

        path = args.results_dir / "findings.json"
        if not path.exists():
            raise SystemExit(f"no findings at {path}. Run 'panbench experiment' first.")
        findings = json.loads(path.read_text())
        truth = findings["truth"]
        print(
            f"{findings['region']}: {truth['n_variants']:,} truth variants, "
            f"{truth['confident_fraction']:.1%} of the chromosome confident"
        )
        for stratum, count in truth["by_stratum"].items():
            print(f"  {stratum:<24} {count:>7,}")
        print(f"\n{'arm':<18} {'F1 restricted':>14} {'F1 unrestricted':>16}")
        for name, arm in findings["simulation"]["arms"].items():
            print(
                f"{name:<18} {arm['restricted']['f1']:>14.4f} "
                f"{arm['unrestricted']['f1']:>16.4f}"
            )
        print(f"\nnot run: {findings['caller_comparison_not_run']}")
        return 0

    raise SystemExit(f"unknown command {args.command!r}")


if __name__ == "__main__":
    raise SystemExit(main())
