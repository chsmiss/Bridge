#!/usr/bin/env python3
"""Extract the first paired FASTQ URLs from an existing GitHub Actions workflow.

The Zymo protein pilot uses the same workflow file that produced the downloaded GFA
artifact, avoiding accidental comparison of a graph and a protein assembly from different
read sets.
"""

from __future__ import annotations

import argparse
import re
import shlex
from pathlib import Path
from typing import Optional, Sequence

URL_RE = re.compile(r"https?://[^\s'\"<>]+?(?:\.fastq\.gz|\.fq\.gz)")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workflow", type=Path, required=True)
    parser.add_argument("--github-env", type=Path)
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    text = args.workflow.read_text(encoding="utf-8")
    urls = []
    for match in URL_RE.finditer(text):
        url = match.group(0).rstrip("),]")
        if url not in urls:
            urls.append(url)
    if len(urls) < 2:
        raise SystemExit(f"could not find a FASTQ pair in {args.workflow}")

    # Prefer explicit R1/R2 or _1/_2 pairing when a workflow contains auxiliary URLs.
    first = next((url for url in urls if re.search(r"(?:R?1|_1)\.f(?:ast)?q\.gz$", url, re.I)), urls[0])
    second = next(
        (
            url
            for url in urls
            if url != first and re.search(r"(?:R?2|_2)\.f(?:ast)?q\.gz$", url, re.I)
        ),
        next(url for url in urls if url != first),
    )
    lines = [f"READ1_URL={first}", f"READ2_URL={second}"]
    if args.github_env:
        with args.github_env.open("a", encoding="utf-8") as handle:
            handle.write("\n".join(lines) + "\n")
    else:
        print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
