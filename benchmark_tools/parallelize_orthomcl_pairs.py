#!/usr/bin/env python3
"""Add opt-in process parallelism to OrthoMCL 1.4 species-pair inference."""

from __future__ import annotations

import argparse
from pathlib import Path


LOOP_START = "for(my $i=0;$i<scalar(@taxa)-1;$i++) {"
INNER_START = "\tfor(my $j=$i+1;$j<scalar(@taxa);$j++) {"
LOOP_END_MARKER = "\n\n%blastquery=();"
PARALLEL_MARKER = "ORTHOMCL_PAIR_WORKERS"


def parallelize_source(source: str) -> str:
    """Replace the serial inter-taxon loop with an equivalent forked loop."""
    if PARALLEL_MARKER in source:
        raise ValueError("OrthoMCL source is already pair-parallel")
    start = source.find(LOOP_START)
    end = source.find(LOOP_END_MARKER, start)
    if start < 0 or end < 0:
        raise ValueError("Could not locate the OrthoMCL inter-taxon loop")

    original = source[start:end].rstrip()
    lines = original.splitlines()
    if lines[:2] != [LOOP_START, INNER_START] or lines[-2:] != ["\t}", "}"]:
        raise ValueError("Unexpected OrthoMCL inter-taxon loop structure")

    body = "\n".join(lines[2:-2])
    body = body.replace("$taxa[$i]", "$ta").replace("$taxa[$j]", "$tb")
    replacement = f"""my $process_intertaxon_pair = sub {{
\tmy ($ta, $tb) = @_;
{body}
\treturn $connect{{$ta.' '.$tb}};
}};

my @intertaxon_pairs;
for(my $i=0;$i<scalar(@taxa)-1;$i++) {{
\tfor(my $j=$i+1;$j<scalar(@taxa);$j++) {{
\t\tpush @intertaxon_pairs, [$taxa[$i], $taxa[$j]];
\t}}
}}

my $pair_workers = $ENV{{'ORTHOMCL_PAIR_WORKERS'}} || 1;
if ($pair_workers <= 1) {{
\tforeach my $pair (@intertaxon_pairs) {{
\t\t$process_intertaxon_pair->(@$pair);
\t}}
}} else {{
\trequire IO::Handle;
\trequire POSIX;
\trequire Storable;
\twrite_log("\\nRunning inter-taxon inference with $pair_workers workers\\n");
\torthomcl_module::LOG->flush();
\torthomcl_module::BBH->flush();

\tmy %children;
\tmy %result_files;
\tmy $next_pair = 0;
\twhile ($next_pair < scalar(@intertaxon_pairs) || scalar(keys %children)) {{
\t\twhile ($next_pair < scalar(@intertaxon_pairs) &&
\t\t       scalar(keys %children) < $pair_workers) {{
\t\t\tmy $pair_id = $next_pair;
\t\t\tmy ($ta, $tb) = @{{$intertaxon_pairs[$pair_id]}};
\t\t\tmy $pair_key = $ta.' '.$tb;
\t\t\tmy $result_file = dirname($bpo_file)."/pair_$pair_id.storable";
\t\t\t$result_files{{$pair_key}} = $result_file;
\t\t\tmy $pid = fork();
\t\t\tdieWithUnexpectedError("fork failed: $!") unless defined $pid;
\t\t\tif ($pid == 0) {{
\t\t\t\topen_bpofile($bpo_file);
\t\t\t\t%ortho = ();
\t\t\t\t$process_intertaxon_pair->($ta, $tb);
\t\t\t\tStorable::nstore(
\t\t\t\t\t[$connect{{$pair_key}}, [keys %ortho]],
\t\t\t\t\t$result_file
\t\t\t\t);
\t\t\t\torthomcl_module::LOG->flush();
\t\t\t\torthomcl_module::BBH->flush();
\t\t\t\tPOSIX::_exit(0);
\t\t\t}}
\t\t\t$children{{$pid}} = $pair_key;
\t\t\t$next_pair++;
\t\t}}

\t\tmy $finished = wait();
\t\tdieWithUnexpectedError("wait failed: $!") if $finished < 0;
\t\tmy $pair_key = delete $children{{$finished}};
\t\tdieWithUnexpectedError("worker failed for $pair_key") if $? != 0;
\t}}

\tforeach my $pair (@intertaxon_pairs) {{
\t\tmy ($ta, $tb) = @$pair;
\t\tmy $pair_key = $ta.' '.$tb;
\t\tmy ($pair_connect, $pair_ortho) =
\t\t\t@{{Storable::retrieve($result_files{{$pair_key}})}};
\t\t$connect{{$pair_key}} = $pair_connect;
\t\tforeach my $gene (@$pair_ortho) {{ $ortho{{$gene}} = 1; }}
\t\tunlink $result_files{{$pair_key}};
\t}}
}}
"""
    return source[:start] + replacement.rstrip() + source[end:]


def parallelize_file(path: Path) -> None:
    path.write_text(parallelize_source(path.read_text()))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("orthomcl_script", type=Path)
    args = parser.parse_args()
    parallelize_file(args.orthomcl_script)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
