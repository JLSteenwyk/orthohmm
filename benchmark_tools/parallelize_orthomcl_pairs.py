#!/usr/bin/env python3
"""Add opt-in process parallelism to OrthoMCL 1.4 species-pair inference."""

from __future__ import annotations

import argparse
from pathlib import Path


LOOP_START = "for(my $i=0;$i<scalar(@taxa)-1;$i++) {"
INNER_START = "\tfor(my $j=$i+1;$j<scalar(@taxa);$j++) {"
LOOP_END_MARKER = "\n\n%blastquery=();"
PARALLEL_MARKER = "ORTHOMCL_PAIR_WORKERS"
FORWARD_LOOKUP = """\t\t\t\t\tif (blastqueryab($nodes1[$k],$nodes2[$l])) {
\t\t\t\t\t\tmy ($s,$pm,$pe,$pi)=(blastqueryab($nodes1[$k],$nodes2[$l]))[0,3,4,5];"""
FORWARD_LOOKUP_ONCE = """\t\t\t\t\tmy @forward_hit=$cached_blastqueryab->($nodes1[$k],$nodes2[$l]);
\t\t\t\t\tif (@forward_hit) {
\t\t\t\t\t\tmy ($s,$pm,$pe,$pi)=@forward_hit[0,3,4,5];"""
REVERSE_LOOKUP = """\t\t\t\t\tif (blastqueryab($nodes2[$l],$nodes1[$k])) {
\t\t\t\t\t\tmy ($s,$pm,$pe,$pi)=(blastqueryab($nodes2[$l],$nodes1[$k]))[0,3,4,5];"""
REVERSE_LOOKUP_ONCE = """\t\t\t\t\tmy @reverse_hit=$cached_blastqueryab->($nodes2[$l],$nodes1[$k]);
\t\t\t\t\tif (@reverse_hit) {
\t\t\t\t\t\tmy ($s,$pm,$pe,$pi)=@reverse_hit[0,3,4,5];"""


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
    if body.count(FORWARD_LOOKUP) != 1 or body.count(REVERSE_LOOKUP) != 1:
        raise ValueError("Could not locate duplicate OrthoMCL co-ortholog lookups")
    body = body.replace(FORWARD_LOOKUP, FORWARD_LOOKUP_ONCE)
    body = body.replace(REVERSE_LOOKUP, REVERSE_LOOKUP_ONCE)
    replacement = f"""my $process_intertaxon_pair = sub {{
\tmy ($ta, $tb) = @_;
\tmy %query_hit_cache;
\tmy $cached_blastqueryab = sub {{
\t\tmy ($query, $subject) = @_;
\t\tunless (exists $query_hit_cache{{$query}}) {{
\t\t\tmy %hits;
\t\t\tif (defined $blastquery{{$query}}) {{
\t\t\t\tmy ($start, $end) = split(";", $blastquery{{$query}});
\t\t\t\tforeach my $line_id ($start..$end) {{
\t\t\t\t\tmy @hit = (getline_from_bpofile($line_id))[0,1,3,5,6,7];
\t\t\t\t\t$hits{{$hit[2]}} = \\@hit unless exists $hits{{$hit[2]}};
\t\t\t\t}}
\t\t\t}}
\t\t\t$query_hit_cache{{$query}} = \\%hits;
\t\t}}
\t\treturn @{{$query_hit_cache{{$query}}{{$subject}}}}
\t\t\tif exists $query_hit_cache{{$query}}{{$subject}};
\t\treturn 0;
\t}};
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
\tmy @pending_pair_ids;
\tfor (my $pair_id=0; $pair_id<scalar(@intertaxon_pairs); $pair_id++) {{
\t\tmy ($ta, $tb) = @{{$intertaxon_pairs[$pair_id]}};
\t\tmy $pair_key = $ta.' '.$tb;
\t\tmy $result_file = dirname($bpo_file)."/pair_$pair_id.storable";
\t\t$result_files{{$pair_key}} = $result_file;
\t\tif (-s $result_file && eval {{ Storable::retrieve($result_file); 1 }}) {{
\t\t\twrite_log("Reusing completed inter-taxon result for $ta and $tb\\n");
\t\t}} else {{
\t\t\tunlink $result_file if -e $result_file;
\t\t\tpush @pending_pair_ids, $pair_id;
\t\t}}
\t}}

\tmy $next_pending = 0;
\twhile ($next_pending < scalar(@pending_pair_ids) || scalar(keys %children)) {{
\t\twhile ($next_pending < scalar(@pending_pair_ids) &&
\t\t       scalar(keys %children) < $pair_workers) {{
\t\t\tmy $pair_id = $pending_pair_ids[$next_pending];
\t\t\tmy ($ta, $tb) = @{{$intertaxon_pairs[$pair_id]}};
\t\t\tmy $pair_key = $ta.' '.$tb;
\t\t\tmy $result_file = $result_files{{$pair_key}};
\t\t\tmy $pid = fork();
\t\t\tdieWithUnexpectedError("fork failed: $!") unless defined $pid;
\t\t\tif ($pid == 0) {{
\t\t\t\topen_bpofile($bpo_file);
\t\t\t\t%ortho = ();
\t\t\t\t$process_intertaxon_pair->($ta, $tb);
\t\t\t\tmy $temporary_result = $result_file.".".$$;
\t\t\t\tStorable::nstore(
\t\t\t\t\t[$connect{{$pair_key}}, [keys %ortho]],
\t\t\t\t\t$temporary_result
\t\t\t\t);
\t\t\t\trename $temporary_result, $result_file or
\t\t\t\t\tdieWithUnexpectedError("cannot publish $result_file: $!");
\t\t\t\torthomcl_module::LOG->flush();
\t\t\t\torthomcl_module::BBH->flush();
\t\t\t\tPOSIX::_exit(0);
\t\t\t}}
\t\t\t$children{{$pid}} = $pair_key;
\t\t\t$next_pending++;
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
