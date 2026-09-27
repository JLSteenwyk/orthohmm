# Base Runtime Dependency Graph Check

The [graph receipt](integrated_base_dependency_graph_20260927.json) expands
the ten Conda owners in the observed search-runtime inventory through their
recorded `depends` fields. Conda's installed `MatchSpec` and `PackageRecord`
APIs check version/build constraints; dependency strings are not interpreted
with a custom parser. All 19 selected metadata-file identities were rechecked.
The receipt retains exact package filenames, URLs, archive hashes, declared
licenses, dependencies, constraints, parser identity and collection command.

## Findings

- Nineteen recorded packages are reachable from the ten roots.
- All 38 dependency edges with matching local records satisfy their specs.
- Three conditional constraints apply to selected packages and pass. Two
  constrain packages absent from the selected set and do not require them.
- The observed glibc 2.39 satisfies the declared `__glibc >=2.17` virtual
  dependency. This does not package glibc or establish cross-host compatibility.
- Python declares `pip`, but no local Conda package record owns a package
  with that name. The graph is therefore explicitly incomplete.

This is an ownership gap, not proof that pip is absent: a separate
`importlib.metadata.distribution("pip")` query reports version 26.0.1 in the
base Python's site-packages. That observation is not archive/payload
verification or a recommendation to reuse it. The previously validated
integrated workflow has a separately supplied installer and pins pip 26.2.1
inside its fresh environments. These identities must remain distinct.

An initial collection attempt omitted the required `build_number` in the
synthetic virtual-package record and stopped before writing a receipt.
Adding `build_number=0` allowed the installed Conda API to check the observed
glibc version. No environment or native execution was retried or changed.

## Bootstrap Implication

A fresh base-runtime recipe needs an explicit, validated installer bootstrap
and dependency policy; this graph must not be advertised as a complete Conda
lock. Preserve the 19 recorded native package identities, separately pin and
validate the installer/wheel overlay, and test the resulting new environment
before claiming reconstruction. Do not repair the running base environment
or let a solver silently choose new scientific dependencies.

Only five selected package archives have so far been newly acquired and
compared in this follow-up. The graph is not proof that all 19 downloads are
available, that an installation succeeds, or that every base-runtime file
matches its original archive. Remaining archive acquisition, restoration,
runtime closure, licensing/security review and publication gates stay open.
