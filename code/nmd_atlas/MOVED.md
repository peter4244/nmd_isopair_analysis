# The atlas moved to its own repository

**`github.com/peter4244/nmd_lung_atlas_site`** — 2026-08-12.

The NMD Lung Atlas static site and its data exporter used to live here at `code/nmd_atlas/`.
They now have their own repository, which is what deploys
[nmd-lungcells.castaldilab.org](https://nmd-lungcells.castaldilab.org/).

**Nothing was lost.** The new repository was created with `git filter-repo`, so it carries all 15
commits that touched this directory, with their authors and dates. The SMG1i column work
(`atlas-smg1i-cpm-column`, `7666dfa`) is in it too, under its original author.

**Why it moved.** The same code existed in *two* repositories — here and in the analysis repo —
byte-identical on both mains, which is what a copy looks like immediately before it drifts. A
change would have landed in whichever one someone opened first. One copy now, and it is the one
that deploys.

**A file kept only as a pointer is still a copy of nothing.** This directory holds this note and
nothing else on purpose: a path that has moved should say where, not 404.
