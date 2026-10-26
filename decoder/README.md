# Decode base-pair probabilities in another project

`IPknot::decoder` is a C++17 library for layered RNA secondary-structure
decoding, including pseudoknots. It takes probabilities supplied by the caller;
it does not link probability engines, ViennaRNA, Python, spdlog, or cxxopts.
Refinement, automatic threshold search, and fixed-structure constraints remain
in the IPknot application.

All solver-specific C++ APIs are confined to `src/ip.cpp`. Structural decoding
and fixed-structure formulations use `IP` and the solver-independent `IPModel` instead of
solver headers or backend macros. `IPModelSolver` handles repeated optimization
with changing objectives, including backend-specific model reuse, behind that
same boundary.

## Use from CMake

Add the decoder from an IPknot checkout (keep the checkout's `src/` and `cmake/`
directories, which contain the shared solver implementations):

```cmake
add_subdirectory(path/to/ipknot/decoder)
target_link_libraries(my_project PRIVATE IPknot::decoder)
```

The decoder subdirectory defaults to solver-free DD. To enable ILP, set
`ENABLE_ILP=ON` and select the same solver options as the IPknot CLI, e.g.
`ENABLE_HIGHS=ON`. The solver is selected at build time; a solver-enabled library
supports both DD and ILP at runtime. Check `ipknot::Decoder::ilp_available()`
before selecting ILP. Requesting an unavailable backend throws
`std::runtime_error`, including for empty input.

You can also build and install just the library:

```sh
cmake -S decoder -B build-decoder -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_ILP=OFF -DIPKNOT_DECODER_BUILD_EXAMPLE=ON
cmake --build build-decoder -j
ctest --test-dir build-decoder --output-on-failure
cmake --install build-decoder --prefix /path/to/install
```

Then configure a consumer with `CMAKE_PREFIX_PATH=/path/to/install` and use:

```cmake
find_package(IPknotDecoder CONFIG REQUIRED)
target_link_libraries(my_project PRIVATE IPknot::decoder)
```

The top-level build also provides this target. Use `IPKNOT_BUILD_CLI=OFF` to
exclude the CLI and probability engines, and `BUILD_TESTING=OFF` to exclude tests.
The library has the same GPL-3.0-or-later license as IPknot; see [COPYING](../COPYING).

## Minimal C++ example

```cpp
#include <ipknot/decoder.h>

ipknot::DecoderOptions options;
options.thresholds = {0.2f, 0.1f};  // two pseudoknot levels
// Improved beam is the default, matching the CLI and recommended LPC settings.
ipknot::Decoder decoder(options);

// One entry per physical pair, with zero-based left < right.
std::vector<ipknot::PairProbability> probabilities{
    {0, 7, 0.9f}, {1, 6, 0.9f}, {3, 10, 0.9f}, {4, 9, 0.9f}};
auto result = decoder.decode(11, probabilities);

// result.bpseq[i]: zero-based partner, or -1 if unpaired
// result.levels[i]: zero-based level, or -1 if unpaired
// Both endpoints have the same level.
// result.objective: sum of the selected pair weights
// result.dd: optional DD diagnostics (iterations, bounds, pruning, stop reason)
```

The result above contains partners `{7,6,-1,10,9,-1,1,0,-1,4,3}` and objective
`1.5`. The sequence length is explicit so trailing unpaired bases are retained.
The decoder can be reused across calls. See
[examples/decode_probabilities.cpp](examples/decode_probabilities.cpp) for a
complete executable; its optional argument is `improved-beam` (the example's
default), `beam`, `nussinov`, `wavefront`, `full-dd`, or `ilp`.

## Linear and nonlinear decoding

| Mode | Settings | Scope |
| --- | --- | --- |
| Improved beam DD (default; recommended for LPC) | `options.backend = ipknot::Backend::DD; options.dd.improved_beam = true;` | Dominance pruning before beam selection; expected O(n+M) time and memory for fixed levels, beam/witness widths, iterations, and recovery budgets |
| Standard beam DD | `options.dd.improved_beam = false` | Finite beam selection without suffix dominance; same fixed-budget complexity |
| DD with exact level DP | `options.dd.nussinov_dp = true` | Unpruned Nussinov per level; worst-case cubic time and quadratic DP memory; crossing/witness limits still apply |
| DD with exact wavefront DP | `options.dd.wavefront_dp = true; options.dd.dp_threads = 4;` | Exact quadratic chart with within-DP parallelism; crossing/witness limits still apply |
| DD with the full crossing graph | Exact level DP, plus `options.dd.crossing_beam = 0; options.dd.witnesses = 0;` | Removes crossing/witness truncation too; no linear complexity guarantee |
| ILP | `options.backend = ipknot::Backend::ILP` | Full layered integer-programming formulation; requires a linked solver; no linear complexity guarantee |

For bounded DD, the defaults are DP beam 100, crossing beam 100, 16 witnesses,
and 50 iterations. Improved beam is the default in both the library and CLI;
`DDOptions` and the low-level `DDNussinov` constructor also default to improved
beam. Set `options.dd.improved_beam = false` to select standard beam explicitly.
Set `options.dd.threads = 2; options.dd.dp_threads = 1;` to reduce single-RNA latency
with OpenMP. Use one level thread for batches of independent RNAs.
Zero beam/crossing/witness widths mean unlimited widths.
The existing DD bound, recovery, exchange, and tracing settings are available
through `options.dd`. `dd.enabled` is unused here;
`options.backend` selects the backend.

Exact level DP applies suffix dominance and omits decomposable closed candidates
without beam pruning. Tied optima may trace back to a different structure than
older revisions. Worst-case cubic time and quadratic chart memory still apply.
Set `options.dd.threads = 2` to run independent level oracles concurrently;
the default is 1, and 0 uses the OpenMP runtime limit capped at the level count.
Serial and parallel executions aggregate scores in the same order. Recovery
remains sequential. CMake discovers OpenMP by default; disable it with
`-DENABLE_OPENMP=OFF`. Builds without OpenMP accept 0/1 and reject larger values.
Installed CMake packages propagate the OpenMP dependency to consumers.
The library compiles the same DD, sparsification, and parallel DP implementations
as the CLI. Tests build a separate consumer against the installed headers and
`IPknot::decoder` target, and verify improved beam, sparse exact DP, wavefront DP,
and serial/parallel consistency through the public API. They also verify that
omitting the beam mode matches explicit improved beam and that explicitly
selecting standard beam remains effective.

Set `options.dd.wavefront_dp = true; options.dd.dp_threads = 4;` for exact
single-DP wavefront parallelism. The beam is ignored. This backend keeps
closed-candidate sparsification and stacking states, and allocates roughly
`6*n*(n+1)` bytes for its score/choice chart plus candidate and stack buffers.
Its optimal score matches sparse exact Nussinov; tied structures can differ.
With finite beams, `dp_threads` instead parallelizes each column's components
and prefix splits with a deterministic merge before beam pruning, preserving
serial output. Small workloads automatically stay serial. `dp_threads=0`
uses the runtime limit with work-size caps for beams. Requests other than 1
take precedence over level threads; no nested teams or global runtime settings
are used. Sparse exact (`beam=0`, `wavefront_dp=false`) remains serial within a
level. Recovery preserves level order and can use the same inner workers.

Exact Nussinov solves each *level subproblem*. DD may still stop with a gap
between its feasible solution and upper bound. Its bound certifies the retained
witness graph; with truncated crossings/witnesses it is not a certificate for
the full ILP. ILP and DD can return different structures and resolve ties
differently. The library preserves the existing stacking rules: DD requires an
immediately adjacent stacked pair on the same level; ILP requires neighboring
paired left endpoints and neighboring paired right endpoints on that level.
Set `options.no_lonely_pairs = false` to disable these rules.

## Probability formats and scoring

Probabilities must be finite and in `[0,1]`; endpoints must lie inside the
sequence and each physical pair must occur only once. Entries need not be sorted.
Duplicate pairs, invalid indices, or invalid options throw `std::invalid_argument`.
Omitted pairs have no candidate variable. The decoder does not impose nucleotide
types or a minimum hairpin length: supply only the pairs allowed by your model.

For level `l`, a pair is retained when `p > thresholds[l]` and receives weight
`(p - thresholds[l]) * weights[l]`. Thresholds must be in `[0,1]`, weights must
be positive and finite, and their lengths must agree. Empty weights default to
`1 / number_of_levels`; the default thresholds are `{0.25f, 0.125f}`. Use one
threshold for pseudoknot-free decoding. Each base has at most one partner;
pairs within a level do not cross; every pair on a higher level must cross a
selected pair on each lower level.

The native IPknot sparse format can be passed without conversion:

```cpp
ipknot::SparsePosterior posterior(12); // n+1 rows; row 0 must be empty
posterior[1] = {{8, 0.9f}};           // one-based row and partner indices
posterior[2] = {{7, 0.9f}};
posterior[4] = {{11, 0.9f}};
posterior[5] = {{10, 0.9f}};
auto result = decoder.decode(posterior);
```

An upper-triangular or symmetric sparse matrix is accepted. Only upper entries
contribute candidates; lower entries are validated and then ignored. Duplicate
upper entries are rejected. `decode(0, {})` or a sparse matrix containing only
the empty row 0 returns an empty structure. An entirely empty sparse matrix is
invalid because it lacks the required row 0.
