Config.sh / CMake-YAML alignment: status and remaining work
*************************************************************

Background
==========

A collaborator building a CMake + YAML configuration pipeline
(``feature/config_to_yamls``) on top of the existing ``Config.sh`` /
Makefile system reported that their tooling failed against
``CosmoConfig.sh``, and raised two problems:

1. ``ADP_ACCRETION`` in :file:`src/io/parameters.c` sits outside any
   ``#ifdef BH_ACTIVE`` guard, alongside other inconsistent ``#ifdef``
   nesting elsewhere in the source tree.
2. There was no single source of truth for which flags depend on which,
   so new tooling had no reliable way to know, for example, that turning
   on ``BONDI_ACCRETION`` should also imply ``BH_ACTIVE``.

This page tracks what has been fixed, what remains, and what a
collaborator picking up the CMake/YAML integration should check before
merging that work into ``main``.

What's done
===========

.. list-table::
   :header-rows: 1

   * - Fix
     - PR
   * - Missing ``endif`` in the Makefile's ``FOF`` block (build-breaking regression)
     - `#44 <https://github.com/solas-sims/arepo_solas/pull/44>`_
   * - Duplicate ``reallocate_memory_maxpartbhs()`` in :file:`allocate.c`; wrong relative include paths in :file:`sfr_starbystar.c`
     - `#45 <https://github.com/solas-sims/arepo_solas/pull/45>`_
   * - Preprocessor self-guards added to 23 files whose bodies weren't wrapped in the ``#ifdef`` matching their own Makefile ``OBJS`` gate
     - `#46 <https://github.com/solas-sims/arepo_solas/pull/46>`_
   * - ``compute_mu()`` called unconditionally outside ``#ifdef USE_GRACKLE`` in two star-formation files
     - `#47 <https://github.com/solas-sims/arepo_solas/pull/47>`_
   * - ``config_flags.yaml`` added as the single source of truth for master-switch and derived flags, with a CI-enforced validator, and ``cmake/parse_yaml.py`` updated to compute the same derived flags the Makefile does
     - `#48 <https://github.com/solas-sims/arepo_solas/pull/48>`_ (open)

See :doc:`config-options` for the resulting three-tier flag model
(primary / master switch / derived) and how ``config_flags.yaml`` is
used.

Once PR #48 lands, the original ``ADP_ACCRETION``/``BH_ACTIVE`` example
is functionally resolved for both build paths: the Makefile already
derived ``BH_ACCRETION_ACTIVE``/``BH_ACTIVE`` correctly whenever
``ADP_ACCRETION`` was set, and ``cmake/parse_yaml.py`` now does the
same for a YAML config with ``adp_accretion: ON``. The code in
:file:`src/io/parameters.c` is no longer nested to match (see
`Remaining work`_, item 5), but it no longer produces an inconsistent
build.

Remaining work
==============

1. **Merge PR #48.** Until it lands, ``cmake/parse_yaml.py`` on
   ``main`` still silently drops every master-switch/derived flag, so
   any CMake build enabling e.g. ``bondi_accretion`` will compile
   without ``BH_ACTIVE`` ever being defined.

2. **No CMake/YAML equivalent of the Makefile's dependency validation.**
   The Makefile has a separate cascade of ``$(error ...)`` checks (e.g.
   *"BH_MERGER requires BH_ACTIVE"*, *"BH_ACTIVE requires BLACKHOLES"*,
   *"STAR_FEEDBACK_ACTIVE requires STARS"*) that reject invalid flag
   combinations at configure time. ``cmake/AREPOCMakeUtils.cmake`` has
   no equivalent logic at all, so an invalid YAML config (e.g.
   ``bh_merger: ON`` without any accretion or feedback model) will
   currently build silently instead of failing fast. This is a
   different tier from the master-switch/derived-flag work in PR #48
   -- it is about *requirements between primary flags*, not
   *derivation of new flags*.

   Recommended approach: extend ``config_flags.yaml`` with a
   ``requires:`` map mirroring the Makefile's error blocks, add the
   equivalent check to ``cmake/AREPOCMakeUtils.cmake``, and extend
   ``scripts/validate_config_flags.py`` to cross-check it against the
   Makefile's ``$(error ...)`` lines the same way it already
   cross-checks the derivation blocks. Not started.

3. **``feature/config_to_yamls`` needs to be rebased onto current
   ``main``** before it can be merged. It branched before the
   self-guard rollout, the ``compute_mu`` fix, and ``config_flags.yaml``
   landed, and it independently made partial ``#ifdef`` edits to some
   of the same files (:file:`bh_feedback.c`, :file:`bh_jet_density.c`,
   :file:`bh_seed.c`, :file:`bh_update.c`, :file:`init.c`,
   :file:`io/parameters.c`, :file:`criterion_derefinement.c`,
   :file:`sfr_AGORA.c`, :file:`sfr_JEANS.c`, :file:`sfr_eEOS.c`). These
   need to be reconciled against what's now on ``main`` so guards
   aren't duplicated or left inconsistent with ``config_flags.yaml``.

4. **``INDIVIDUAL_STAR_BY_STAR_FORMATION`` is still broken and
   untested**, independent of the pipeline work: :file:`sfr_starbystar.c`
   references ``SP[].Hsml`` without a ``STAR_FEEDBACK_ACTIVE`` guard and
   ``SP[].Metallicity`` without a ``METALS`` guard, and
   :file:`sfr_starbystar.c`/:file:`sfr_eEOS.c`  both define
   ``cooling_and_starformation()``, causing a duplicate-symbol link
   error when both are compiled in. This was deliberately parked rather
   than fixed opportunistically while chasing the ``compute_mu`` gap.
   It doesn't block the CMake/YAML integration unless an example
   config enables this flag.

5. **Cosmetic**: the ``ADP_ACCRETION``/``TORQUE_ACCRETION``/
   ``BH_ACCRETION_ACTIVE`` parameter blocks in
   :file:`src/io/parameters.c` (around line 680) could be nested inside
   a single ``#ifdef BH_ACTIVE`` for readability, matching the pattern
   used a few lines above for ``BhRadius``/``BhDesNgb``. No longer a
   functional bug after PR #48; low priority.

6. **Example YAML configs.** ``template_config_complete.yaml`` and the
   per-example ``Config.sh.yaml`` files added on
   ``feature/config_to_yamls`` live outside the repo root, so
   ``scripts/validate_config_flags.py``'s current glob
   (``template_config*.yaml`` at the repo root) won't catch a derived
   flag accidentally set in one of them. Either broaden the glob once
   that branch merges, or confirm none of those files set a derived
   flag directly.

Suggested order
================

#48 -> reconcile/rebase ``feature/config_to_yamls`` -> item 2
(dependency validation for CMake) -> item 6 (broaden validator glob) ->
merge the CMake/YAML pipeline into ``main``. Items 4 and 5 are
independent and can happen at any point.
