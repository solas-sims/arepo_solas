Arepo_SOLAS version
====================

Modified from the AREPO public version for the SOLAS project.

Documentation is in `documentation/source/` (Sphinx):

- [SOLAS parameter reference](documentation/source/solas-parameters.md)
- [SIDM](documentation/source/sidm.md)
- [Examples](documentation/source/examples.md)
- [Technical reference](documentation/source/technical-reference.md)

The [wiki](https://github.com/solas-sims/arepo_solas/wiki) has an older high-level description ([Code](https://github.com/solas-sims/arepo_solas/wiki/Code), [Parameters](https://github.com/solas-sims/arepo_solas/wiki/Parameters)); where it disagrees with the in-repo documentation, the in-repo documentation is current.

How to contribute
=================

1. Fork from or create a `new branch` at the [official repository](https://github.com/solas-sims/arepo_solas).
2. Make your changes in a new, appropriately named branch (e.g. `starformation` or `bh_feedback`).
    * All git commit messages should loosely follow the [standard format outlined here](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html).
    * Ensure your new code is commented enough so that others can understand what it's doing. If possible, update the [wiki](https://github.com/solas-sims/arepo_solas/wiki) as well.
3. Issue a [pull request](https://help.github.com/en/articles/about-pull-requests) to the official repository.

### IMPORTANT

1. Direct pushes to the `main` have been disabled; one has to merge through a PR. Code review is not enforced but having a PR will help you check whether you had inadvertently changed something.
2. Delete stale branches after merging.
