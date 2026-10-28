Code development
************************


Scientific software development is a crucial aspect of computational
(astro-)physics.  With the progess made over the past decades,
numerical methods as well as implementations have evolved and
increased in complexity. The scope of simulation codes has increased
accordingly, such that the question on how to organize development
becomes important for further scientific progress.

This is the SOLAS fork of Arepo, which does extend the public code (star formation and feedback,
radiation, black holes, SIDM and more). New physics is welcome here; follow the contribution workflow in the
repository `README.md` (a branch per change, merged through a pull request, stale branches deleted after merging).
The text below about the public release applies to the upstream code, not to this fork.

The public version of Arepo is intended as a basis, providing an efficient code structure and state of the art
numerical techniques for gravity and hydrodynamics, and does not accept extensions in functionality: only
bugfixes, documentation and examples are included in revisions of the public release.

Issue reporting
===============

A code of the scope of Arepo will inevitably have a number of issues
and so far undiscovered bugs.  Detection of all issues is very
complicated and time-consuming work. This means in practice that we
rely on users to actually report all bugs and issues they come across,
helping us to improve the quality of the code. We therefore encourage
users to report all issues they have, including things that remain
unclear after reading the documentation. This way we hope to
constantly improve the qualtiy of the code in all its aspects.
 
Issues with the **public** Arepo release go to Rainer Weinberger (rweinberger@aip.de). Issues specific to this fork
should be raised on the [arepo_solas issue tracker](https://github.com/solas-sims/arepo_solas/issues).


Code extensions 
===============

Extensions should be hosted in separate repositories or branches of
the main repository.  We highly welcome such extensions and encourage
developers to make them publicly available (under GNU GPL v3). While
this can be done independently of the authors, we would encourage
developers to inform the authors once there is a publicly available
module they are willing to share, in order to have a list of available
extensions on the code homepage.

Some guidelines for code extensions and new modules:

  * Strive for modularity
 
  * Minimize the number of changes in existing source code files to a
    few function calls and structured variables within existing
    elements.
	
  * The source code of a new modue should largely be in separate
    (i.e. new) files.
  
  * Document the module
 
  * All parameter and config-options should be clearly explained in
    this documentation.  Feel free to add an addional page to this
    documentation explaining how your module works.
	
  * Document what each function does in the source code.
	
  * The Template-Config.sh file should have a short explanation for
    each added flag.
  
  * Verification and examples
 
  * If possible, create one or more addional examples illustrating and
    testing the module.


