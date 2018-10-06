# BabrahamLinkON 0.2

Most changes pertain to the short pipeline. Moved main scripts into bin folder.

## Precleaning

* The sequence beyond anchor is no longer removed from output reads.
* J end UMI moved from deduplication, so that all UMI extraction is now preformed at precleaning stage.
* Fixed comments in R script that caused J germline plots not to be produced.
* Check if germline plot is produced.
* Fixed an incorrect mispriming correction offset for the kappa locus.

#### Mispriming error estimate
* Extended error estimation for human sequences.

#### For short reads:
* json file written with identity of reads assembled and unassembled to allow merging after deduplication.

## Deduplication
* Option to output sequences with ambiguous N nucleotides (normally these are filled in with basepairs from sequences best matching the consensus).
* Fixed bug in UMI report tables (most likely cause by pandas update).
* Simplified options, made MSA available to all pipelines.
* Moved general functions into deduplicattion_general.py and split deduplicate_bundle_parallel for clarity.
* Updated UMI correction to work with latest version.

## Annotation and clone assembly
* Added json option so output reads can be marked assembled or unassembled.
* DJ reads are also filtered by the V end IgBlast calls, instead on just J end.

## Setup

* Moved package data into MANIFEST.in to enable `pip --editable` mode to work with pkg_resources.resource_filename.


# BabrahamLinkON 0.1

* Initial release
