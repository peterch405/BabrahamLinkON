#BabrahamLinkON 0.2

Most changes pertain to the short pipeline.

##Precleaning

* The sequence beyond anchor is no longer removed from output reads.
* J end UMI moved from deduplication, so that all UMI extraction is now preformed at precleaning stage.
* Fixed comments in R script that caused J germline plots not to be produced.
* Check if germline plot is produced.
* Fixed an incorrect mispriming correction offset for the kappa locus.

####Mispriming error estimate
* Extended error estimation for human sequences.

####For short reads:
* json file written with identity of reads assembled and unassembled to allow merging after deduplication.

##Deduplication
* Option to output sequences with ambiguous N nucleotides (normally these are filled in with basepairs from sequences best matching the consensus).
* Fixed bug in UMI report tables (most likely cause by pandas update).
* Simplified options, made MSA available to all pipelines.

##Annotation and clone assembly
* Added json option so output reads can be marked assembled or unassembled.
* DJ reads are also filtered by the V end IgBlast calls, instead on just J end.


#BabrahamLinkON 0.1

* Initial release
