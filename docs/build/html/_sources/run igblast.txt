Run IgBlast with ChangeO
========================

For more information on how to set up IgBlast and run ChangeO please refer to the ChangeO documentation:
https://changeo.readthedocs.io/en/version-0.3.3---makedb-expansion/examples/igblast.html

To run igblast on low umi sequences::

    cat dedup.fasta low_umi.fasta > all.fasta


Select how many cores to use with IgBlast::

    cores=8
    pyfasta split -n $cores all.fasta


    for f in all.*.fasta; do
	echo $f
	igblastn \
	    -germline_db_V IgBlast_database/Mus_musculus_IGHV \
	    -germline_db_D IgBlast_database/Mus_musculus_IGHD \
	    -germline_db_J IgBlast_database/Mus_musculus_IGHJ \
	    -auxiliary_data igblast/optional_file/mouse_gl.aux \
	    -domain_system imgt -ig_seqtype Ig -organism mouse \
	    -outfmt '7 std qseq sseq btop' \
	    -query $f \
	    -out ${f}.fmt7 \
	    -num_threads 1 &
	
    done


    FAIL=0
    for job in `jobs -p`
    do
    echo $job
        wait $job || let "FAIL+=1"
    done

    echo $FAIL

    if [ "$FAIL" == "0" ];
    then
    echo "All finished succesfully!"
    else
    echo "FAIL! ($FAIL)"
    fi

Run ChangeO MakeDb.py to pasrse IgBlast results::

    for o in all.*.fmt7; do
	MakeDb.py igblast -i $o \
	-s ${o%.*} \
	-r IgBlast_database/Mus_musculus_IGH[JDV].fasta \
	--regions --scores
	
	#Delete all intermediate files
	rm $o
	rm ${o%.*}
    done

To assemble the clones into clonotypes::

    assemble_clones.py --tab_file 'all.*.fasta_db-pass.tab' --plot

