Plan for assembly is to try a few methods and see what comes out well

The ONT pipeline is established and designed to work with nanopore

Wtdbg2 designed to run faster and more efficiently than canu, with the space to polish later

Falcon designed to be aware of heterozygosity in genomes - might be useful as jellyfish indicates a high level of het (2 distinct k-mer peaks)

Devs reccomend between 10 and 20 tb of physical space

Can alter number of cores - how many should be used here?


suggested paramaters to pass to canu

genome size = 1.5g

corMaxEvidenceErate=0.15 to correct for repetitive sequence (suggested for plants by devs)

rawErrorRate= ?? is the maximum expected difference in an alignment of two _uncorrected_ reads. It is a meta-parameter that sets other parameters. Default for now?

minReadLength and minOverlapLength. The defaults are to discard reads shorter than 1000bp and to not look for overlaps shorter than 500bp. Increasing minReadLength can improve run time, and increasing minOverlapLength can improve assembly quality by removing false overlaps. However, increasing either too much will quickly degrade assemblies by either omitting valuable reads or missing true overlaps. Default for now

corOutCoverage= default is to expect 40X - devs suggest incressing this to allow more contigs to obtain two haplotypes in polypolid genomes - could do this here for heterozygosity






