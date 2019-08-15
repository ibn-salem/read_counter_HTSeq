""" 
Counts paired-end fragments from SAM/BAM files for each region (e.g. genes, transcripts or exons) in a .bed file.
"""
epilog="""Jonas Ibn-Salem <ibnsalem@molgen.mpg.de> 29.09.12"""

import argparse     # to parse commandline arguemnts
import HTSeq        # to parse regions, SAM/BAM file and count reads 
import collections  # for dict with default zero counts

def commandline():
    """ returns the args from commandline. """
    command_parser = argparse.ArgumentParser(description=__doc__, epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
    command_parser.add_argument('-i','--input_files', type=str, nargs='*', required=True, help='input SAM/BAM file(s).')
    command_parser.add_argument('-b','--bed_file', type=str, required=True, help='Bed file.')
    command_parser.add_argument('-l','--labels', type=str, nargs='*', default=[], help='labels for each sam in same number and order as input_files')
    command_parser.add_argument('-c','--use_chrom_name', action="store_true",  help='In case off mapping against transcriptome, use the chrom-name as feature for counting.'),
    command_parser.add_argument('-ft','--file_type', type=str, default="sam", choices=["sam", "bam"],  help='Read-alignment file format.'),
    command_parser.add_argument('-r','--order', type=str, default="name", choices=["name", "pos"],  help="Sorting order of <alignment_file> (default: name). Paired-end sequencing data must be sorted either by position or by read name, and the sorting order must be specified. Ignored for single-end data."),
    command_parser.add_argument('-o','--output_file', type=str, required=True, help='output file in bed format.')
    args =  command_parser.parse_args()
    return args

def getRegNamesFromBed(bed_file):
    """ return only the name column as list """
    return [line.strip().split("\t")[0] for line in open(bed_file)]

def getGenomicarrayOfSetsAndNames(bed_file):
    """ Returns a GenomicArrayOfSets of all regions and a list of region names """
    
    # build parser for rgions
    regionParser = HTSeq.BED_Reader( bed_file )    
    
    # build GenomicArrayOfSets for all regions
    regions = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
    
    # add all region to the GenomicArrayOfSets
    for feature in regionParser:
       regions[ feature.iv ] += feature.name
    
    region_names = [feature.name for feature in regionParser]
    
    return regions, region_names

def count_PE_reads(sam_files, labels, regions, file_type="sam", use_chrom_name=False, order="name"):
    """ counts fragments (PE read pairs) for each region from all SAM/BAM files """
    
    assert len(sam_files) == len(labels)
    if use_chrom_name:
        print "INFO: Running in mode for counting per chromosome name."
    
    m = len(sam_files)

    # initialize a list with default zero counts
    all_counts = [collections.Counter( ) for i in range(m)]
    
    # iterate over all sam/bam files
    for j in range(m):
        
        print "INFO: Start to count reads in", sam_files[j], "..."
        
        if file_type == "sam":
            almnt_file = HTSeq.SAM_Reader( sam_files[j] )
        else:
            almnt_file = HTSeq.BAM_Reader( sam_files[j] )
        
        # pair alignment records according to PE pairs and iterate over pairs
        if order == "name":
            print "INFO: Assuming SAM/BAM file ordered by read name."    
            alignmentIterator = HTSeq.pair_SAM_alignments( almnt_file )
        else:
            print "INFO: Assuming SAM/BAM file ordered by position"
            alignmentIterator = HTSeq.pair_SAM_alignments_with_buffer( almnt_file ,max_buffer_size=100*3000000 )
        
        for pair in alignmentIterator: 
            
            first_almnt, second_almnt = pair  # extract pair
            
            # check if both pairs are mapped
            if first_almnt == None or second_almnt == None or not ( first_almnt.aligned and second_almnt.aligned ):
                all_counts[j][ "_unmapped" ] += 1
                continue
            
            
            # potential speed up for transcript fragments as reference
            if use_chrom_name:
                
                if first_almnt.iv.chrom == second_almnt.iv.chrom:
                    all_counts[j][ first_almnt.iv.chrom ] += 1
                else:
                    all_counts[j][ "_no_feature" ] += 1

            else:
                # build set for all regions overalapping with the reads
                gene_ids_first = set()
                gene_ids_second = set()

                # extract all region names that overlap with the reads and add them to set
                for iv, val in regions[ first_almnt.iv ].steps():
                    gene_ids_first |= val
                for iv, val in regions[ second_almnt.iv ].steps():
                    gene_ids_second |= val
            
                # take only those genes that are common for first and second read
                gene_ids = gene_ids_first & gene_ids_second
    
                # handle read-pairs not mapped to a feature
                if len(gene_ids) == 0:
                    all_counts[j][ "_no_feature" ] += 1
                                
                # if pair maps to a unique gene count it
                else:
                    # add increase counter for all genes
                    for gene_id in list(gene_ids):
                        all_counts[j][ gene_id ] += 1
                
    # return counts 
    return(all_counts)

    
def write_counts(all_counts, all_names, labels, outFile):

    m = len(labels)
    
    with open(outFile, 'w') as outHandle:
                
        outHandle.write("\t".join(["region"] + labels) + '\n')
        
        for gene in all_names:
            
            # get counts for the gene from all samples
            geneCounts = [str(all_counts[j][gene]) for j in range(m)]
                                    
            # write tab separated output line
            outHandle.write(gene + "\t" + "\t".join( geneCounts ) + "\n")
    

def write_unmapped_counts(all_counts, labels, outFile):
    
    m = len(labels)
    
    with open(outFile, 'w') as outHandle:
        
        outHandle.write("\t".join(["mapping_type"] + labels) + '\n')
        
        for gene in ["_unmapped", "_no_feature"]:
            
            # get counts for the gene from all samples
            geneCounts = [str(all_counts[j][gene]) for j in range(m)]
                
            # write tab separated output line
            outHandle.write(gene + "\t" + "\t".join( geneCounts ) + "\n")

                

if __name__ == "__main__":
    args = commandline()
    
    if not args.labels:
        from os import path
        args.labels = [path.split(sam)[1] for sam in args.input_files]
        
    assert len(args.input_files) == len(args.labels)
    
    print "INFO: Start parsing of regions in: ", args.bed_file
    if args.use_chrom_name:

        region_names = getRegNamesFromBed(args.bed_file)
        regions=None
    else:
        regions, region_names = getGenomicarrayOfSetsAndNames( args.bed_file )
    print "INFO: Finished parsing of regions."
    
    all_counts = count_PE_reads(args.input_files, args.labels, regions, args.file_type, args.use_chrom_name, args.order)
    
    print "INFO: Finsihed counting of reads."

    print "INFO: Start writing output file:", args.output_file
    write_counts(all_counts, region_names, args.labels, args.output_file)
    
    print "INFO: Start writing of unmapped read counts to file:", args.output_file + ".unmapped_stats"
    write_unmapped_counts(all_counts, args.labels, args.output_file + ".unmapped_stats")

